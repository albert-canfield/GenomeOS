# SPDX-License-Identifier: AGPL-3.0-or-later
"""Where a polygenic score falls among people: a reference distribution from public allele frequencies.

`polygenic.score_person` gives a raw weighted sum, which says nothing about high or low until it is
set against the same sum in a population. The PGS Catalog publishes weights, not that distribution,
so it is modelled here. Under Hardy-Weinberg and independence between the score's variants, a score
S = sum_i w_i d_i with effect-allele frequency p_i in a group has

    mean  = sum_i 2 p_i w_i            variance = sum_i 2 p_i (1 - p_i) w_i^2

and a person's percentile is the normal CDF of (S - mean) / sd. The frequencies are gnomAD v4
genome allele frequencies per genetic-ancestry group (and 1000 Genomes phase 3 superpopulations,
kept as a cross-check), read variant by variant from Ensembl's REST API and cached under
data/knowledge/pgs (git-ignored). The mean needs no independence; the variance does: variants in
linkage disequilibrium add covariance terms the formula leaves out. How much is measured, not
assumed, on the HPRC panel the project already holds (attribution/human_panel.py: 88 haplotype
assemblies of 44 people aligned to hg38): within each ancestry group the variance of the people's
actual scores is set against the formula with the panel's own frequencies (`panel_check`). The
rules that read that comparison are frozen in `RULES` before it is run.

A person is compared only against their own declared group. When the person declares none, or the
group is not one the source covers, the percentile is UNKNOWN; nothing is borrowed from a
neighbouring group. Both sides of the comparison use the same variants: those with a frequency in
the group and a genotype in the person. Variants the person lacks but the group has are not
dropped silently: their expected contribution and spread give the `missing_variants` interval.

A percentile is a position among people of a group on one published score. It is not a risk, not
an absolute probability of the trait, and it inherits the score's own ancestry limits.
"""

from __future__ import annotations

import gzip
import json
import math
import time
import urllib.error
import urllib.request
from bisect import bisect_left
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from statistics import NormalDist
from typing import Any

ENSEMBL = "https://rest.ensembl.org"
NCBI_DATASETS = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{ids}/dataset_report?page_size=100"
HPRC_METADATA = (
    "https://raw.githubusercontent.com/human-pangenomics/HPP_Year1_Data_Freeze_v1.0/main/"
    "sample_metadata/hprc_year1_sample_metadata.txt"
)
KNOWLEDGE = Path("data/knowledge/pgs")
PANEL = Path("data/knowledge/human_panel")
RESULT = Path("data/results/pgs_reference.json")
REGISTRATION = Path("data/results/pgs_reference_registration.json")
USER_AGENT = "GenomeOS/0.1"

GNOMAD = ("afr", "ami", "amr", "asj", "eas", "fin", "mid", "nfe", "sas")
GNOMAD_LABEL = {
    "afr": "African/African American",
    "ami": "Amish",
    "amr": "Admixed American",
    "asj": "Ashkenazi Jewish",
    "eas": "East Asian",
    "fin": "Finnish",
    "mid": "Middle Eastern",
    "nfe": "non-Finnish European",
    "sas": "South Asian",
}
THOUSAND = ("AFR", "AMR", "EAS", "EUR", "SAS")
POPULATIONS = tuple(f"gnomADg:{g}" for g in GNOMAD) + tuple(f"1000GENOMES:phase_3:{g}" for g in THOUSAND)
MAX_FETCH_VARIANTS = 20_000  # per-variant lookups above this are refused: a dense score needs a panel
BATCH = 200  # Ensembl's POST /variation limit

# Ancestry the donors declared, for the people GIAB publishes; an imported person declares theirs with
# an `ancestry` field (a gnomAD group code) in their manifest.json.
DECLARED = {"HG002": "asj", "HG003": "asj", "HG004": "asj"}
DECLARED_SOURCE = (
    "GIAB Ashkenazi Jewish trio (HG002 son NA24385, HG003 father NA24149, HG004 mother NA24143): the "
    "ancestry the donors declared to the Personal Genome Project, not inferred from the genotypes"
)

# Frozen before the panel is read (registered in data/results/pgs_reference_registration.json).
RULES: dict[str, Any] = {
    "coverage_floor": 0.80,
    "coverage_floor_reads": "a percentile needs the variants with both a group frequency and a person "
    "genotype to carry at least this share of the score's total |weight|; below it the percentile is "
    "UNKNOWN",
    "panel_min_people": 10,
    "ld": "R1. Per score, within each HPRC ancestry group with at least panel_min_people people, the "
    "variance of the people's scores (over SNVs placed in every haplotype of the group) is divided by "
    "the independence variance from the same people's allele frequencies (heterozygosity corrected "
    "by 2n/(2n-1)); groups are pooled by degrees of freedom and a 95% interval taken from the "
    "chi-square distribution. Independence HOLDS for the score when that interval contains 1, and "
    "FAILS otherwise. When it holds, a person's percentile interval spans the interval of the ratio. "
    "When it fails, the person's percentile is UNKNOWN and the independence number is kept only as a "
    "diagnostic, because no LD correction measured in the panel's groups is borrowed for another group",
    "calibration": "R2. For each HPRC group with at least panel_min_people people, each person's score "
    "(over SNVs placed in every haplotype of the group and with a gnomAD frequency in the matching "
    "gnomAD group: AFR->afr, AMR->amr) is turned into z = (score - mean) / sd under the gnomAD model. "
    "The model is CONSISTENT with the group when the 95% t interval of the mean z contains 0 and the "
    "95% chi-square interval of the variance of z contains 1, and NOT CONSISTENT otherwise. This is "
    "descriptive: the HPRC groups are admixed Caribbean, West African and Latin American cohorts, "
    "not the gnomAD groups themselves, and it changes no person's percentile",
    "author_frequencies": "R3. Where the scoring file carries the authors' effect-allele frequency, "
    "the gnomAD nfe frequency of each matched variant is compared with it; a variant whose nfe "
    "frequency sits nearer 1 - author than author by more than 0.1 is counted as a possible allele "
    "flip. Descriptive",
}

NOTE = (
    "a percentile is a position among people of the stated group on this one score, modelled from "
    "public allele frequencies under Hardy-Weinberg and independence between variants; it is not a "
    "risk, not an absolute probability of the trait, and it says nothing outside the group it names"
)
EVIDENCE = (
    "curated: gnomAD v4 genome allele frequencies through Ensembl REST; derived: the moments of the "
    "score under Hardy-Weinberg and independence, and the person's position in them"
)
_N = NormalDist()


# -- small statistics, standard library only ------------------------------------------------------
def chi2_quantile(q: float, df: int) -> float:
    """Wilson-Hilferty approximation to the chi-square quantile (good to ~1% for df >= 5)."""
    z = _N.inv_cdf(q)
    a = 2.0 / (9.0 * df)
    return df * max(1.0 - a + z * math.sqrt(a), 0.0) ** 3


def t_quantile(q: float, df: int) -> float:
    """Cornish-Fisher expansion of Student's t quantile about the normal one."""
    z = _N.inv_cdf(q)
    v = float(df)
    return (
        z
        + (z**3 + z) / (4 * v)
        + (5 * z**5 + 16 * z**3 + 3 * z) / (96 * v**2)
        + (3 * z**7 + 19 * z**5 + 17 * z**3 - 15 * z) / (384 * v**3)
    )


def moments(items) -> dict[str, float]:
    """Mean, variance and the sampling error of the mean of sum w*d over (w, p, an) under HWE."""
    mu = var = se2 = 0.0
    n = 0
    for w, p, an in items:
        mu += 2.0 * p * w
        var += 2.0 * p * (1.0 - p) * w * w
        if an:
            se2 += 4.0 * w * w * p * (1.0 - p) / an
        n += 1
    return {"variants": n, "mean": mu, "var": var, "sd": math.sqrt(var), "se_mean": math.sqrt(se2)}


def percentile(z: float) -> float:
    return 100.0 * _N.cdf(z)


# -- frequencies ------------------------------------------------------------------------------------
def ensembl_alleles(effect: str, other: str) -> tuple[str, list[str]]:
    """VCF-style alleles in Ensembl's notation: a shared leading base dropped, an empty allele '-'."""
    others = [o for o in other.split("/") if o] if other else []
    if len(effect) == 1 and all(len(o) == 1 for o in others):
        return effect, others
    out = []
    eff = effect
    for o in others:
        k = 0
        while k < min(len(effect), len(o)) and effect[k] == o[k]:
            k += 1
        eff = effect[k:] or "-"
        out.append(o[k:] or "-")
    return eff, out


def effect_frequency(entries: list[dict], population: str, effect: str) -> tuple[float, int] | None:
    """(effect-allele frequency, allele number) in one population from Ensembl's list, or None.

    An allele Ensembl does not list for a population that it does list is at 1 minus the listed ones.
    """
    rows = [e for e in entries if e.get("population") == population]
    if not rows:
        return None
    an = 0
    for e in rows:
        f, c = e.get("frequency"), e.get("allele_count")
        if f and c:
            an = int(round(float(c) / float(f)))
            break
    hit = [e for e in rows if e.get("allele") == effect]
    if hit:
        p = float(hit[0].get("frequency") or 0.0)
    else:
        seen: dict[str, float] = {}
        for e in rows:
            seen.setdefault(str(e.get("allele")), float(e.get("frequency") or 0.0))
        p = max(0.0, 1.0 - sum(seen.values()))
    return min(max(p, 0.0), 1.0), an


def match_record(row: tuple, record: dict | None, allele: str | None = None) -> tuple[str, dict[str, list]]:
    """(status, {population: [p, an]}) for one scoring row against Ensembl's record of its id.

    `allele` is Ensembl's spelling of the effect allele when the haplotype comparison of
    `indel_allele` has already shown the record to be the same variant; position and spelling are
    then not checked again.
    """
    chrom, pos, effect, other = row[0], row[1], row[2], row[3]
    if record is None:
        return "no_record", {}
    eff, others = ensembl_alleles(effect, other)
    if allele is not None:
        eff = allele
    snv = len(effect) == 1 and all(len(o) == 1 for o in others)
    c = chrom.removeprefix("chr")
    maps = [
        m
        for m in record.get("mappings") or []
        if str(m.get("seq_region_name")) == c and m.get("assembly_name", "GRCh38") == "GRCh38"
    ]
    if not maps:
        return "no_grch38_mapping", {}
    if allele is None:
        ok = [m for m in maps if (m.get("start") == pos if snv else abs(int(m.get("start", 0)) - pos) <= 1)]
        if not ok:
            return "position_mismatch", {}
        alleles = str(ok[0].get("allele_string", "")).split("/")
        if eff not in alleles:
            return "effect_allele_not_at_id", {}
    pops = record.get("populations") or []
    freqs = {}
    for pop in POPULATIONS:
        f = effect_frequency(pops, pop, eff)
        if f is not None:
            freqs[pop] = [round(f[0], 6), f[1]]
    return "matched", freqs


INDEL_WINDOW = 30  # bases either side of an indel within which Ensembl's placement is looked for


def haplotype(window: str, w0: int, start: int, end: int, allele: str) -> str:
    """`window` (1-based, first base at w0) with [start, end] replaced by `allele` ('-' is empty).

    Ensembl writes an insertion with end = start - 1, which this reads as a zero-length slot.
    """
    return window[: start - w0] + ("" if allele == "-" else allele) + window[end - w0 + 1 :]


def indel_allele(
    row: tuple, window: str, w0: int, candidate: dict, ref_allele: str | None = None
) -> str | None:
    """Ensembl's spelling of the row's effect allele when `candidate` (an overlap hit) is the same
    variant, or None.

    Indels in repeats are written at different positions and in different trims by the scoring file
    and by Ensembl, so the test is the sequence each produces: the row's effect and other alleles,
    put on the hg38 window with `ref_allele` as the reference one, must give the same two haplotypes
    as two of the candidate's alleles. Without `ref_allele` every reading the window allows is tried.
    """
    pos, effect = row[1], row[2]
    others = [o for o in row[3].split("/") if o]
    try:
        s, e = int(candidate["start"]), int(candidate["end"])
    except (KeyError, TypeError, ValueError):
        return None
    if not others or s < w0 or e > w0 + len(window) - 1:
        return None
    haps = {a: haplotype(window, w0, s, e, a) for a in candidate.get("alleles") or []}
    if window not in haps.values():
        return None
    for r in [ref_allele] if ref_allele else readings(row, window, w0):

        def put(x: str, r=r) -> str:
            return window[: pos - w0] + x + window[pos - w0 + len(r) :]

        he = put(effect)
        ho = [put(o) for o in others]
        if he in ho:
            continue
        spelled = [a for a, h in haps.items() if h == he]
        if spelled and all(h in haps.values() for h in ho):
            return spelled[0]
    return None


def readings(row: tuple, window: str, w0: int) -> list[str]:
    """The row's alleles that match hg38 at its position: each is a possible reference allele.

    Inside a repeat both can match (C/CA before an A is an insertion of A or a deletion of A), and a
    scoring file without a reference column does not say which.
    """
    pos = row[1]
    alleles = [*(o for o in row[3].split("/") if o), row[2]]
    return [a for a in alleles if window[pos - w0 : pos - w0 + len(a)] == a]


def _reference_windows():
    """fetch(chrom, start1, end1) over the local hg38 sequence, and a closer; None where not fetched."""
    from genomeos.coords import Locus
    from genomeos.genome import IndexedGenome, reference_fasta

    open_: dict[str, Any] = {}

    def fetch(chrom: str, a: int, b: int) -> str | None:
        if chrom not in open_:
            fa = reference_fasta(chrom, Path("data/reference"))
            open_[chrom] = IndexedGenome(str(fa)) if fa.exists() else None
        g = open_[chrom]
        return str(g.fetch(Locus(chrom, a - 1, b))).upper() if g is not None else None

    def close() -> None:
        for g in open_.values():
            if g is not None:
                g.close()

    return fetch, close


def _indel_lookup(row: tuple, window_at) -> tuple[str, str, int]:
    """(rsID, Ensembl's effect allele, bytes) for an indel row placed by haplotype, or ("", "", bytes)."""
    chrom, pos = row[0], row[1]
    longest = max([len(row[2])] + [len(o) for o in row[3].split("/") if o])
    w0, w1 = max(1, pos - INDEL_WINDOW), pos + longest + INDEL_WINDOW
    window = window_at(chrom, w0, w1)
    if not window:
        return "", "", 0
    d, nbytes = _request(
        f"/overlap/region/human/{chrom.removeprefix('chr')}:{w0 + 1}-{w1 - 1}?feature=variation"
    )
    cands = [v for v in d or [] if str(v.get("id", "")).startswith("rs")]
    hits = {}
    for r in readings(row, window, w0):
        found = [(v["id"], a) for v in cands if (a := indel_allele(row, window, w0, v, r)) is not None]
        if found:
            hits[r] = found[0]
    if len(hits) != 1:  # nothing placed, or placed under two readings of which allele is hg38's
        return "", "", nbytes
    rs, a = next(iter(hits.values()))
    return rs, a, nbytes


def _request(path: str, payload: dict | None = None, tries: int = 12) -> tuple[Any, int]:
    url = path if path.startswith("http") else ENSEMBL + path
    data = json.dumps(payload).encode() if payload is not None else None
    headers = {"Accept": "application/json", "User-Agent": USER_AGENT}
    if data is not None:
        headers["Content-Type"] = "application/json"
    for k in range(tries):
        req = urllib.request.Request(url, data=data, headers=headers, method="POST" if data else "GET")
        try:
            with urllib.request.urlopen(req, timeout=300) as r:  # noqa: S310
                body = r.read()
            return json.loads(body), len(body)
        except urllib.error.HTTPError as ex:
            if ex.code in (429, 500, 502, 503, 504) and k < tries - 1:
                time.sleep(min(float(ex.headers.get("Retry-After") or 2**k), 60))
                continue
            if ex.code in (400, 404):
                return None, 0
            raise
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError):
            if k < tries - 1:
                time.sleep(min(2**k, 60))
                continue
            raise
    return None, 0


def _id_at(row: tuple) -> tuple[str, int]:
    """The rsID Ensembl places at a row's position with its alleles, for rows the file gives none."""
    chrom, pos, effect, other = row[0], row[1], row[2], row[3]
    eff, others = ensembl_alleles(effect, other)
    longest = max([len(effect)] + [len(o) for o in other.split("/") if o])
    end = pos + longest - 1
    d, nbytes = _request(f"/overlap/region/human/{chrom.removeprefix('chr')}:{pos}-{end}?feature=variation")
    for v in d or []:
        alleles = v.get("alleles") or []
        if str(v.get("id", "")).startswith("rs") and eff in alleles and all(o in alleles for o in others):
            return v["id"], nbytes
    return "", nbytes


def af_path(pgs_id: str, knowledge: Path = KNOWLEDGE) -> Path:
    return knowledge / f"{pgs_id}_af.json.gz"


def load_frequencies(pgs_id: str, knowledge: Path = KNOWLEDGE) -> dict | None:
    p = af_path(pgs_id, knowledge)
    if not p.exists():
        return None
    with gzip.open(p, "rt") as fh:
        return json.load(fh)


def fetch_frequencies(
    pgs_id: str, knowledge: Path = KNOWLEDGE, progress=None, max_variants: int = MAX_FETCH_VARIANTS
) -> dict | None:
    """Every scoring row's effect-allele frequency per population, fetched once and cached.

    None when the score has more variants than `max_variants`: a genome-wide score of a million
    variants is not looked up one by one (and its variants sit in LD, so the independence variance
    would not hold for it anyway).
    """
    have = load_frequencies(pgs_id, knowledge)
    if have is not None:
        return have
    from genomeos.genome.polygenic import scoring_file, weight_table

    with gzip.open(scoring_file(pgs_id, knowledge), "rt") as fh:
        _, rows = weight_table(fh)
    if len(rows) > max_variants:
        return None
    t0 = time.time()
    fetched = requests = 0
    ids = [r[5] if r[5].startswith("rs") else "" for r in rows]
    # Progress is kept in a partial file, so a run that Ensembl's 503s stop resumes where it ended.
    partial = knowledge / f"{pgs_id}_af_partial.json.gz"
    state: dict[str, dict] = {"ids": {}, "records": {}}
    if partial.exists():
        with gzip.open(partial, "rt") as fh:
            state = json.load(fh)

    def save() -> None:
        knowledge.mkdir(parents=True, exist_ok=True)
        with gzip.open(partial, "wt") as fh:
            json.dump(state, fh)

    # rows the file gives no id are placed by position, a few lookups at a time: Ensembl answered one
    # overlap request in 6 to 25 seconds on 2026-10-10, and its limit is 15 requests a second
    no_id = [i for i, x in enumerate(ids) if not x]
    need = [i for i in no_id if str(i) not in state["ids"]]
    failed = 0
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = {pool.submit(_id_at, rows[i]): i for i in need}
        for n, fut in enumerate(as_completed(futures), 1):
            try:
                found, nb = fut.result()
            except (urllib.error.URLError, TimeoutError, OSError):
                failed += 1
                continue
            state["ids"][str(futures[fut])] = found
            fetched += nb
            requests += 1
            if n % 25 == 0:
                save()
    save()
    if failed:
        raise RuntimeError(f"{pgs_id}: {failed} position lookups failed after retries; run again to resume")
    for i in no_id:
        ids[i] = state["ids"][str(i)]
    by_position = len(no_id)
    if progress and no_id:
        progress(f"{pgs_id}: {sum(1 for i in no_id if ids[i]):,} of {len(no_id):,} position-only rows placed")
    unique = sorted({x for x in ids if x} - set(state["records"]))
    for k in range(0, len(unique), BATCH):
        chunk = unique[k : k + BATCH]
        d, nb = _request("/variation/homo_sapiens?pops=1", {"ids": chunk})
        fetched += nb
        requests += 1
        for x in chunk:
            state["records"][x] = (d or {}).get(x)
        save()
        if progress:
            progress(f"{pgs_id}: {min(k + BATCH, len(unique)):,} of {len(unique):,} ids")
        time.sleep(0.1)
    records = {x: v for x, v in state["records"].items() if v is not None}
    first = []
    for i, r in enumerate(rows):
        status, freqs = match_record(r, records.get(ids[i]) if ids[i] else None)
        first.append(("no_id_at_position" if not ids[i] else status, freqs))
    # indels the first pass could not match are placed by the haplotypes they produce
    indels = [
        i
        for i, r in enumerate(rows)
        if first[i][0] != "matched"
        and not (len(r[2]) == 1 and all(len(o) == 1 for o in r[3].split("/") if o))
    ]
    state.setdefault("indels", {})
    todo = [i for i in indels if str(i) not in state["indels"]]
    window_at, close = _reference_windows()
    try:
        for n, i in enumerate(todo, 1):
            rs, allele, nb = _indel_lookup(rows[i], window_at)
            state["indels"][str(i)] = [rs, allele]
            fetched += nb
            requests += 1
            if n % 10 == 0:
                save()
    finally:
        close()
    save()
    extra = sorted({v[0] for v in state["indels"].values() if v[0]} - set(state["records"]))
    for k in range(0, len(extra), BATCH):
        chunk = extra[k : k + BATCH]
        d, nb = _request("/variation/homo_sapiens?pops=1", {"ids": chunk})
        fetched += nb
        requests += 1
        for x in chunk:
            state["records"][x] = (d or {}).get(x)
        save()
    records = {x: v for x, v in state["records"].items() if v is not None}
    by_haplotype = 0
    variants = []
    for i, r in enumerate(rows):
        status, freqs = first[i]
        rid = ids[i]
        placed = state["indels"].get(str(i)) if i in indels else None
        if placed and placed[0]:
            s2, f2 = match_record(r, records.get(placed[0]), allele=placed[1])
            if s2 == "matched":
                status, freqs, rid = "matched_by_haplotype", f2, placed[0]
                by_haplotype += 1
        variants.append([r[0], r[1], r[2], r[3], r[4], rid, r[6], status, freqs])
    table = {
        "pgs_id": pgs_id,
        "source": {
            "service": ENSEMBL,
            "endpoints": [
                "POST /variation/homo_sapiens?pops=1",
                "GET /overlap/region/human (rows with no id)",
            ],
            "populations": list(POPULATIONS),
            "gnomad": "gnomAD v4 genomes as Ensembl serves them (group codes with `remaining`), CC0 1.0",
            "thousand_genomes": "1000 Genomes phase 3 superpopulations as Ensembl serves them",
        },
        "date": time.strftime("%Y-%m-%d"),
        "requests": requests,
        "bytes": fetched,
        "ids_found_by_position": by_position,
        "indels_matched_by_haplotype": by_haplotype,
        "seconds": round(time.time() - t0, 1),
        "columns": ["chrom", "pos", "effect", "other", "weight", "id", "author_af", "status", "freqs"],
        "variants": variants,
    }
    knowledge.mkdir(parents=True, exist_ok=True)
    with gzip.open(af_path(pgs_id, knowledge), "wt") as fh:
        json.dump(table, fh)
    partial.unlink(missing_ok=True)
    return table


def aligned(table: dict, rows: list[tuple]) -> bool:
    """The cached table describes exactly these scoring rows, in order."""
    v = table.get("variants") or []
    return len(v) == len(rows) and all(
        (a[0], a[1], a[2]) == (b[0], b[1], b[2]) for a, b in zip(v, rows, strict=True)
    )


# -- a population's distribution --------------------------------------------------------------------
def group_distribution(table: dict, population: str) -> dict[str, Any]:
    """The full-score reference in one population: coverage, mean and sd over variants with a frequency."""
    v = table["variants"]
    total_abs = sum(abs(r[4]) for r in v) or 1.0
    total_sq = sum(r[4] * r[4] for r in v) or 1.0
    have = [(r[4], r[8][population][0], r[8][population][1]) for r in v if population in r[8]]
    m = moments(have)
    return {
        "variants_with_frequency": len(have),
        "variants_in_score": len(v),
        "coverage": round(len(have) / len(v), 4) if v else None,
        "weight_share_with_frequency": round(sum(abs(w) for w, _, _ in have) / total_abs, 4),
        "weight_sq_share_without_frequency": round(1.0 - sum(w * w for w, _, _ in have) / total_sq, 4),
        "mean": round(m["mean"], 6),
        "sd": round(m["sd"], 6),
        "se_mean": round(m["se_mean"], 6),
        "median_allele_number": sorted(an for _, _, an in have)[len(have) // 2] if have else None,
    }


def author_check(table: dict, population: str = "gnomADg:nfe") -> dict[str, Any] | None:
    """R3: the authors' effect-allele frequencies against gnomAD's for the same alleles."""
    pairs = [
        (r[6], r[8][population][0], r[2], r[3])
        for r in table["variants"]
        if r[6] is not None and population in r[8]
    ]
    if not pairs:
        return None
    n = len(pairs)
    ma = sum(a for a, _, _, _ in pairs) / n
    mb = sum(b for _, b, _, _ in pairs) / n
    sab = sum((a - ma) * (b - mb) for a, b, _, _ in pairs)
    saa = sum((a - ma) ** 2 for a, _, _, _ in pairs)
    sbb = sum((b - mb) ** 2 for _, b, _, _ in pairs)
    flips = [(a, b, e, o) for a, b, e, o in pairs if abs(b - (1 - a)) + 0.1 < abs(b - a)]
    pal = {frozenset("AT"), frozenset("CG")}
    return {
        "population": population,
        "variants_compared": n,
        "variants_with_author_frequency": sum(1 for r in table["variants"] if r[6] is not None),
        "pearson_r": round(sab / math.sqrt(saa * sbb), 4) if saa and sbb else None,
        "mean_abs_difference": round(sum(abs(a - b) for a, b, _, _ in pairs) / n, 4),
        "possible_flips": len(flips),
        "possible_flips_palindromic": sum(1 for _, _, e, o in flips if frozenset(e + o) in pal),
    }


# -- the person ---------------------------------------------------------------------------------------
def declared_ancestry(name: str, root: Path | None = None) -> dict[str, Any]:
    from genomeos.genome.individuals import ROOT

    root = root or ROOT
    m = root / name / "manifest.json"
    declared = None
    source = None
    if m.exists():
        try:
            declared = json.loads(m.read_text()).get("ancestry")
        except (OSError, json.JSONDecodeError):
            declared = None
        if declared:
            source = "the person's manifest.json, field `ancestry`"
    if not declared and name in DECLARED:
        declared, source = DECLARED[name], DECLARED_SOURCE
    if not declared:
        return {
            "declared": None,
            "group": None,
            "source": "none declared: an `ancestry` field in the person's manifest.json, set to one of "
            + ", ".join(GNOMAD)
            + ", names the group",
        }
    g = str(declared).strip().lower()
    return {
        "declared": declared,
        "group": f"gnomADg:{g}" if g in GNOMAD else None,
        "label": GNOMAD_LABEL.get(g),
        "source": source,
    }


def ld_reading(pgs_id: str, result: Path = RESULT) -> dict[str, Any] | None:
    """R1's pooled reading for a score from the committed panel result, or None when not measured."""
    if not result.exists():
        return None
    try:
        d = json.loads(result.read_text())
    except (OSError, json.JSONDecodeError):
        return None
    s = (d.get("scores") or {}).get(pgs_id) or {}
    return ((s.get("panel") or {}).get("ld") or {}).get("pooled")


def person_reference(
    name: str,
    pgs_id: str,
    rows: list[tuple],
    dosage: list[int | None],
    root: Path | None = None,
    knowledge: Path = KNOWLEDGE,
    fetch: bool = False,
    progress=None,
    result: Path = RESULT,
) -> dict[str, Any]:
    """The person's percentile in their declared group's modelled distribution, or UNKNOWN and why."""
    anc = declared_ancestry(name, root)
    base: dict[str, Any] = {
        "ancestry": anc,
        "rules": {"coverage_floor": RULES["coverage_floor"]},
        "note": NOTE,
    }

    def unknown(reason: str, **extra) -> dict[str, Any]:
        return {"status": "UNKNOWN", "percentile": None, "reason": reason, **base, **extra}

    if anc["group"] is None:
        if anc["declared"]:
            return unknown(f"the declared ancestry {anc['declared']!r} is not a group the source covers")
        return unknown("no ancestry declared for this person, and no group is borrowed")
    table = load_frequencies(pgs_id, knowledge)
    if table is None and fetch:
        try:
            table = fetch_frequencies(pgs_id, knowledge, progress)
        except (urllib.error.URLError, TimeoutError, OSError) as ex:
            return unknown(f"the allele frequencies could not be fetched: {ex}")
    if table is None:
        if len(rows) > MAX_FETCH_VARIANTS:
            return unknown(
                f"{len(rows):,} variants: a dense genome-wide score is not looked up variant by variant, "
                "and its variants sit in linkage disequilibrium, so the independence variance does not "
                "hold for it; a genotyped reference panel of the person's group is needed"
            )
        return unknown("no allele-frequency table cached for this score (`scripts/pgs_reference.py --fetch`)")
    if not aligned(table, rows):
        return unknown("the cached frequency table does not match this scoring file row for row")
    pop = anc["group"]
    total_abs = sum(abs(r[4]) for r in rows) or 1.0
    comp: list[tuple[float, float, int]] = []
    miss: list[tuple[float, float, int]] = []
    x = 0.0
    no_af = no_af_but_genotyped = 0
    no_af_abs = 0.0
    for r, d in zip(table["variants"], dosage, strict=True):
        f = r[8].get(pop)
        if f is None:
            no_af += 1
            no_af_abs += abs(r[4])
            no_af_but_genotyped += d is not None
            continue
        if d is None:
            miss.append((r[4], f[0], f[1]))
            continue
        comp.append((r[4], f[0], f[1]))
        x += r[4] * d
    share = sum(abs(w) for w, _, _ in comp) / total_abs
    counts = {
        "in_score": len(rows),
        "with_frequency_in_group": len(rows) - no_af,
        "comparable": len(comp),
        "genotyped_without_frequency": no_af_but_genotyped,
        "frequency_without_genotype": len(miss),
        "without_frequency": no_af,
        "weight_share_comparable": round(share, 4),
        "weight_share_without_frequency": round(no_af_abs / total_abs, 4),
    }
    m = moments(comp)
    if not m["sd"]:
        return unknown("the comparable variants carry no variance in this group", variants=counts)
    mm = moments(miss)
    z = (x - m["mean"]) / m["sd"]
    ld = ld_reading(pgs_id, result)
    diag = {
        "group": pop,
        "person_partial_score": round(x, 6),
        "mean": round(m["mean"], 6),
        "sd": round(m["sd"], 6),
        "z": round(z, 4),
        "percentile_under_independence": round(percentile(z), 2),
    }
    if share < RULES["coverage_floor"]:
        return unknown(
            f"the comparable variants carry {share:.1%} of the score's |weight|, below the floor of "
            f"{RULES['coverage_floor']:.0%}",
            variants=counts,
            diagnostic=diag,
        )
    if ld is not None and ld.get("reading") == "fails":
        return unknown(
            "R1: on the HPRC panel the score's variance departs from the independence model "
            f"(ratio {ld.get('ratio')}, 95% interval {ld.get('ci')}), and no LD correction is borrowed",
            variants=counts,
            diagnostic=diag,
            ld=ld,
        )
    r_lo, r_hi = (ld["ci"][0], ld["ci"][1]) if ld is not None else (1.0, 1.0)
    num = x - m["mean"]

    def pct(numer: float, var: float) -> float:
        return percentile(numer / math.sqrt(var))

    af_iv = sorted(pct(num + s * 1.96 * m["se_mean"], m["var"]) for s in (-1, 1))
    ld_iv = sorted(pct(num, m["var"] * r) for r in (r_lo, r_hi)) if ld is not None else None
    full = m["var"] + mm["var"]
    miss_iv = sorted(pct(num + s * 1.96 * math.sqrt(mm["var"]), full) for s in (-1, 1))
    spread = 1.96 * math.sqrt(m["se_mean"] ** 2 + mm["var"])
    corners = [pct(num + s * spread, m["var"] * r + mm["var"]) for s in (-1, 1) for r in (r_lo, r_hi)]
    return {
        "status": "ok",
        **base,
        "group": pop,
        "group_label": anc.get("label"),
        "variants": counts,
        "person_partial_score": round(x, 6),
        "mean": round(m["mean"], 6),
        "sd": round(m["sd"], 6),
        "z": round(z, 4),
        "percentile": round(percentile(z), 2),
        "interval": {
            "af_sampling": [round(v, 2) for v in af_iv],
            "ld": [round(v, 2) for v in ld_iv] if ld_iv else None,
            "missing_variants": [round(v, 2) for v in miss_iv],
            "combined": [round(min(corners), 2), round(max(corners), 2)],
        },
        "ld": ld,
        "ld_note": None
        if ld is not None
        else "LD not measured on the panel for this score: the interval omits the variance-model term",
        "missing_variants_mean_contribution": round(mm["mean"], 6),
        "evidence": EVIDENCE,
    }


# -- the HPRC panel: what LD and the model do in real genomes --------------------------------------
def hprc_people(assemblies: list[str], knowledge: Path = KNOWLEDGE) -> dict[str, dict[str, str]]:
    """Accession (without version) -> {sample, superpopulation, subpopulation} for the panel's assemblies."""
    cache = knowledge / "hprc_samples.json"
    if cache.exists():
        return json.loads(cache.read_text())
    acc = sorted({a.split(".")[0] for a in assemblies if a.startswith("GCA_")})
    d, _ = _request(NCBI_DATASETS.format(ids=",".join(acc)))
    sample_of = {
        r["accession"].split(".")[0]: (r.get("organism") or {}).get("infraspecific_names", {}).get("isolate")
        for r in (d or {}).get("reports", [])
    }
    req = urllib.request.Request(HPRC_METADATA, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=120) as r:  # noqa: S310
        text = r.read().decode()
    pops: dict[str, tuple[str, str]] = {}
    for line in text.splitlines()[1:]:
        f = line.split("\t")
        if len(f) >= 6:
            sub, sup = f[4].strip(), f[5].strip()
            if not sup and "MAASAI" in line.upper():
                sub, sup = "MKK", "AFR"  # HapMap Maasai in Kinyawa, Kenya: an African cohort
            pops.setdefault(f[0], (sub, sup))
    out = {
        a: {
            "sample": s,
            "subpopulation": pops.get(s, ("", ""))[0],
            "superpopulation": pops.get(s, ("", ""))[1],
        }
        for a, s in sample_of.items()
        if s
    }
    knowledge.mkdir(parents=True, exist_ok=True)
    cache.write_text(json.dumps(out, indent=1, sort_keys=True))
    return out


def _haplotype_pairs(assemblies: list[str], people: dict[str, dict[str, str]]) -> dict[str, list[int]]:
    pairs: dict[str, list[int]] = {}
    for k, a in enumerate(assemblies):
        p = people.get(a.split(".")[0])
        if p:
            pairs.setdefault(p["sample"], []).append(k)
    return {s: ix for s, ix in pairs.items() if len(ix) == 2}


def panel_dosages(
    rows: list[tuple], panel: Path = PANEL, people=None, knowledge: Path = KNOWLEDGE, progress=None
):
    """Yield (row index, {sample: effect dosage or None}) for each SNV row, chromosome by chromosome.

    A haplotype's allele at an hg38 position is the site's recorded base, hg38's base where the site
    says `=` or where no site is recorded inside a block the assembly is present in, and unplaced
    (None) where the assembly is absent, deleted or gapped there. A person's dosage needs both
    haplotypes placed.
    """
    from genomeos.genome.individuals import _reference_base_reader

    by_chrom: dict[str, dict[int, list[int]]] = {}
    for i, r in enumerate(rows):
        if len(r[2]) == 1 and len(r[3]) == 1 and r[2] in "ACGT" and r[3] in "ACGT":
            by_chrom.setdefault(r[0], {}).setdefault(r[1] - 1, []).append(i)
    for chrom, wanted in by_chrom.items():
        d = panel / chrom
        if not (d / "meta.json").exists():
            continue
        assemblies = json.loads((d / "meta.json").read_text())["assemblies"]
        if people is None:
            people = hprc_people(assemblies, knowledge)
        pairs = _haplotype_pairs(assemblies, people)
        sites: dict[int, tuple[str, str]] = {}
        with gzip.open(d / "sites.tsv.gz", "rt") as fh:
            for line in fh:
                p = int(line[: line.index("\t")])
                if p in wanted:
                    _, b, a = line.rstrip("\n").split("\t")
                    sites[p] = (b, a)
        keys = sorted(wanted)
        status: dict[int, str] = {}
        with gzip.open(d / "blocks.tsv.gz", "rt") as fh:
            for line in fh:
                s, e, st, _ = line.rstrip("\n").split("\t")
                s0, e0 = int(s), int(e)
                j = bisect_left(keys, s0)
                while j < len(keys) and keys[j] < e0:
                    status[keys[j]] = st
                    j += 1
        base_at, genome = _reference_base_reader(chrom)
        try:
            for p0 in keys:
                if p0 in sites:
                    ref, a = sites[p0]

                    def allele(k: int, ref=ref, a=a) -> str | None:
                        c = a[k]
                        return ref if c == "=" else (None if c in ".-N" else c)

                elif p0 in status and base_at is not None:
                    ref, st = base_at(p0 + 1), status[p0]

                    def allele(k: int, ref=ref, st=st) -> str | None:
                        return ref if st[k] == "s" else None

                else:

                    def allele(k: int) -> str | None:
                        return None

                for i in wanted[p0]:
                    eff = rows[i][2]
                    out: dict[str, int | None] = {}
                    for sample, (h1, h2) in pairs.items():
                        a1, a2 = allele(h1), allele(h2)
                        out[sample] = None if a1 is None or a2 is None else (a1 == eff) + (a2 == eff)
                    yield i, out
        finally:
            if genome is not None:
                genome.close()
        if progress:
            progress(f"panel {chrom}: {sum(len(v) for v in wanted.values()):,} SNV rows")


def panel_check(
    rows: list[tuple],
    genotypes,
    people: dict[str, dict[str, str]],
    table: dict | None = None,
    min_people: int = RULES["panel_min_people"],
) -> dict[str, Any]:
    """R1 (variance ratio) and R2 (calibration against gnomAD) from the panel's dosages.

    `genotypes` yields (row index, {sample: dosage or None}); nothing per person is returned, only
    group summaries.
    """
    group_of = {p["sample"]: p["superpopulation"] or "unassigned" for p in people.values()}
    members: dict[str, list[str]] = {}
    for s, g in sorted(group_of.items()):
        members.setdefault(g, []).append(s)
    gnomad_for = {"AFR": "gnomADg:afr", "AMR": "gnomADg:amr", "EAS": "gnomADg:eas", "SAS": "gnomADg:sas"}
    acc = {
        g: {
            "n": len(ss),
            "variants": 0,
            "variants_monomorphic": 0,
            "score": dict.fromkeys(ss, 0.0),
            "analytic": 0.0,
            "cal_variants": 0,
            "cal_score": dict.fromkeys(ss, 0.0),
            "cal_mean": 0.0,
            "cal_var": 0.0,
        }
        for g, ss in members.items()
    }
    snv_rows = 0
    for i, dos in genotypes:
        snv_rows += 1
        w = rows[i][4]
        freqs = table["variants"][i][8] if table is not None else {}
        for g, ss in members.items():
            vals = [dos.get(s) for s in ss]
            if any(v is None for v in vals):
                continue
            a = acc[g]
            a["variants"] += 1
            n2 = 2 * len(ss)
            p = sum(vals) / n2
            if p in (0.0, 1.0):
                a["variants_monomorphic"] += 1
            if n2 > 1:
                a["analytic"] += 2 * p * (1 - p) * w * w * n2 / (n2 - 1)
            for s, v in zip(ss, vals, strict=True):
                a["score"][s] += w * v
            pop = gnomad_for.get(g)
            if pop and pop in freqs:
                q = freqs[pop][0]
                a["cal_variants"] += 1
                a["cal_mean"] += 2 * q * w
                a["cal_var"] += 2 * q * (1 - q) * w * w
                for s, v in zip(ss, vals, strict=True):
                    a["cal_score"][s] += w * v
    groups: dict[str, Any] = {}
    num = den = 0.0
    df = 0
    for g, a in acc.items():
        n = a["n"]
        xs = list(a["score"].values())
        emp = sum((x - sum(xs) / n) ** 2 for x in xs) / (n - 1) if n > 1 else None
        ratio = emp / a["analytic"] if emp is not None and a["analytic"] else None
        out: dict[str, Any] = {
            "people": n,
            "snv_rows_placed_in_every_haplotype": a["variants"],
            "monomorphic_in_group": a["variants_monomorphic"],
            "empirical_var": round(emp, 8) if emp is not None else None,
            "independence_var": round(a["analytic"], 8),
            "ratio": round(ratio, 4) if ratio is not None else None,
            "pooled": n >= min_people and ratio is not None,
        }
        if n >= min_people and ratio is not None:
            num += (n - 1) * emp
            den += (n - 1) * a["analytic"]
            df += n - 1
        pop = gnomad_for.get(g)
        if pop and a["cal_variants"] and a["cal_var"] > 0 and n > 1:
            zs = [(x - a["cal_mean"]) / math.sqrt(a["cal_var"]) for x in a["cal_score"].values()]
            mz = sum(zs) / n
            vz = sum((z - mz) ** 2 for z in zs) / (n - 1)
            cal: dict[str, Any] = {
                "gnomad_group": pop,
                "variants": a["cal_variants"],
                "mean_z": round(mz, 4),
                "var_z": round(vz, 4),
            }
            if n >= min_people:
                half = t_quantile(0.975, n - 1) * math.sqrt(vz / n)
                vlo = vz * (n - 1) / chi2_quantile(0.975, n - 1)
                vhi = vz * (n - 1) / chi2_quantile(0.025, n - 1)
                cal["mean_z_ci"] = [round(mz - half, 4), round(mz + half, 4)]
                cal["var_z_ci"] = [round(vlo, 4), round(vhi, 4)]
                cal["reading"] = (
                    "consistent" if (mz - half) <= 0 <= (mz + half) and vlo <= 1 <= vhi else "not consistent"
                )
            else:
                cal["reading"] = f"not judged: fewer than {min_people} people"
            out["calibration"] = cal
        groups[g] = out
    pooled: dict[str, Any] | None = None
    if df and den:
        r = num / den
        lo = r * df / chi2_quantile(0.975, df)
        hi = r * df / chi2_quantile(0.025, df)
        pooled = {
            "groups": [g for g, o in groups.items() if o["pooled"]],
            "df": df,
            "ratio": round(r, 4),
            "ci": [round(lo, 4), round(hi, 4)],
            "reading": "holds" if lo <= 1 <= hi else "fails",
        }
    return {"snv_rows": snv_rows, "groups": groups, "pooled": pooled}


# -- words for the CLI and the dossier ------------------------------------------------------------------
def describe(reference: dict[str, Any] | None) -> str:
    """One line: the percentile with its interval, group and denominators, or UNKNOWN and why."""
    if not reference:
        return "percentile UNKNOWN: no reference computed"
    if reference.get("status") != "ok":
        return f"percentile UNKNOWN: {reference.get('reason')}"
    v = reference["variants"]
    lo, hi = reference["interval"]["combined"]
    ld = "" if reference.get("ld") else "; LD not measured"
    return (
        f"percentile {reference['percentile']:.1f} (interval {lo:.1f} to {hi:.1f}{ld}) among gnomAD "
        f"{reference.get('group_label') or reference['group']} genomes, the declared group, over "
        f"{v['comparable']:,} of {v['in_score']:,} variants ({v['weight_share_comparable']:.0%} of the "
        "|weight|); a position among people, not a risk"
    )


def cell(reference: dict[str, Any] | None) -> str:
    """The dossier's short form: `62.1 (48.0-74.9), asj` or `UNKNOWN`."""
    if not reference or reference.get("status") != "ok":
        return "UNKNOWN"
    lo, hi = reference["interval"]["combined"]
    return f"{reference['percentile']:.1f} ({lo:.1f}-{hi:.1f}), {reference['group'].split(':')[-1]}"
