# SPDX-License-Identifier: AGPL-3.0-or-later
"""CPTAC phosphosites on cBioPortal, mapped to UniProt residues and joined to the curated sites.

Two stages, each cached under data/cache/cptac/ (git-ignored):

    uv run python scripts/cptac_phospho.py map    # entities -> single sites -> UniProt accession + residue
    uv run python scripts/cptac_phospho.py join   # curated join, values fetched, per-site summary written

``map`` reads the generic-assay entities of the phospho profiles, keeps single localised sites, maps
each RefSeq protein to a reviewed UniProt accession (UniProt REST ID mapping), fetches the RefSeq
sequences (NCBI E-utilities) and carries positions over (identical sequence, or a unique 15-residue
window), requiring the residue letter to match at both ends. Every drop is counted with its reason.
``join`` is scored against the registration in docs/PROTEIN.md, "Relative abundance in tumours
(2026-09-28)", and the ``CPTAC_*`` constants in genomeos/molecules/ptm.py. The committed result holds
per-site, per-study summaries only (tumours with a value, median log2 ratio); no per-patient value is
written anywhere outside the cache, and the cached value downloads are deleted after the summary.
"""

from __future__ import annotations

import gzip
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from genomeos.cancer.cbioportal import CBioPortal  # noqa: E402
from genomeos.molecules import ptm  # noqa: E402

CACHE = Path("data/cache/cptac")
UA = {"User-Agent": "GenomeOS/0.1"}
DOWNLOAD_CAP = 2 * 1024**3

PROFILES = (
    # RefSeq-keyed: the accession is in the entity
    "luad_cptac_2020_phosphoproteome",
    "gbm_cptac_2021_phosphoproteome",
    "paad_cptac_2021_phosphoproteome",
    "brain_cptac_2020_phosphoprotein",
    # gene-keyed: gene symbol and a position on an unnamed RefSeq isoform
    "brca_cptac_2020_phosphoproteome",
    "ucec_cptac_2020_phosphoproteome",
    "coad_cptac_2019_phosphoprotein_quantification",
    "brca_tcga_phosphoprotein_quantification",
    "ov_tcga_phosphoprotein_quantification",
    "brca_tcga_pan_can_atlas_2018_phosphoprotein_quantification",
    "ov_tcga_pan_can_atlas_2018_phosphoprotein_quantification",
    # gene-level aggregate named "acetylprotein": excluded at parse time
    "lusc_cptac_2021_phosphoproteome",
)

_downloaded = 0


def _count(n: int) -> None:
    global _downloaded
    _downloaded += n
    if _downloaded > DOWNLOAD_CAP:
        raise SystemExit(f"download cap reached: {_downloaded} bytes")


def _http(url: str, data: bytes | None = None, headers: dict | None = None, tries: int = 4) -> bytes:
    for i in range(tries):
        try:
            req = urllib.request.Request(url, data=data, headers={**UA, **(headers or {})})
            with urllib.request.urlopen(req, timeout=300) as r:  # noqa: S310
                raw = r.read()
                _count(len(raw))
                while raw[:2] == b"\x1f\x8b":
                    raw = gzip.decompress(raw)
                return raw
        except Exception:  # noqa: BLE001
            if i == tries - 1:
                raise
            time.sleep(5 * (i + 1))
    raise RuntimeError("unreachable")


def meta(client: CBioPortal, profile: str) -> list[dict]:
    p = CACHE / f"meta_{profile}.json"
    if not p.exists():
        url = f"{client.base}/generic-assay-meta/{profile}?projection=DETAILED"
        p.write_bytes(_http(url, headers={"Accept-Encoding": "gzip"}))
    return json.loads(p.read_text())


def definitions() -> dict[str, dict[str, Any]]:
    """accession -> {gene, sequence, sites: [(position, class)]} from the compiled definitions."""
    out = {}
    for p in sorted(ptm.CACHE.glob("*.json")):
        if p.name.startswith("_"):
            continue
        try:
            d = json.loads(p.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        if "sections" not in d:
            continue
        acc = ptm.accession_of(d)
        if not acc:
            continue
        ident = d["sections"].get("identity", {}).get("items", {})
        out[acc] = {
            "gene": ident.get("gene") or d.get("gene"),
            "sequence": ptm.sequence_of(d),
            "sites": [(s["start"], s["class"]) for s in ptm.sites(d) if isinstance(s.get("start"), int)],
        }
    return out


def _paged(url: str) -> str:
    """All pages of a UniProt REST result (the stream endpoint fails on large ID-mapping jobs)."""
    lines: list[str] = []
    while url:
        for i in range(4):
            try:
                req = urllib.request.Request(url, headers={**UA, "Accept-Encoding": "gzip"})
                with urllib.request.urlopen(req, timeout=300) as r:  # noqa: S310
                    raw = r.read()
                    link = r.headers.get("Link") or ""
                break
            except Exception:  # noqa: BLE001
                if i == 3:
                    raise
                time.sleep(5 * (i + 1))
        _count(len(raw))
        while raw[:2] == b"\x1f\x8b":
            raw = gzip.decompress(raw)
        page = raw.decode().splitlines()
        lines.extend(page if not lines else page[1:])
        m = re.search(r'<([^>]+)>;\s*rel="next"', link)
        url = m.group(1) if m else ""
    return "\n".join(lines)


def curated_file_counts() -> tuple[Counter, int]:
    """Definition files per accession, and curated phospho sites counted per file (the 41,661 basis)."""
    per_acc: Counter = Counter()
    n = 0
    for p in sorted(ptm.CACHE.glob("*.json")):
        if p.name.startswith("_"):
            continue
        try:
            d = json.loads(p.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        if "sections" not in d:
            continue
        per_acc[ptm.accession_of(d)] += 1
        n += sum(1 for s in ptm.sites(d) if s["class"] == "phospho")
    return per_acc, n


def uniprot_idmapping(ids: list[str]) -> dict[str, list[list[Any]]]:
    p = CACHE / "idmapping.json"
    cached = json.loads(p.read_text()) if p.exists() else {}
    todo = [i for i in ids if i not in cached]
    for k in range(0, len(todo), 500):
        chunk = todo[k : k + 500]
        body = urllib.parse.urlencode(
            {"from": "RefSeq_Protein", "to": "UniProtKB", "ids": ",".join(chunk)}
        ).encode()
        job = json.loads(_http("https://rest.uniprot.org/idmapping/run", data=body))["jobId"]
        while True:
            try:
                st = json.loads(_http(f"https://rest.uniprot.org/idmapping/status/{job}", tries=1))
            except Exception:  # noqa: BLE001  (a finished job redirects to its results, which can fail)
                break
            if st.get("jobStatus") in (None, "FINISHED") or "results" in st:
                break
            time.sleep(3)
        tsv = _paged(
            f"https://rest.uniprot.org/idmapping/uniprotkb/results/{job}?format=tsv&fields=accession,reviewed&size=500"
        )
        for i in chunk:
            cached.setdefault(i, [])
        for line in tsv.splitlines()[1:]:
            f = line.split("\t")
            if len(f) >= 3 and f[0] in cached:
                cached[f[0]].append([f[1], f[2] == "reviewed"])
        p.write_text(json.dumps(cached))
        print(f"  idmapping {k + len(chunk)}/{len(todo)}", flush=True)
    return cached


def refseq_sequences(ids: list[str]) -> dict[str, str]:
    p = CACHE / "refseq.json"
    cached = json.loads(p.read_text()) if p.exists() else {}
    todo = [i for i in ids if i not in cached]
    for k in range(0, len(todo), 200):
        chunk = todo[k : k + 200]
        body = urllib.parse.urlencode(
            {"db": "protein", "rettype": "fasta", "retmode": "text", "id": ",".join(chunk)}
        ).encode()
        fasta = _http("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", data=body).decode()
        got: dict[str, list[str]] = {}
        cur = None
        for line in fasta.splitlines():
            if line.startswith(">"):
                cur = line[1:].split()[0]
                got[cur] = []
            elif cur:
                got[cur].append(line.strip())
        for i in chunk:
            cached[i] = "".join(got.get(i, []))
        if k % 2000 == 0:
            p.write_text(json.dumps(cached))
            print(f"  refseq {k + len(chunk)}/{len(todo)}", flush=True)
        time.sleep(0.4)
    p.write_text(json.dumps(cached))
    return cached


def map_sites() -> dict[str, Any]:
    client = CBioPortal()
    defs = definitions()
    by_gene: dict[str, list[str]] = defaultdict(list)
    for acc, d in defs.items():
        if d["gene"]:
            by_gene[d["gene"]].append(acc)
    parsed: dict[str, list[tuple[str, dict]]] = {}
    drops: dict[str, Counter] = {}
    for prof in PROFILES:
        c: Counter = Counter()
        rows = []
        for e in meta(client, prof):
            r = ptm.parse_cptac_entity(prof, e)
            if "drop" in r:
                c["parse: " + r["drop"]] += 1
            else:
                rows.append((e["stableId"], r))
        parsed[prof] = rows
        drops[prof] = c
    refseqs = sorted({r["refseq"] for rows in parsed.values() for _, r in rows if r["refseq"]})
    print(f"{len(refseqs)} RefSeq accessions", flush=True)
    idmap = uniprot_idmapping(refseqs)
    # A superseded RefSeq version (NP_x.1 when NP_x.2 is current) has no UniProt cross-reference; the
    # unversioned accession names the same protein record, and the sequence of the exact version used by
    # CPTAC is still what is compared below, so this widens the mapping without loosening the check.
    unversioned = sorted({rs.split(".")[0] for rs in refseqs if not idmap.get(rs)})
    by_base = uniprot_idmapping(unversioned) if unversioned else {}
    via_base = set()
    for rs in refseqs:
        if not idmap.get(rs) and by_base.get(rs.split(".")[0]):
            idmap[rs] = by_base[rs.split(".")[0]]
            via_base.add(rs)
    seqs = refseq_sequences(refseqs)

    how: Counter = Counter()
    refseq_status: Counter = Counter()
    ref_to_acc: dict[str, str | None] = {}
    for rs in refseqs:
        hits = sorted({a for a, rev in idmap.get(rs, []) if rev and a in defs})
        ref_to_acc[rs] = None
        if not idmap.get(rs):
            refseq_status["no UniProt entry"] += 1
        elif not hits:
            refseq_status["no reviewed UniProt entry among the compiled definitions"] += 1
        elif len(hits) > 1:
            refseq_status["several compiled UniProt entries"] += 1
        elif not seqs.get(rs):
            refseq_status["RefSeq sequence not retrieved"] += 1
        else:
            ref_to_acc[rs] = hits[0]
            same = seqs[rs] == defs[hits[0]]["sequence"]
            refseq_status["identical to UniProt canonical" if same else "differs from UniProt canonical"] += 1
            if rs in via_base:
                refseq_status["of which mapped through the unversioned accession"] += 1

    mapped: dict[str, list[list[Any]]] = {}
    for prof, rows in parsed.items():
        c = drops[prof]
        out = []
        for sid, r in rows:
            if r["refseq"]:
                acc = ref_to_acc.get(r["refseq"])
                if acc is None:
                    c["map: RefSeq protein not mapped to one compiled UniProt entry"] += 1
                    continue
                src = seqs[r["refseq"]]
                if not 0 < r["position"] <= len(src) or src[r["position"] - 1] != r["residue"]:
                    c["map: residue letter disagrees with the RefSeq sequence"] += 1
                    continue
                pos, why = ptm.transfer_position(src, defs[acc]["sequence"], r["position"])
                if pos is None:
                    c["map: " + why] += 1
                    continue
                tier = "refseq"
            else:
                accs = by_gene.get(r["gene"], [])
                if len(accs) != 1:
                    c["map: gene symbol not one compiled UniProt entry"] += 1
                    continue
                acc, pos, why = accs[0], r["position"], "gene symbol, position as given"
                tier = "gene"
            seq = defs[acc]["sequence"]
            if not 0 < pos <= len(seq) or seq[pos - 1] != r["residue"]:
                c[f"map: residue letter disagrees with UniProt ({tier}-keyed)"] += 1
                continue
            how[f"{tier}: {why}"] += 1
            out.append([sid, acc, pos, r["residue"], tier])
        n = Counter((x[1], x[2]) for x in out)
        keep = [x for x in out if n[(x[1], x[2])] == 1]
        if len(keep) < len(out):
            c["map: several entities on one UniProt site in the profile"] += len(out) - len(keep)
        c["mapped"] = len(keep)
        mapped[prof] = keep
    summary = {
        "entities": {p: len(meta(client, p)) for p in PROFILES},
        "refseq_accessions": len(refseqs),
        "refseq_status": dict(refseq_status.most_common()),
        "how_positions_were_carried": dict(how.most_common()),
        "per_profile": {p: dict(drops[p].most_common()) for p in PROFILES},
        "downloaded_bytes_this_run": _downloaded,
    }
    (CACHE / "mapped.json").write_text(json.dumps({"summary": summary, "mapped": mapped}))
    return summary


def fetch_summaries(client: CBioPortal, profile: str, entities: list[str]) -> dict[str, tuple[int, float]]:
    """entity -> (tumours with a value, median log2 ratio). Per-tumour values stay in memory only."""
    study = client._get(f"/molecular-profiles/{profile}")["studyId"]
    per: dict[str, list[float]] = defaultdict(list)
    for k in range(0, len(entities), 200):
        body = json.dumps({"genericAssayStableIds": entities[k : k + 200], "sampleListId": f"{study}_all"})
        raw = _http(
            f"{client.base}/generic_assay_data/{profile}/fetch?projection=SUMMARY",
            data=body.encode(),
            headers={
                "Content-Type": "application/json",
                "Accept": "application/json",
                "Accept-Encoding": "gzip",
            },
        )
        for rec in json.loads(raw):
            try:
                v = float(rec.get("value"))
            except (TypeError, ValueError):
                continue
            if v == v and abs(v) != float("inf"):
                per[rec["stableId"]].append(v)
        time.sleep(client.sleep)
    return {e: (len(vs), round(median(vs), 3)) for e, vs in per.items() if vs}


def join() -> None:
    # Superseded before its result was committed: join_committed adds the per-profile check that the
    # values are log2 ratios (the pancreatic profile holds log2 intensities). Kept as first registered.
    return join_committed()
    from genomeos import manifest as mf
    from genomeos.results import save_result

    client = CBioPortal()
    m = json.loads((CACHE / "mapped.json").read_text())
    msum, mapped = m["summary"], m["mapped"]
    defs = definitions()
    curated_phospho = {(a, p) for a, d in defs.items() for p, cls in d["sites"] if cls == "phospho"}
    curated_other = {(a, p): cls for a, d in defs.items() for p, cls in d["sites"] if cls != "phospho"}
    n_curated_all = sum(len(d["sites"]) for d in defs.values())

    studies: list[str] = []
    joined: dict[tuple[str, int], dict[str, Any]] = {}
    per_profile: dict[str, dict[str, Any]] = {}
    other_class: Counter = Counter()
    for prof in PROFILES:
        rows = mapped.get(prof, [])
        hit = [x for x in rows if (x[1], x[2]) in curated_phospho]
        for x in rows:
            if (x[1], x[2]) in curated_other:
                other_class[curated_other[(x[1], x[2])]] += 1
        stats = {"mapped_sites": len(rows), "on_curated_phospho_site": len(hit)}
        if hit:
            vals = fetch_summaries(client, prof, [x[0] for x in hit])
            stats["with_a_value"] = len(vals)
            if vals:
                si = len(studies)
                studies.append(prof)
                for sid, acc, pos, res, tier in hit:
                    if sid in vals:
                        n, med = vals[sid]
                        j = joined.setdefault((acc, pos), {"residue": res, "tiers": set(), "rows": []})
                        j["tiers"].add(tier)
                        j["rows"].append([si, n, med])
        per_profile[prof] = stats
        print(prof, stats, f"downloaded {_downloaded / 1e6:.0f} MB", flush=True)

    def frac(keys: set) -> float:
        return round(len(keys) / len(curated_phospho), 4) if curated_phospho else 0.0

    refseq_keys = {k for k, v in joined.items() if "refseq" in v["tiers"]}
    all_keys = set(joined)
    off = [
        f"{a}:{v['residue']}{p}" for (a, p), v in joined.items() if v["residue"] not in ptm.PHOSPHO_ACCEPTORS
    ]

    def mismatch_rate(tier: str) -> float | None:
        bad = good = 0
        for prof, c in msum["per_profile"].items():
            bad += c.get(f"map: residue letter disagrees with UniProt ({tier}-keyed)", 0)
            if tier == "refseq":
                bad += c.get("map: residue letter disagrees with the RefSeq sequence", 0)
            good += sum(1 for x in mapped.get(prof, []) if x[4] == tier)
        return round(bad / (bad + good), 4) if bad + good else None

    reg = ptm.CPTAC_REGISTRATION
    mm_ref, mm_gene = mismatch_rate("refseq"), mismatch_rate("gene")
    gene_withheld = mm_gene is None or mm_gene >= reg["gene_keyed_mismatch_bound"]
    if gene_withheld:
        keep = refseq_keys
        for k in list(joined):
            if k not in keep:
                del joined[k]
            else:
                joined[k]["rows"] = [r for r in joined[k]["rows"] if studies[r[0]] in PROFILES[:4]]
    akt1 = joined.get(reg["check_reported_only"][:2])
    egfr = joined.get(reg["check_present_lung"][:2])
    egfr_ok = bool(egfr) and any(studies[r[0]].startswith("luad") for r in egfr["rows"])
    npm1_ok = reg["check_present_any"][:2] in joined
    tp53_absent = reg["check_absent"] not in joined
    sites: dict[str, list[list[Any]]] = defaultdict(list)
    for (acc, pos), v in sorted(joined.items()):
        sites[acc].append([pos, v["residue"], sorted(v["rows"])])
    n_studies = [len(v["rows"]) for v in joined.values()]
    meds = [r[2] for v in joined.values() for r in v["rows"]]
    payload = {
        "field": ptm.CPTAC_FIELD,
        "meaning": ptm.CPTAC_MEANING,
        "negatives_first": {
            "curated_phospho_sites_without_a_value": len(curated_phospho) - len(joined),
            "mapping_drops_by_profile": msum["per_profile"],
            "refseq_proteins_not_mapped": {
                k: v for k, v in msum["refseq_status"].items() if not k.startswith(("identical", "differs"))
            },
            "mapped_sites_on_a_curated_site_of_another_class": dict(other_class.most_common()),
            "gene_keyed_profiles_withheld": gene_withheld,
        },
        "curated_sites": n_curated_all,
        "curated_phospho_sites": len(curated_phospho),
        "joined_refseq_keyed": len(refseq_keys),
        "joined_any_profile": len(all_keys),
        "joined_committed": len(joined),
        "fraction_of_curated_phospho_refseq_keyed": frac(refseq_keys),
        "fraction_of_curated_phospho_any_profile": frac(all_keys),
        "proteins_with_joined_sites": len(sites),
        "studies_per_joined_site": {
            "median": median(n_studies) if n_studies else None,
            "one": sum(1 for n in n_studies if n == 1),
            "five_or_more": sum(1 for n in n_studies if n >= 5),
        },
        "median_log2_ratio_across_site_study_pairs": {
            "median": round(median(meds), 3) if meds else None,
            "min": min(meds) if meds else None,
            "max": max(meds) if meds else None,
        },
        "per_profile": per_profile,
        "mapping": {
            k: msum[k]
            for k in ("entities", "refseq_accessions", "refseq_status", "how_positions_were_carried")
        },
        "registration": {
            **{k: (list(v) if isinstance(v, tuple) else v) for k, v in reg.items()},
            "refseq_keyed_within": reg["expected_refseq_keyed"][0]
            <= frac(refseq_keys)
            <= reg["expected_refseq_keyed"][1],
            "any_profile_within": reg["expected_any_profile"][0]
            <= frac(all_keys)
            <= reg["expected_any_profile"][1],
            "egfr_y1092_in_lung": egfr_ok,
            "npm1_s125_present": npm1_ok,
            "akt1_s473_reported": [studies[r[0]] for r in akt1["rows"]] if akt1 else [],
            "median_log2_ratio_within_bound": bool(meds)
            and reg["median_log2_ratio_within"][0] <= median(meds) <= reg["median_log2_ratio_within"][1],
            "tp53_m1_absent": tp53_absent,
            "off_acceptor_joined": off,
            "residue_mismatch_rate_refseq_keyed": mm_ref,
            "residue_mismatch_refseq_keyed_under_bound": mm_ref is not None
            and mm_ref < reg["refseq_keyed_mismatch_bound"],
            "residue_mismatch_rate_gene_keyed": mm_gene,
        },
        "studies": studies,
        "sites_format": (
            "accession -> [[position, residue, [[study index, tumours with a value, median log2 ratio]]]]"
        ),
        "sites": dict(sites),
        "licence": ptm.CPTAC_LICENCE,
        "downloaded_bytes": _downloaded,
    }
    manifest = {
        "sources": [
            {
                "accession": f"cBioPortal:{p}",
                "version": "public API, fetched 2026-09-28",
                "licence": "ODbL-1.0",
            }
            for p in studies
        ]
        + [
            {
                "accession": "UniProt REST idmapping RefSeq_Protein->UniProtKB",
                "version": "2026_03",
                "licence": "CC BY 4.0",
            },
            {
                "accession": "NCBI RefSeq protein (E-utilities efetch, versioned accessions)",
                "version": "as versioned",
                "licence": "NCBI public data",
            },
        ],
        "inputs": [mf.input_entry(str(CACHE / "mapped.json"), None), mf.input_entry(str(ptm.CACHE), None)],
        "assembly": "n/a: protein residue positions on UniProt canonical sequences",
        "coordinates": {"base": 1, "interval": "closed"},
        "code": {},
        "parameters": {"flank": 7, "profiles": list(PROFILES), "download_cap_bytes": DOWNLOAD_CAP},
        "exclusions": [
            "multi-site and unlocalised entities",
            "gene-level lusc aggregate",
            "entities whose residue letter disagrees",
            "several entities on one site in one profile",
        ],
        "partitions": "n/a: descriptive join, no fitted model",
    }
    path = save_result("phosphosite_tumour_abundance", payload, manifest=manifest)
    print(json.dumps({k: v for k, v in payload.items() if k not in ("sites",)}, indent=1)[:12000])
    print(path)


def join_committed() -> None:
    """The join as committed: the registered join plus a per-profile units check."""
    from genomeos import manifest as mf
    from genomeos.results import save_result

    client = CBioPortal()
    m = json.loads((CACHE / "mapped.json").read_text())
    msum, mapped = m["summary"], m["mapped"]
    defs = definitions()
    curated_phospho = {(a, p) for a, d in defs.items() for p, cls in d["sites"] if cls == "phospho"}
    curated_other = {(a, p): cls for a, d in defs.items() for p, cls in d["sites"] if cls != "phospho"}
    n_curated_all = sum(len(d["sites"]) for d in defs.values())
    files_per_acc, n_curated_phospho_files = curated_file_counts()

    reg = ptm.CPTAC_REGISTRATION
    lo, hi = reg["median_log2_ratio_within"]
    per_profile: dict[str, dict[str, Any]] = {}
    other_class: Counter = Counter()
    found: list[tuple[tuple[str, int], str, str, str, int, float]] = []  # key, residue, profile, tier, n, med
    ratio_units: dict[str, bool] = {}
    for prof in PROFILES:
        rows = mapped.get(prof, [])
        hit = [x for x in rows if (x[1], x[2]) in curated_phospho]
        for x in rows:
            if (x[1], x[2]) in curated_other:
                other_class[curated_other[(x[1], x[2])]] += 1
        stats: dict[str, Any] = {"mapped_sites": len(rows), "on_curated_phospho_site": len(hit)}
        if hit:
            cache = CACHE / f"summary_{prof}.json"  # per-site counts and medians only, never per tumour
            if cache.exists():
                vals = {k: tuple(v) for k, v in json.loads(cache.read_text()).items()}
            else:
                vals = fetch_summaries(client, prof, [x[0] for x in hit])
                cache.write_text(json.dumps(vals))
            stats["with_a_value"] = len(vals)
            meds = [v[1] for v in vals.values()]
            if meds:
                pm = round(median(meds), 3)
                stats["median_of_site_medians"] = pm
                # units check: ratios to a pooled reference centre near zero; log2 intensities do not
                ratio_units[prof] = lo <= pm <= hi
                stats["log2_ratio_units"] = ratio_units[prof]
            for sid, acc, pos, res, tier in hit:
                if sid in vals:
                    found.append(((acc, pos), res, prof, tier, int(vals[sid][0]), float(vals[sid][1])))
        per_profile[prof] = stats
        print(prof, stats, f"downloaded {_downloaded / 1e6:.0f} MB", flush=True)

    def frac(keys: set) -> float:
        return round(len(keys) / len(curated_phospho), 4) if curated_phospho else 0.0

    def frac_files(keys: set) -> float:
        return round(sum(files_per_acc[a] for a, _ in keys) / n_curated_phospho_files, 4)

    refseq_keys = {f[0] for f in found if f[3] == "refseq"}
    all_keys = {f[0] for f in found}

    def mismatch_rate(tier: str) -> float | None:
        bad = good = 0
        for prof, c in msum["per_profile"].items():
            bad += c.get(f"map: residue letter disagrees with UniProt ({tier}-keyed)", 0)
            if tier == "refseq":
                bad += c.get("map: residue letter disagrees with the RefSeq sequence", 0)
            good += sum(1 for x in mapped.get(prof, []) if x[4] == tier)
        return round(bad / (bad + good), 4) if bad + good else None

    mm_ref, mm_gene = mismatch_rate("refseq"), mismatch_rate("gene")
    gene_withheld = mm_gene is None or mm_gene >= reg["gene_keyed_mismatch_bound"]
    committed_profiles = [
        p for p in PROFILES if ratio_units.get(p) and (p in PROFILES[:4] or not gene_withheld)
    ]
    studies = committed_profiles
    joined: dict[tuple[str, int], dict[str, Any]] = {}
    for key, res, prof, _tier, n, med in found:
        if prof in committed_profiles:
            j = joined.setdefault(key, {"residue": res, "rows": []})
            j["rows"].append([studies.index(prof), n, round(med, 3)])
    off = [
        f"{a}:{v['residue']}{p}" for (a, p), v in joined.items() if v["residue"] not in ptm.PHOSPHO_ACCEPTORS
    ]
    akt1 = sorted({f[2] for f in found if f[0] == reg["check_reported_only"][:2]})
    egfr_ok = any(f[0] == reg["check_present_lung"][:2] and f[2].startswith("luad") for f in found)
    npm1_ok = reg["check_present_any"][:2] in all_keys
    tp53_absent = reg["check_absent"] not in all_keys
    sites: dict[str, list[list[Any]]] = defaultdict(list)
    for (acc, pos), v in sorted(joined.items()):
        sites[acc].append([pos, v["residue"], sorted(v["rows"])])
    n_studies = [len(v["rows"]) for v in joined.values()]
    meds = [r[2] for v in joined.values() for r in v["rows"]]
    payload = {
        "field": ptm.CPTAC_FIELD,
        "meaning": ptm.CPTAC_MEANING,
        "negatives_first": {
            "curated_phospho_sites_without_a_value": len(curated_phospho) - len(joined),
            "mapping_drops_by_profile": msum["per_profile"],
            "refseq_proteins_not_mapped": {
                k: v
                for k, v in msum["refseq_status"].items()
                if not k.startswith(("identical", "differs", "of which"))
            },
            "profiles_not_in_log2_ratio_units": sorted(p for p, ok in ratio_units.items() if not ok),
            "mapped_sites_on_a_curated_site_of_another_class": dict(other_class.most_common()),
            "gene_keyed_profiles_withheld": gene_withheld,
        },
        "curated_sites": n_curated_all,
        "curated_phospho_sites": len(curated_phospho),
        "curated_phospho_sites_counted_per_definition_file": n_curated_phospho_files,
        "joined_refseq_keyed": len(refseq_keys),
        "joined_any_profile": len(all_keys),
        "joined_committed": len(joined),
        "fraction_of_curated_phospho_refseq_keyed": frac(refseq_keys),
        "fraction_of_curated_phospho_any_profile": frac(all_keys),
        "fraction_of_curated_phospho_committed": frac(set(joined)),
        "fraction_per_definition_file": {
            "refseq_keyed": frac_files(refseq_keys),
            "any_profile": frac_files(all_keys),
            "committed": frac_files(set(joined)),
        },
        "committed_profiles": committed_profiles,
        "proteins_with_joined_sites": len(sites),
        "studies_per_joined_site": {
            "median": median(n_studies) if n_studies else None,
            "one": sum(1 for n in n_studies if n == 1),
            "all_committed_profiles": sum(1 for n in n_studies if n == len(studies)),
        },
        "median_log2_ratio_across_site_study_pairs": {
            "median": round(median(meds), 3) if meds else None,
            "min": min(meds) if meds else None,
            "max": max(meds) if meds else None,
        },
        "per_profile": per_profile,
        "mapping": {
            k: msum[k]
            for k in ("entities", "refseq_accessions", "refseq_status", "how_positions_were_carried")
        },
        "registration": {
            **{k: (list(v) if isinstance(v, tuple) else v) for k, v in reg.items()},
            "refseq_keyed_within": reg["expected_refseq_keyed"][0]
            <= frac_files(refseq_keys)
            <= reg["expected_refseq_keyed"][1],
            "any_profile_within": reg["expected_any_profile"][0]
            <= frac_files(all_keys)
            <= reg["expected_any_profile"][1],
            "egfr_y1092_in_lung": egfr_ok,
            "npm1_s125_present": npm1_ok,
            "akt1_s473_reported_in": akt1,
            "median_log2_ratio_refseq_keyed_as_registered": round(
                median([f[5] for f in found if f[3] == "refseq"]), 3
            ),
            "median_log2_ratio_within_bound_committed": bool(meds)
            and reg["median_log2_ratio_within"][0] <= median(meds) <= reg["median_log2_ratio_within"][1],
            "tp53_m1_absent": tp53_absent,
            "off_acceptor_joined": off,
            "residue_mismatch_rate_refseq_keyed": mm_ref,
            "residue_mismatch_refseq_keyed_under_bound": mm_ref is not None
            and mm_ref < reg["refseq_keyed_mismatch_bound"],
            "residue_mismatch_rate_gene_keyed": mm_gene,
        },
        "studies": studies,
        "sites_format": (
            "accession -> [[position, residue, [[study index, tumours with a value, median log2 ratio]]]]"
        ),
        "sites": dict(sites),
        "licence": ptm.CPTAC_LICENCE,
        "downloaded_bytes": _downloaded,
    }
    manifest = {
        "sources": [
            {
                "accession": f"cBioPortal:{p}",
                "version": "public API, fetched 2026-09-28",
                "licence": "ODbL-1.0",
            }
            for p in studies
        ]
        + [
            {
                "accession": "UniProt REST idmapping RefSeq_Protein->UniProtKB",
                "version": "2026_03",
                "licence": "CC BY 4.0",
            },
            {
                "accession": "NCBI RefSeq protein (E-utilities efetch, versioned accessions)",
                "version": "as versioned",
                "licence": "NCBI public data",
            },
        ],
        "inputs": [mf.input_entry(str(CACHE / "mapped.json"), None), mf.input_entry(str(ptm.CACHE), None)],
        "assembly": "n/a: protein residue positions on UniProt canonical sequences",
        "coordinates": {"base": 1, "interval": "closed"},
        "code": {},
        "parameters": {"flank": 7, "profiles": list(PROFILES), "download_cap_bytes": DOWNLOAD_CAP},
        "exclusions": [
            "multi-site and unlocalised entities",
            "gene-level lusc aggregate",
            "profiles whose median of site medians is outside the registered log2-ratio bound",
            "entities whose residue letter disagrees",
            "several entities on one site in one profile",
        ],
        "partitions": "n/a: descriptive join, no fitted model",
    }
    path = save_result("phosphosite_tumour_abundance", payload, manifest=manifest)
    print(json.dumps({k: v for k, v in payload.items() if k not in ("sites",)}, indent=1)[:12000])
    print(path)


def main() -> None:
    CACHE.mkdir(parents=True, exist_ok=True)
    stage = sys.argv[1] if len(sys.argv) > 1 else "map"
    if stage == "map":
        print(json.dumps(map_sites(), indent=1))
    elif stage == "join":
        join()
    else:
        raise SystemExit("usage: cptac_phospho.py map|join")


if __name__ == "__main__":
    main()
