# SPDX-License-Identifier: AGPL-3.0-or-later
"""A neural and embryonic biosample for reader v1: fetal spinal cord (area B, lane-readerneural, 2026-10-10).

    uv run python scripts/reader_neural.py register   # the choice and the sanity check, before any fetch
    uv run python scripts/reader_neural.py fetch      # marks and DNase peaks; the layer's manifest row
    uv run python scripts/reader_genome_wide.py "spinal cord"   # the reader, as for every biosample
    uv run python scripts/reader_neural.py check      # the registered sanity check, run once

Reader v1 read thirteen biosamples. Its neural ones were a neuroblastoma line (SK-N-SH) and cultured
astrocytes, its embryonic one a stem cell line (H1), and its only tissues were adult testis and ovary:
no primary tissue of the nervous system and nothing embryonic that was not a cell line. This adds one
biosample, chosen in writing before anything was fetched (`register`), read by the unchanged reader
(`scripts/reader_genome_wide.py`, the same `read_chromosome`, the same thresholds), and checked once
against a sanity test registered before any of its readings existed (`check`).

The sanity test asks one question: does the new biosample's reader call the promoters of a closed,
cited list of neural tube genes read more often than K562's, beyond what each biosample's genome-wide
read share predicts? It checks the biosample's label. Accessibility is not function: it says nothing
about what any of these genes does.

What this script does not change: `genomeos/attribution/context_evidence.READER_TERMS`, and through
it `compile.reader_cells` and `readergate.CELLS`, stay at thirteen. Those feed registered tests over
thirteen cells (readergate, multrna); widening them is a decision for their owners, and the
registration records what it would take.

No model is called and no key is used. Data: ENCODE portal metadata and released GRCh38 peak files
(ENCODE data use policy).
"""

from __future__ import annotations

import json
import math
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402

ENTRY = "scripts/reader_neural.py"
OWN_CODE = (
    ENTRY,
    "scripts/reader_genome_wide.py",
    "genomeos/genome/reader.py",
    "genomeos/genome/epigenome.py",
)
REGISTRATION = "reader_neural_registration"
RESULT = "reader_neural_result"

BIOSAMPLE = "spinal cord"  # ENCODE biosample_ontology.term_name, as reader.find_dnase_file queries it
TERM = "UBERON:0002240"
COMPARATOR = "K562"
ALPHA = 0.05
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

#: The DNase file the reader's own rule (reader.find_dnase_file: the largest released GRCh38 DNase
#: narrowPeak for the term) picked on 2026-10-10. `register` and `fetch` refuse if the portal moved.
EXPECTED_DNASE = "ENCFF041JKB"

#: The registered gene list, closed before any spinal cord reading existed. Human symbols (HGNC) with
#: the chromosome GENCODE v50 places each on. The seven progenitor homeodomain proteins of the ventral
#: neural tube code and the markers of the five ventral neuron classes they specify.
GENES: tuple[tuple[str, str, str], ...] = (
    ("PAX7", "chr1", "progenitor, class I (dorsal-ventral code)"),
    ("IRX3", "chr16", "progenitor, class I"),
    ("DBX1", "chr11", "progenitor, class I"),
    ("DBX2", "chr12", "progenitor, class I"),
    ("PAX6", "chr11", "progenitor, class I"),
    ("NKX6-1", "chr4", "progenitor, class II"),
    ("NKX2-2", "chr20", "progenitor, class II"),
    ("EVX1", "chr7", "V0 interneuron marker"),
    ("EN1", "chr2", "V1 interneuron marker"),
    ("VSX2", "chr14", "V2 interneuron marker (Chx10)"),
    ("MNX1", "chr7", "motor neuron marker (HB9)"),
    ("ISL1", "chr5", "motor neuron marker"),
    ("SIM1", "chr6", "V3 interneuron marker"),
)
CITATIONS = (
    "Briscoe J, Pierani A, Jessell TM, Ericson J (2000). A homeodomain protein code specifies "
    "progenitor cell identity and neuronal fate in the ventral neural tube. Cell 101:435-445. "
    "doi:10.1016/S0092-8674(00)80853-3 (class I PAX7, IRX3, DBX1, DBX2, PAX6; class II NKX6-1, NKX2-2)",
    "Jessell TM (2000). Neuronal specification in the spinal cord: inductive signals and "
    "transcriptional codes. Nat Rev Genet 1:20-29. doi:10.1038/35049541 (the progenitor code and the "
    "markers of the ventral neuron classes V0 EVX1, V1 EN1, V2 CHX10/VSX2, MN HB9/MNX1 and ISL1, V3 SIM1)",
)

#: The three verdicts, word for word; `check` fills the numbers and nothing else.
WORDING = {
    "PASS": (
        "Sanity check passed. Fetal spinal cord's reader calls {a} of {n} registered ventral neural "
        "tube genes read and K562's {k} of {n} ({b} read only in spinal cord, {c} only in K562; "
        "one-sided exact sign test p = {p}). Over the genome the two read {sc_read:,} of {sc_coding:,} "
        "and {k_read:,} of {k_coding:,} coding genes, so the excess at these genes is {e_sc:+.3f} in "
        "spinal cord against {e_k:+.3f} in K562. The new biosample's promoter calls look like neural "
        "tissue's at these genes. This checks the biosample's label; accessibility is not function, "
        "and it says nothing about what any of these genes does."
    ),
    "FAIL": (
        "Sanity check failed ({why}). Fetal spinal cord's reader calls {a} of {n} registered ventral "
        "neural tube genes read and K562's {k} of {n} ({b} read only in spinal cord, {c} only in K562; "
        "one-sided exact sign test p = {p}); over the genome the two read {sc_read:,} of {sc_coding:,} "
        "and {k_read:,} of {k_coding:,} coding genes (excess {e_sc:+.3f} against {e_k:+.3f}). The new "
        "biosample's reader results are not to be used in any comparison until the failure is explained."
    ),
    "UNDECIDED": (
        "Sanity check undecided. The direction is as registered: fetal spinal cord's reader calls {a} "
        "of {n} registered ventral neural tube genes read and K562's {k} of {n}, and the excess over "
        "each genome-wide read share ({sc_read:,} of {sc_coding:,}; {k_read:,} of {k_coding:,}) is "
        "{e_sc:+.3f} against {e_k:+.3f}; but {b} genes read only in spinal cord against {c} only in "
        "K562 cannot be told from chance at the registered 0.05 (one-sided exact sign test p = {p}). "
        "The biosample may be used with this stated beside any use. Accessibility is not function."
    ),
}


# ---- the registered statistic ------------------------------------------------------------------


def sign_test_p(b: int, c: int) -> float:
    """One-sided exact sign test: P(X >= b) for X ~ Binomial(b + c, 1/2); 1.0 with no discordant gene."""
    n = b + c
    if n == 0:
        return 1.0
    return sum(math.comb(n, i) for i in range(b, n + 1)) / 2**n


def verdict(counts: dict[str, int]) -> dict[str, Any]:
    """The registered decision from the counts: a, k (read in spinal cord, in K562) of n; b, c (read in
    one only); and each biosample's genome-wide read and coding totals."""
    n = counts["n"]
    p = sign_test_p(counts["b"], counts["c"])
    share_sc = counts["sc_read"] / counts["sc_coding"]
    share_k = counts["k_read"] / counts["k_coding"]
    e_sc = counts["a"] / n - share_sc
    e_k = counts["k"] / n - share_k
    if counts["b"] <= counts["c"]:
        outcome, why = "FAIL", "spinal cord does not read more of the listed genes than K562"
    elif e_sc <= e_k:
        outcome, why = "FAIL", "the excess over the genome-wide read share is not larger in spinal cord"
    elif p > ALPHA:
        outcome, why = "UNDECIDED", ""
    else:
        outcome, why = "PASS", ""
    fields = {**counts, "p": round(p, 5), "e_sc": e_sc, "e_k": e_k, "why": why}
    return {
        "outcome": outcome,
        "p_one_sided": round(p, 5),
        "excess_spinal_cord": round(e_sc, 4),
        "excess_K562": round(e_k, 4),
        "genome_wide_read_share_spinal_cord": round(share_sc, 4),
        "genome_wide_read_share_K562": round(share_k, 4),
        "reading": WORDING[outcome].format(**fields),
    }


def gene_reads(reading: dict[str, Any], symbol: str) -> bool:
    """The reader's own per-gene call: a protein-coding gene of the chromosome's reading is read unless
    `silent_genes` (closed and poised promoters) lists it."""
    return symbol not in set(reading["silent_genes"])


# ---- helpers -----------------------------------------------------------------------------------


def _git(*args: str) -> str:
    return subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, text=True, check=False
    ).stdout.strip()


def _manifest(sources: list[dict], inputs: list[dict], parameters: dict, partitions: Any) -> dict:
    return {
        "sources": sources,
        "inputs": inputs,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "code": {},
        "parameters": parameters,
        "exclusions": [],
        "partitions": partitions,
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


ENCODE_SOURCE = {
    "accession": "ENCODE portal (https://www.encodeproject.org): DNase-seq and Histone ChIP-seq metadata "
    "and released GRCh38 peak files",
    "version": "queried 2026-10-10",
    "licence": "ENCODE data use policy (https://www.encodeproject.org/help/citing-encode/)",
}


def _committed(name: str) -> str:
    """The commit that last touched a tracked result, refusing one that is untracked or modified."""
    rel = f"data/results/{name}.json"
    if not _git("ls-files", rel):
        raise SystemExit(f"{rel} is not tracked: commit the registration before this step")
    if _git("status", "--porcelain", "--", rel):
        raise SystemExit(f"{rel} differs from HEAD: the step runs only against the committed registration")
    return _git("log", "-1", "--format=%H", "--", rel)


# ---- register ----------------------------------------------------------------------------------


def _survey() -> dict[str, Any]:
    """What the portal held on 2026-10-10 for the neural candidates, as the choice weighed it."""
    return {
        "brain": "DNase peak files for the term mix adult (53 to 90+ years), fetal (56 to 142 days) and "
        "organoid samples, and the term has no H3K27ac experiment, so the reader would run without "
        "its marks on a file the size rule draws from that mixture",
        "dorsolateral prefrontal cortex": "adult cortex with DNase and all three marks, but 114 of the "
        "181 Histone ChIP-seq experiments for the three marks carry a diagnosis (Alzheimer's disease, "
        "mild or other cognitive impairment) and the reader's rules (largest DNase file; marks from the "
        "biosample summary covering most marks) cannot exclude a diagnosed donor",
        "head of caudate nucleus, posterior cingulate gyrus, frontal cortex, middle frontal gyrus": "DNase "
        "without an H3K27ac experiment",
        "caudate nucleus": "adult striatum with DNase and all three marks (Roadmap, donors aged 75 and "
        "81 years; the DNase from other adult donors): adult brain but not early development. It is the "
        "recorded next candidate if an adult brain biosample is wanted",
        "cerebellum, midbrain, medulla oblongata, globus pallidus, Ammon's horn, superior temporal gyrus, "
        "occipital lobe, inferior parietal cortex, pons, putamen, cerebellar cortex": "no human GRCh38 "
        "H3K27ac or H3K27me3 experiment for the term (cerebellum has one human H3K4me3)",
        "neural progenitor cell (from H9), neuronal stem cell (from H1), neural cell (from H1)": "all "
        "three marks, but in vitro differentiated from a stem line rather than tissue, and two derive "
        "from H1, a line the reader already holds; neural cell has no DNase",
        "neurosphere": "cultured, no DNase peak file",
    }


def register() -> None:
    from genomeos.genome import epigenome as ep
    from genomeos.genome import reader
    from genomeos.results import load_result, save_result

    if any((ROOT / "data/results").glob(f"dnase_{reader.slug(BIOSAMPLE)}_chr*.bed.gz")):
        raise SystemExit("spinal cord peaks are already on disk: a registration must precede the fetch")
    panel = list((load_result("reader_genome_wide") or {}).get("cell_types") or [])
    if BIOSAMPLE in panel:
        raise SystemExit("spinal cord is already on the reader's panel: a registration must precede it")
    dnase = reader.find_dnase_file(BIOSAMPLE)
    if dnase["accession"] != EXPECTED_DNASE:
        raise SystemExit(f"the reader's DNase rule now picks {dnase['accession']}, not {EXPECTED_DNASE}")
    marks = ep.pick_marks({m: ep.mark_options(BIOSAMPLE, m) for m in ep.MARKS})
    methylation = ep.resolve_methylation(BIOSAMPLE)
    peaks_bytes = int(dnase["size"] or 0) + sum(
        int((e.get("peaks") or {}).get("size") or 0) for e in marks.values()
    )
    free = shutil.disk_usage(ROOT).free
    registration = {
        "registered": "2026-10-10, before any spinal cord file was fetched and before any spinal cord "
        "reading existed; written by scripts/reader_neural.py register",
        "item": "ROADMAP area B: a neural or embryonic biosample beyond SK-N-SH and H1",
        "choice": {
            "biosample": BIOSAMPLE,
            "ontology_term": TERM,
            "classification": "tissue; the term's DNase-seq experiments (6) and its Histone ChIP-seq "
            "experiments for the three reader marks (3) are all from human embryos, 59 to 113 days",
            "why": [
                "the reader held no primary tissue of the nervous system and nothing embryonic that was "
                "not a cell line: its neural biosamples were SK-N-SH (a neuroblastoma line) and astrocyte "
                "(cultured primary cells), its embryonic one H1 (a stem line), its tissues adult testis "
                "and ovary",
                "of the human nervous system tissue terms with released GRCh38 DNase-seq peaks, three "
                "also have released GRCh38 Histone ChIP-seq for all three marks the reader uses (H3K27me3, "
                "H3K27ac, H3K4me3): spinal cord (all embryonic), caudate nucleus (adult) and dorsolateral "
                "prefrontal cortex (adult, mostly diagnosed donors). Spinal cord is the one that is neural "
                "and embryonic at once, so one biosample answers both halves of the item",
            ],
            "alternatives_considered": _survey(),
            "dnase": {
                **{k: dnase[k] for k in ("accession", "size", "candidates")},
                "rule": "reader.find_dnase_file: the largest released GRCh38 DNase narrowPeak for the term",
                "experiment": "ENCSR910OQF",
                "replicates": "biological replicate 1 of the experiment's 2 (replicate 1: female embryo, "
                "59 days, donor ENCDO294EUZ; replicate 2, not read: male embryo, 72 days)",
                "assay": "DNase-seq, peaks (bed narrowPeak), GRCh38",
            },
            "marks": {
                m: {
                    "experiment": e.get("experiment"),
                    "biosample": e.get("biosample"),
                    "peaks": e.get("peaks"),
                    "value": e.get("value"),
                }
                for m, e in marks.items()
            },
            "methylation": {k: methylation.get(k) for k in ("experiment", "biosample", "value", "reason")},
            "caveats": [
                "the DNase comes from a 59-day embryo and the marks from a 108-day (H3K27me3, H3K27ac, "
                "H3K4me1, H3K9me3) and a 113-day embryo (H3K4me3): the poised and read-by-marks calls "
                "combine chromatin of different embryos at different ages",
                "every mark is one biological replicate with pseudoreplicated peaks",
                "bulk tissue: a mixture of progenitors, neurons, glia and non-neural cells",
                "the DNase donor is female, so chrY carries only stray peaks",
            ],
            "licence": "ENCODE data use policy: released data, cite ENCODE and the Roadmap Epigenomics "
            "Consortium",
            "size_check": {
                "peak_bytes_to_fetch": peaks_bytes,
                "not_fetched": "fold-change signal bigWigs and WGBS (about 1.7 GB and 0.8 GB); the reader "
                "reads peaks only",
                "limit_bytes": 1_500_000_000,
                "disk_free_bytes_at_registration": free,
            },
        },
        "run": "scripts/reader_genome_wide.py 'spinal cord': the unchanged read_chromosome over the "
        "local chromosomes, the same promoter window (TSS +/- 1 kb), the same poised and read-by-marks "
        "rules, the same node thresholds; the per-chromosome readings stay local except chr21, as for "
        "every biosample",
        "sanity_check": {
            "question": "does the new biosample's reader call the promoters of a closed, cited list of "
            "ventral neural tube genes read more often than K562's, beyond each biosample's genome-wide "
            "read share?",
            "genes": [{"symbol": s, "chrom": c, "role": r} for s, c, r in GENES],
            "citations": list(CITATIONS),
            "comparator": f"{COMPARATOR}, the reader's existing reading "
            "(data/results/reader_K562_<chrom>.json)",
            "read_call": "the reader's own per-gene call, unchanged: a gene that is protein_coding once on "
            "its registered chromosome in GENCODE v50 is read when that chromosome's reading does not "
            "list it in silent_genes (closed and poised promoters); read by its marks counts as read",
            "denominator": "n = the listed genes found once as protein_coding on their registered "
            "chromosome; a gene not found is excluded and named",
            "counts": "a, k = genes read in spinal cord, in K562; b = read in spinal cord and not K562; "
            "c = read in K562 and not spinal cord",
            "genome_wide_read_share": "sum of genes_read over sum of coding_genes over the 24 per-chromosome "
            "readings of each biosample",
            "excess": "e = a/n minus spinal cord's genome-wide read share; e_K = k/n minus K562's",
            "test": f"one-sided exact sign test on the discordant genes, P(X >= b), X ~ Binomial(b + c, "
            f"1/2), alpha {ALPHA}",
            "decision": {
                "PASS": "b > c, e > e_K and p <= 0.05",
                "FAIL (the falsifier)": "b <= c, or e <= e_K",
                "UNDECIDED": "b > c and e > e_K, but p > 0.05",
            },
            "wording": WORDING,
            "on_fail": "the new biosample's reader results are not used in any comparison until the "
            "failure is explained; nothing is re-run with another list, comparator or threshold",
            "descriptive_only": "the same counts for every other biosample on the panel are reported "
            "beside the verdict, labelled descriptive; they are not a test",
        },
        "not_claimed": [
            "accessibility is not function: an open or marked promoter is necessary for transcription, "
            "not proof of it",
            "the check tests the biosample's label at thirteen genes, not the reader's accuracy in it",
        ],
        "consumers_left_unchanged": {
            "context_evidence.READER_TERMS": "thirteen; spinal cord would enter as "
            f"'{BIOSAMPLE}': '{TERM}', a term AlphaGenome's track metadata carries (embryonic spinal cord "
            "DNase and RNA-seq tracks)",
            "compile.reader_cells and readergate.CELLS": "derive from READER_TERMS, so they follow it; "
            "readergate's registered test is over thirteen cells and needs a measured RNA experiment per "
            "biosample (ENCODE holds polyA plus RNA-seq for embryonic spinal cord); widening it is a new "
            "registration, not an edit",
        },
        "panel_before": panel,
    }
    reg_inputs = [mf.input_entry("data/results/reader_genome_wide.json")]
    m = _manifest(
        [ENCODE_SOURCE],
        reg_inputs,
        {"alpha": ALPHA, "genes": len(GENES), "comparator": COMPARATOR, "expected_dnase": EXPECTED_DNASE},
        "n/a: a registration, no evaluation split",
    )
    save_result(REGISTRATION, registration, manifest=m)
    print(
        f"registered: {BIOSAMPLE}, DNase {dnase['accession']}, {peaks_bytes:,} peak bytes, "
        f"{free / 2**30:.1f} GiB free",
        flush=True,
    )


# ---- fetch -------------------------------------------------------------------------------------


def fetch() -> None:
    from genomeos.genome import default_gencode, reader
    from genomeos.genome import epigenome as ep
    from genomeos.results import load_result, save_result

    reg_sha = _committed(REGISTRATION)
    reg = load_result(REGISTRATION)
    free = shutil.disk_usage(ROOT).free
    if free < 15 * 2**30:
        raise SystemExit(f"{free / 2**30:.1f} GiB free, under the 15 GiB floor: paused")
    picked = ep.pick_marks({m: ep.mark_options(BIOSAMPLE, m) for m in ep.MARKS})
    for mark, e in picked.items():
        want = (reg["choice"]["marks"].get(mark) or {}).get("peaks") or {}
        got = e.get("peaks") or {}
        if want.get("accession") != got.get("accession"):
            raise SystemExit(
                f"{mark}: the portal now picks {got.get('accession')}, registered {want.get('accession')}"
            )
    for mark, e in picked.items():
        if e.get("peaks"):
            kept = ep.fetch_mark_peaks({"cell_type": BIOSAMPLE, **e})
            print(f"{BIOSAMPLE} {mark}: {e['peaks']['accession']} {sum(kept.values()):,} peaks", flush=True)
    # the chromosomes the roll-up reads, computed as it computes them, so it finds every file cached
    local = [c for c in CHROMS if (ROOT / f"data/reference/{c}.fa.gz").exists() and default_gencode({c})]
    info = reader.fetch_peaks(BIOSAMPLE, set(local))
    if info["accession"] != EXPECTED_DNASE:
        raise SystemExit(f"DNase {info['accession']} fetched, registered {EXPECTED_DNASE}")
    print(
        f"{BIOSAMPLE} DNase {info['accession']}: {sum(info['kept'].values()):,} peaks on {len(local)} "
        "chromosomes",
        flush=True,
    )
    # the layer's manifest gains one row, built by the functions build_manifest uses; the other rows
    # are carried as they are, and the coverage table is recomputed only after checking that the
    # recomputation reproduces the stored one
    before = ep.load_manifest()
    if BIOSAMPLE in before["cell_types"]:
        raise SystemExit("the epigenome manifest already has a spinal cord row")
    if ep.coverage_table(before["cell_types"]) != before["coverage"]:
        raise SystemExit("the stored coverage table is not what coverage_table gives for the stored rows")
    cells = dict(before["cell_types"])
    cells[BIOSAMPLE] = {
        "dnase": ep.dnase_accession(BIOSAMPLE),
        "marks": {m: {"cell_type": BIOSAMPLE, **e} for m, e in picked.items()},
        "methylation": ep.resolve_methylation(BIOSAMPLE),
    }
    payload = {"source": before["source"], "rule": before["rule"], "cell_types": cells}
    payload["coverage"] = ep.coverage_table(cells)
    fetched = sorted((ROOT / "data/results").glob(f"dnase_{reader.slug(BIOSAMPLE)}_chr*.bed.gz"))
    fetched += sorted((ROOT / ep.CACHE / "peaks").glob(f"{ep.slug(BIOSAMPLE)}_*_chr*.bed.gz"))
    m = _manifest(
        [ENCODE_SOURCE],
        [
            mf.input_entry("data/results/epigenome_manifest.json"),
            mf.files_entry(f"{BIOSAMPLE} DNase and mark peaks as fetched", fetched),
        ],
        {"biosample": BIOSAMPLE, "registration_commit": reg_sha},
        "n/a: a manifest of fetched files, no evaluation split",
    )
    save_result(ep.MANIFEST, payload, manifest=m)
    print(f"epigenome manifest: {len(cells)} biosamples", flush=True)


# ---- check -------------------------------------------------------------------------------------


def check() -> None:
    from genomeos.genome import Annotation, default_gencode
    from genomeos.genome.reader import slug
    from genomeos.results import load_result, save_result

    if load_result(RESULT) is not None:
        raise SystemExit("the sanity check has a result already; it is run once")
    reg_sha = _committed(REGISTRATION)
    reg = load_result(REGISTRATION)
    registered = [(g["symbol"], g["chrom"]) for g in reg["sanity_check"]["genes"]]
    if registered != [(s, c) for s, c, _ in GENES] or reg["sanity_check"]["wording"] != WORDING:
        raise SystemExit("the code's gene list or wording differs from the committed registration")
    reading_sha = _committed(f"reader_{slug(BIOSAMPLE)}_chr21")
    after = subprocess.run(
        ["git", "merge-base", "--is-ancestor", reg_sha, reading_sha], cwd=ROOT, check=False
    )
    if after.returncode != 0 or reg_sha == reading_sha:
        raise SystemExit("the spinal cord reading was not committed after the registration")
    gw = load_result("reader_genome_wide")
    panel = list(gw["cell_types"])
    if BIOSAMPLE not in panel:
        raise SystemExit("spinal cord is not on the reader's panel yet: run the roll-up first")

    read_files: list[Path] = []

    def reading(cell: str, chrom: str) -> dict[str, Any] | None:
        p = ROOT / f"data/results/reader_{slug(cell)}_{chrom}.json"
        if not p.exists():
            return None
        read_files.append(p)
        return json.loads(p.read_text())

    # the denominator: listed genes found once as protein_coding on their registered chromosome
    gffs, found, excluded = [], [], []
    for symbol, chrom, _ in GENES:
        gff = default_gencode({chrom})
        gffs.append(gff)
        ann = Annotation.from_gff3(gff, {chrom})
        hits = [g for g in ann.genes.values() if g.symbol == symbol and g.locus.chrom == chrom]
        if len(hits) == 1 and hits[0].type == "protein_coding":
            found.append((symbol, chrom))
        else:
            excluded.append({"symbol": symbol, "chrom": chrom, "found": [h.type for h in hits]})

    def calls(cell: str) -> dict[str, Any] | None:
        per_gene, cache = {}, {}
        for symbol, chrom in found:
            if chrom not in cache:
                cache[chrom] = reading(cell, chrom)
            if cache[chrom] is None:
                return None
            per_gene[symbol] = gene_reads(cache[chrom], symbol)
        totals = [reading(cell, c) for c in CHROMS]
        if any(t is None for t in totals):
            return None
        return {
            "genes": per_gene,
            "read": sum(per_gene.values()),
            "genome_read": sum(t["genes_read"] for t in totals),
            "genome_coding": sum(t["coding_genes"] for t in totals),
        }

    sc, k = calls(BIOSAMPLE), calls(COMPARATOR)
    if sc is None or k is None:
        raise SystemExit("a per-chromosome reading of spinal cord or K562 is missing on this machine")
    b = sum(1 for s, _ in found if sc["genes"][s] and not k["genes"][s])
    c = sum(1 for s, _ in found if k["genes"][s] and not sc["genes"][s])
    counts = {
        "n": len(found),
        "a": sc["read"],
        "k": k["read"],
        "b": b,
        "c": c,
        "sc_read": sc["genome_read"],
        "sc_coding": sc["genome_coding"],
        "k_read": k["genome_read"],
        "k_coding": k["genome_coding"],
    }
    v = verdict(counts)
    descriptive = {}
    for cell in panel:
        if cell in (BIOSAMPLE, COMPARATOR):
            continue
        d = calls(cell)
        descriptive[cell] = (
            None
            if d is None
            else {
                "read_of_listed": f"{d['read']} of {len(found)}",
                "genome_wide": f"{d['genome_read']:,} of {d['genome_coding']:,}",
                "excess": round(d["read"] / len(found) - d["genome_read"] / d["genome_coding"], 4),
            }
        )
    result = {
        "registration": f"data/results/{REGISTRATION}.json at {reg_sha}",
        "reading_committed_at": reading_sha,
        "biosample": BIOSAMPLE,
        "comparator": COMPARATOR,
        "counts": counts,
        "per_gene": {s: {"spinal_cord": sc["genes"][s], "K562": k["genes"][s]} for s, _ in found},
        "excluded_genes": excluded,
        **v,
        "descriptive_other_biosamples": {
            "note": "not a test and not part of the verdict: the same counts for the rest of the panel",
            "rows": descriptive,
        },
    }
    m = _manifest(
        [
            ENCODE_SOURCE,
            {"accession": "GENCODE v50 comprehensive annotation (GFF3)", "version": "v50"},
            {
                "accession": "reader v1 per-chromosome readings "
                "(data/results/reader_<biosample>_<chrom>.json)",
                "version": "on disk at the run; each file's sha256 is in inputs",
            },
        ],
        [
            mf.input_entry(f"data/results/{REGISTRATION}.json"),
            mf.input_entry("data/results/reader_genome_wide.json"),
            mf.files_entry("GENCODE v50 GFF3 of the listed genes' chromosomes", sorted(set(map(str, gffs)))),
            mf.files_entry("reader per-chromosome readings", sorted(set(map(str, read_files)))),
        ],
        {"alpha": ALPHA, "genes": len(GENES), "comparator": COMPARATOR},
        "n/a: a registered check on fixed genes, no evaluation split",
    )
    save_result(RESULT, result, manifest=m)
    print(v["reading"], flush=True)


def main() -> None:
    steps = {"register": register, "fetch": fetch, "check": check}
    if len(sys.argv) != 2 or sys.argv[1] not in steps:
        raise SystemExit(__doc__)
    steps[sys.argv[1]]()


if __name__ == "__main__":
    main()
