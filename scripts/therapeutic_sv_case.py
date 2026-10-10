# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""The eleventh benchmark case: a structural variant, registered before it was run.

The roadmap row for cancer and therapeutics has kept three items open since
milestone 1.2, and this is one of them: "a scored structural-variant case". The
ten-case benchmark does contain a fusion, EML4-ALK, and it is the arrangement
where the answer is already settled — ALK is the 3' partner, the product begins
with EML4 and keeps none of ALK's outward-facing part, and the score says 0.0
with that reason. The arrangement the benchmark has never scored is the other
one, where the gene under test is the fusion's 5' partner: its own N-terminus,
its signal peptide and the start of its ectodomain ARE in the product, and
whether its membrane anchor is depends on where the junction falls.

The case, chosen and written down before the pipeline was run on it:

  FGFR3-TACC3 in urothelial carcinoma, sample TCGA-CF-A3MG-01 of
  blca_tcga_pan_can_atlas_2018 (TCGA Bladder Urothelial Carcinoma, PanCancer
  Atlas), as the cBioPortal public REST API returns it. FGFR3 is site1 and so
  the 5' partner, TACC3 site2 and the 3'; chr4 1,806,934 to 1,739,702 on
  GRCh37; in-frame; somatic. Erdafitinib, a small-molecule pan-FGFR inhibitor,
  carries an FDA label for locally advanced or metastatic urothelial carcinoma
  with susceptible FGFR3 alterations, and the approved companion diagnostic's
  susceptible alterations include the FGFR3-TACC3 fusion — so the structural
  call is the measurement that selects the patient, the way the amplification
  is for trastuzumab. No antibody against FGFR3 is approved.

Two commands, in this order, and the registration is committed between them:

  register  writes data/results/therapeutic_sv_case_registration.json: what the
            pipeline must produce for the case to count as recovered, what
            counts as a miss, that the ten existing rows must not move, and
            that the case is added to the benchmark's case list whatever the
            outcome.
  run       runs the case the way the benchmark runs a case (`run_case`, with
            the network, which is the registered mode), re-scores the ten
            existing cases offline with the structural-variant reader as it was
            before this lane's change and as it is now, and compares every
            field of the ten against each other and against the committed
            data/results/therapeutic_benchmark_v2.json. It writes
            data/results/therapeutic_sv_case_result.json and does not rewrite
            the committed benchmark.

Why the ten are re-scored offline and the new case with the network: the
providers answer from the local knowledge stores first and fetch what a store
lacks from live services that publish no version a request can pin, so a
field-for-field comparison has to hold the network out of it. The new case is a
new row and has nothing to be compared with, so it is run in the mode the
benchmark's own result records.

Run: `uv run python scripts/therapeutic_sv_case.py register`
     `uv run python scripts/therapeutic_sv_case.py run [--offline]`
"""

from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.cancer import alterations as alt  # noqa: E402
from genomeos.cancer.alterations import GeneAlteration  # noqa: E402
from genomeos.results import load_result, save_result  # noqa: E402
from genomeos.therapeutics import pipeline as pl  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
ENTRY = "scripts/therapeutic_sv_case.py"
OWN_CODE = frozenset(
    {
        ENTRY,
        "genomeos/cancer/alterations.py",
        "scripts/therapeutic_benchmark.py",
        "tests/test_therapeutic_sv_case.py",
    }
)

BENCHMARK = Path(__file__).with_name("therapeutic_benchmark.py")
COMMITTED = "therapeutic_benchmark_v2"
REGISTRATION = "therapeutic_sv_case_registration"
RESULT = "therapeutic_sv_case_result"
CASE_GENE = "FGFR3"
SAMPLE = "TCGA-CF-A3MG-01"
STUDY = "blca_tcga_pan_can_atlas_2018"
FIXTURES = ROOT / "data/demo/alterations"

#: Every row field but these is compared between the before and after runs.
#: `seconds` is wall-clock time and says nothing about the pipeline's answer.
NOT_COMPARED = frozenset({"seconds"})


def _bench():
    spec = importlib.util.spec_from_file_location("therapeutic_benchmark", BENCHMARK)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _case(bench) -> dict:
    return next(c for c in bench.CASES if c["gene"] == CASE_GENE and c["call"] == "structural variant")


def _ten(bench) -> list[dict]:
    case = _case(bench)
    return [c for c in bench.CASES if c is not case]


# ---- the reader as it was before this lane's change ---------------------------------


def read_sv_table_positional(path: str) -> list[GeneAlteration]:
    """`read_sv_table` as `0ebc4a9` had it: four columns in a fixed order.

    Kept here, rather than read out of git at run time, so that the before side
    of the comparison is a function this file can be read against. It is the
    body of the positional branch and nothing else.
    """
    p = Path(path)
    if not p.exists():
        return []
    out: list[GeneAlteration] = []
    for parts in alt.rows(p):
        gene = parts[0].upper()
        if gene in ("GENE", "GENE1", "SITE1"):
            continue
        partner = parts[1].upper() if len(parts) > 1 else ""
        orientation = ""
        if len(parts) > 3:
            raw = parts[3].strip().rstrip("'").rstrip("’")
            orientation = {"5": "5'", "3": "3'"}.get(raw, "")
        a = GeneAlteration(
            gene,
            "fusion",
            partner=partner,
            detail=parts[2] if len(parts) > 2 else "",
            orientation=orientation,
            source=f"patient structural-variant table ({p.name})",
        )
        a.evidence.append(f"{a.label()} in this tumour's structural-variant table")
        out.append(a)
    return out


class _Reader:
    """The structural-variant reader the pipeline uses, as built or as it was."""

    def __init__(self) -> None:
        self.built = pl.read_sv_table

    def __call__(self, built: bool) -> None:
        pl.read_sv_table = self.built if built else read_sv_table_positional


# ---- the registration ---------------------------------------------------------------


def registration(bench) -> dict:
    case = _case(bench)
    out = {
        "result": REGISTRATION,
        "date": "2026-10-10",
        "lane": "lane-svcase",
        "roadmap_item": 'area F, "a scored structural-variant case" (docs/ROADMAP.md section 3)',
        "registered_before_the_run": True,
        "case": {
            "case": case["case"],
            "gene": CASE_GENE,
            "partner": "TACC3",
            "sample": SAMPLE,
            "study": STUDY,
            "study_name": "TCGA Bladder Urothelial Carcinoma, PanCancer Atlas (2018)",
            "source": (
                "cBioPortal public REST API, structural-variant profile "
                f"{STUDY}_structural_variants and mutation profile {STUDY}_mutations, read 2026-10-10"
            ),
            "terms": (
                "TCGA open-access somatic data, served by cBioPortal, whose terms ask that the study "
                "and the portal be cited (Cerami et al. 2012; Gao et al. 2013); no controlled-access "
                "data and no patient identifier beyond the TCGA sample barcode"
            ),
            "record": {
                "site1_gene": "FGFR3",
                "site1_end": "5'",
                "site1": "chr4:1,806,934 (GRCh37)",
                "site2_gene": "TACC3",
                "site2_end": "3'",
                "site2": "chr4:1,739,702 (GRCh37)",
                "frame": "in-frame (site2EffectOnFrame)",
                "status": "SOMATIC",
                "event": "FGFR3-TACC3",
                "cohort": "8 FGFR3-TACC3 records in this study, 6 of them in-frame",
            },
            "vcf": {
                "variants": [
                    "ARID1A R1287T, chr1:26,773,490 G>C (GRCh38), allele fraction 0.208 from 5 alt "
                    "and 19 ref reads",
                    "TSC1 R190P, chr9:132,921,913 C>G (GRCh38), allele fraction 0.800 from 36 alt "
                    "and 9 ref reads",
                ],
                "lift": (
                    "the record is GRCh37 and the other benchmark fixtures are GRCh38, so both "
                    "positions were lifted through the Ensembl REST assembly map; each lifted "
                    "position carries the record's reference allele and annotates to the protein "
                    "change the record names (checked against Ensembl VEP REST before the fixture "
                    "was committed)"
                ),
                "no_fgfr3_variant": (
                    "FGFR3 is not mutated in this sample's 125 somatic mutations, so the structural "
                    "call is the only route to it"
                ),
            },
            "published_target": (
                "erdafitinib, a small-molecule pan-FGFR tyrosine kinase inhibitor: FDA accelerated "
                "approval 12 April 2019 (BLC2001) and conversion to full approval 19 January 2024 "
                "(THOR) for locally advanced or metastatic urothelial carcinoma with susceptible "
                "FGFR3 genetic alterations, determined by an approved companion diagnostic "
                "(therascreen FGFR RGQ RT-PCR), whose susceptible alterations include the "
                "FGFR3-TACC3 fusion. No antibody against FGFR3 is approved"
            ),
            "modality": "small_molecule",
            "expect": "out_of_scope_expected",
            "different_in_kind": (
                "the benchmark's existing fusion case has the gene under test as the 3' partner, "
                "where the product keeps none of its outward-facing part and the answer is settled at "
                "0.0. Here the gene under test is the 5' partner, so the product carries its own "
                "N-terminus and the curated surface figure may or may not describe the protein the "
                "tumour makes. It is also the first benchmark case read from a caller's own export "
                "rather than from a four-column file written by hand, and the first that is one "
                "patient's record rather than a hand-written file or a cohort summary"
            ),
        },
        "registered": {
            "recovered": {
                "quantity": "recovered, and rank among the ranked candidates",
                "counts_as_recovered": (
                    "FGFR3 appears among the ranked candidates, reached by the structural call alone "
                    "(the VCF carries no FGFR3 variant and there is no copy-number table)"
                ),
                "counts_as_a_miss": "FGFR3 absent from the candidates: TARGET LOST",
                "rank_is_recorded_not_predicted": (
                    "the rank is published and no threshold is registered for it; what is registered "
                    "is that a target reached only through a structural call must be reachable at all"
                ),
            },
            "route": {
                "quantity": "verdict and pass, against expect=out_of_scope_expected",
                "counts_as_correct": (
                    "the pipeline does not claim a surface route: surface_accessibility is None or at "
                    "most 0.5, which is the benchmark's own rule for an out-of-scope case"
                ),
                "counts_as_a_miss": (
                    "surface_accessibility above 0.5 — a surface route offered against a product whose "
                    "extent nothing read, for a target with no approved antibody"
                ),
                "predicted_more_strongly": (
                    "surface_accessibility is None, not 0.0 and not a curated full-length value: 0.0 "
                    "is the 3' partner's answer and FGFR3 is the 5' partner. A published 0.0 here "
                    "would be a mis-reading of the arrangement even though the benchmark's pass rule "
                    "would accept it, and it is registered as a finding rather than as a pass"
                ),
            },
            "top_mechanism": {
                "quantity": "best_mechanism, best_mechanism_established and mechanism_sane",
                "counts_as_correct": (
                    "mechanism_sane is true: no mechanism needing an extracellular epitope heads the "
                    "list while the product's extent is unknown, no agonist heads it against an "
                    "activating driver, and a mechanism that heads it on an unanswered requirement is "
                    "a concern rather than a pass"
                ),
                "counts_as_a_miss": "mechanism_sane false",
                "expected": (
                    "no established preferred mechanism, with the open requirement named in "
                    "requirements_unanswered, since the modelled modalities are antibody-like and the "
                    "surface question is unanswered for this product"
                ),
            },
            "regression_guard": {
                "quantity": "every field of the ten existing benchmark rows",
                "requirement": (
                    "the ten rows are identical, field for field, with the structural-variant reader "
                    "as it was before this lane's change and as it is now; `seconds` is excluded "
                    "because it is wall-clock time"
                ),
                "also_checked": (
                    "the offline before-run is compared with the committed "
                    "data/results/therapeutic_benchmark_v2.json row by row, so that the baseline is "
                    "shown to reproduce the committed record rather than only to agree with itself"
                ),
                "counts_as_a_miss": "any compared field of any of the ten rows differing",
            },
        },
        "eleventh_row": {
            "added": "whatever the outcome",
            "where": (
                "scripts/therapeutic_benchmark.py CASES, committed with this registration and before "
                "the case was scored"
            ),
            "committed_result_not_rewritten": (
                "data/results/therapeutic_benchmark_v2.json keeps the ten rows of the 2026-10-02 run. "
                "Rewriting it means a fresh network run of every case, and the ten rows are this "
                "lane's regression guard; the eleventh row and the guard's measurement are published "
                "in data/results/therapeutic_sv_case_result.json instead, and the next full "
                "benchmark run writes the eleventh row into the benchmark file itself"
            ),
        },
        "falsifiers": [
            "FGFR3 absent from the ranked candidates",
            "any compared field of any of the ten existing rows differing between the before and after runs",
            "a surface route claimed for FGFR3 (surface_accessibility above 0.5)",
            "mechanism_sane false",
        ],
        "not_claimed": (
            "a recovered target is a retrospective consistency between the pipeline's ranking and "
            "what is already approved for that alteration. It is not evidence that any therapy "
            "works, and no sequence, construct, dose or protocol is produced (docs/ROADMAP.md "
            "section 1)"
        ),
        "reader_change": {
            "why_before_the_run": (
                "the record as cBioPortal exports it could not be read at all: the positional reader "
                "took the first column, the sample identifier, for a gene and the fourth, a "
                "chromosome, for the end of the product. The named-column reading was committed "
                "first, with the ALK fixture's parse pinned field for field, and this registration "
                "follows it"
            ),
            "what_it_does_not_do": (
                "it does not read a breakpoint into exon content. The record states both breakpoints "
                "and nothing in the pipeline turns them into the residues the product keeps, which is "
                "why the registered answer for this case's surface dimension is unknown"
            ),
        },
    }
    out["result_manifest"] = _registration_manifest(case)
    return out


def _sources() -> list[dict]:
    return [
        {
            "accession": (
                f"cBioPortal {STUDY} (TCGA Bladder Urothelial Carcinoma, PanCancer Atlas): "
                f"sample {SAMPLE}, the FGFR3-TACC3 structural variant and two of its somatic mutations"
            ),
            "version": (
                "public REST API as read 2026-10-10; the fixture files in data/demo/alterations are "
                "the pin, and their sha256 is in inputs"
            ),
        },
        {
            "accession": "Ensembl REST assembly map and VEP (GRCh37 to GRCh38, fixture check)",
            "version": "live service as read 2026-10-10, unversioned per request",
        },
        {
            "accession": (
                "FDA label for erdafitinib (BALVERSA), urothelial carcinoma with susceptible FGFR3 "
                "alterations"
            ),
            "version": "accelerated approval 12 April 2019; full approval 19 January 2024",
        },
    ]


def _registration_manifest(case: dict) -> dict:
    return {
        "sources": _sources(),
        "inputs": [mf.input_entry(FIXTURES / case["vcf"]), mf.input_entry(FIXTURES / case["sv"])],
        "assembly": (
            "GRCh38 for the VCF; the structural-variant record is the source's own GRCh37, which "
            "nothing in the pipeline reads as a coordinate"
        ),
        "coordinates": {"base": 1, "interval": "closed"},
        "parameters": {"registered_before_the_run": True, "eleventh_row": "added whatever the outcome"},
        "exclusions": [
            "no score is in this file: it is the registration, written before the case was run",
        ],
        "partitions": "n/a: a retrospective case registered before it was run, nothing fitted to it",
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


# ---- the run ------------------------------------------------------------------------


def _compare(before: list[dict], after: list[dict], left: str, right: str) -> list[dict]:
    """Every field of every row held against the other run's."""
    out: list[dict] = []
    for b, a in zip(before, after, strict=True):
        if b["case"] != a["case"]:
            out.append({"case": b["case"], "field": "case", left: b["case"], right: a["case"]})
            continue
        for field in sorted(set(b) | set(a)):
            if field in NOT_COMPARED:
                continue
            if b.get(field) != a.get(field):
                out.append({"case": b["case"], "field": field, left: b.get(field), right: a.get(field)})
    return out


def _rows(bench, cases: list[dict], net: bool, reader_built: bool, reader: _Reader) -> list[dict]:
    reader(reader_built)
    with tempfile.TemporaryDirectory(prefix="genomeos-vep-appends-") as appends:
        return [bench.run_case(c, net, sys.stdout, Path(appends)) for c in cases]


def _manifest(net: bool, bench) -> dict:
    case = _case(bench)
    inputs = [mf.input_entry(FIXTURES / case["vcf"]), mf.input_entry(FIXTURES / case["sv"])]
    for c in _ten(bench):
        folder = ROOT / "data/demo" / c.get("dir", "benchmark")
        for key in ("vcf", "cnv", "sv", "rna"):
            if c.get(key):
                inputs.append(mf.input_entry(folder / c[key], case=c["case"]))
    committed = ROOT / "data/results" / f"{COMMITTED}.json"
    if committed.exists():
        inputs.append(mf.input_entry(committed))
    # The distilled cohort the alterations are graded against, declared as the benchmark
    # declares it: `grade` reads both files for every candidate, so they are inputs to this
    # run as much as the fixtures are.
    for name in ("cancer_msk_impact_2017", "cancer_alterations_msk_impact_2017"):
        p = ROOT / "data/results" / f"{name}.json"
        if p.exists():
            inputs.append(mf.input_entry(p))
    knowledge, optional_absent = bench.declare_knowledge()
    inputs += knowledge
    sources = _sources() + [
        {
            "accession": "GenomeOS federated protein compiler: UniProt, Ensembl, InterPro, PDB, "
            "AlphaFold, Reactome, STRING, HPA",
            "version": "live services, unversioned per request; the local store's sha256 is the pin",
        },
        {
            "accession": "the ten benchmark tumours (data/demo) and their committed result",
            "version": "this repository at the recorded commit",
        },
    ]
    return {
        "sources": sources,
        "inputs": inputs,
        "assembly": (
            "GRCh38 for the VCF and the pipeline; the structural-variant record is the source's own "
            "GRCh37, which nothing in the pipeline reads as a coordinate"
        ),
        "coordinates": {
            "base": 1,
            "interval": "closed",
            "note": (
                "VCF positions are 1-based; the breakpoints are the record's own 1-based GRCh37 "
                "positions and are carried verbatim into the alteration's junction field"
            ),
        },
        "parameters": {
            "case_network": net,
            "ten_rows_network": False,
            "top_genes": 8,
            "indirect": True,
            "fields_not_compared": sorted(NOT_COMPARED),
        },
        "exclusions": [
            "small-molecule routes are not modelled: this case passes by not claiming a surface route",
            "the committed benchmark result is not rewritten by this run",
            "no breakpoint is read into exon content, so the product's extent stays unknown",
        ],
        "partitions": "n/a: a retrospective case registered before it was run, nothing fitted to it",
        "inputs_optional_absent": optional_absent,
        "offline_reproducible": not net,
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


def run(net: bool) -> int:
    bench = _bench()
    reader = _Reader()
    case = _case(bench)
    ten = _ten(bench)
    print(f"ten rows, offline, reader as before this lane's change ({len(ten)} cases)")
    before = _rows(bench, ten, False, False, reader)
    print("ten rows, offline, reader as built")
    after = _rows(bench, ten, False, True, reader)
    differ = _compare(before, after, "before", "after")

    committed = load_result(COMMITTED) or {}
    by_case = {r["case"]: r for r in committed.get("rows", [])}
    shared = [r for r in before if r["case"] in by_case]
    against_committed = _compare([by_case[r["case"]] for r in shared], shared, "committed", "offline_before")

    print(f"the case, network={'on' if net else 'off'}")
    reader(True)
    with tempfile.TemporaryDirectory(prefix="genomeos-vep-appends-") as appends:
        row = bench.run_case(case, net, sys.stdout, Path(appends))

    out = {
        "result": RESULT,
        "date": "2026-10-10",
        "lane": "lane-svcase",
        "registration": REGISTRATION,
        "roadmap_item": 'area F, "a scored structural-variant case"',
        "row": row,
        "structural_variant_record": [a.to_dict() for a in alt.read_sv_table(str(FIXTURES / case["sv"]))],
        "regression_guard": {
            "rows_compared": len(before),
            "fields_that_differ": differ,
            "fields_not_compared": sorted(NOT_COMPARED),
            "against_committed": against_committed,
            "rows_against_committed": len(shared),
            "verdict": (
                "the ten rows are identical in every field but seconds"
                if not differ
                else f"{len(differ)} field differences in the ten rows"
            ),
            "committed_verdict": (
                "the offline before-run reproduces the committed rows field for field"
                if not against_committed
                else f"{len(against_committed)} differences against the committed rows"
            ),
        },
        "note": (
            "The eleventh benchmark case, registered on 2026-10-10 before it was run "
            f"({REGISTRATION}) and added to the benchmark's case list whatever the outcome. Three "
            "quantities are scored apart, as the benchmark scores them: whether the target is "
            "recovered, whether the route is called correctly, and whether the mechanism ranked "
            "first is defensible. The ten existing rows are the regression guard and are re-scored "
            "offline with the structural-variant reader as it was before this lane's change and as "
            "it is now. A recovered target is a retrospective consistency with what is already "
            "approved for the alteration, never evidence that a therapy works."
        ),
    }
    out["result_manifest"] = _manifest(net, bench)
    print(save_result(RESULT, out))
    print(
        json.dumps(
            {
                k: row.get(k)
                for k in (
                    "gene",
                    "recovered",
                    "rank",
                    "score",
                    "surface_accessibility",
                    "best_mechanism",
                    "mechanism_sane",
                    "verdict",
                )
            },
            indent=1,
        )
    )
    return 0 if row.get("pass") and not differ else 1


def main(argv: list[str]) -> int:
    if len(argv) < 2 or argv[1] not in ("register", "run"):
        print(__doc__)
        return 2
    if argv[1] == "register":
        print(save_result(REGISTRATION, registration(_bench())))
        return 0
    return run("--offline" not in argv)


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
