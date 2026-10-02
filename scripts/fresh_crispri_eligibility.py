# SPDX-License-Identifier: AGPL-3.0-or-later
"""Eligibility of a fresh, non-K562 CRISPRi enhancer-gene set, at 0 model requests.

Authorised by Albert, 2026-10-02: "Use the freed slot for a 0-request search for a fresh non-K562
CRISPRi set (>=20 loci, >=30 positives, never read here)."

The floors are `genomeos/attribution/fresh.py`'s, and the locus convention is
`genomeos/attribution/cell2.py`'s, imported rather than restated. Both were fixed before any
candidate table was opened and neither is moved here.

What this script reads from a candidate's published table is restricted in code: every column name
passes `fresh.permit` before a workbook is opened, so an effect size, a fold change, a p-value, an
FDR or a model score is refused rather than loaded. Counting positives from a published hit-call
column is the lane's purpose; reading the magnitude of an effect is not permitted.

No AlphaGenome or other model request is made, nothing is scored, no money is spent.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import openpyxl

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import cell2, fresh  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "fresh_crispri_eligibility"
ROOT = Path(__file__).resolve().parents[1]
CACHE = Path("data/cache/fresh")

#: The files this lane wrote, so the stamp can say which uncommitted code is somebody else's. Three
#: other lanes are live in this one checkout, so a dirty stamp may be no act of this lane; a referee
#: can accept a named foreign file that the counting path does not import, but not a bare dirty flag.
#: This script, as the entry whose transitive import closure is the counting path of its results
#: (genomeos.manifest.counting_path, the one implementation). Named rather than derived from
#: __file__ so the closure is the same however the script is invoked.
ENTRY = "scripts/fresh_crispri_eligibility.py"

OWN_CODE = (
    "genomeos/attribution/fresh.py",
    "scripts/fresh_crispri_eligibility.py",
    "tests/test_fresh_crispri.py",
)

#: The candidate that cleared the floors, with its cached supplementary tables. Downloaded from the
#: publisher's open-access store, no login and no payment; bytes and sha256 recorded in the manifest.
ASTRO = {
    "doi": "10.1038/s41593-025-02154-3",
    "pmid": "41413662",
    "pmcid": "PMC12971484",
    "citation": (
        "Green NFO, Sutton GJ, Perez-Burillo J, et al. CRISPRi screening in cultured human "
        "astrocytes uncovers distal enhancers controlling genes dysregulated in Alzheimer's "
        "disease. Nature Neuroscience 2025;29(3):703-716"
    ),
    "files": {
        "supplementary_table_1": {
            "path": CACHE / "astro_MOESM3.xlsx",
            "url": "https://static-content.springer.com/esm/art%3A10.1038%2Fs41593-025-02154-3"
            "/MediaObjects/41593_2025_2154_MOESM3_ESM.xlsx",
            "sheet": "1B_CandidatePeakAnnotation",
        },
        "supplementary_table_3": {
            "path": CACHE / "astro_MOESM5.xlsx",
            "url": "https://static-content.springer.com/esm/art%3A10.1038%2Fs41593-025-02154-3"
            "/MediaObjects/41593_2025_2154_MOESM5_ESM.xlsx",
            "sheet": "3A_EnhancerGenePairDE",
            "power_sheet": "3F_PowerCalculations",
        },
    },
}


def read_sheet(path: Path, sheet: str, columns: list[str]) -> list[dict]:
    """Rows of `sheet`, restricted to `columns`, every one of which must pass `fresh.permit` first.

    Wholly blank trailing rows are counted and reported rather than dropped in silence: a row quietly
    skipped would change a count without changing what the count claims.
    """
    fresh.permit(columns)
    workbook = openpyxl.load_workbook(path, read_only=True)
    try:
        worksheet = workbook[sheet]
        rows = worksheet.iter_rows(values_only=True)
        header = list(next(rows))
        missing = [c for c in columns if c not in header]
        if missing:
            raise KeyError(f"{path.name}:{sheet} has no column {missing}")
        index = {c: header.index(c) for c in columns}
        out, blank = [], 0
        for row in rows:
            if not any(cell is not None for cell in row):
                blank += 1
                continue
            if len(row) <= max(index.values()):
                raise ValueError(f"{path.name}:{sheet} has a short non-blank row: {row!r}")
            out.append({c: row[index[c]] for c in columns})
    finally:
        workbook.close()
    print(f"  {path.name}:{sheet}: {len(out)} rows, {blank} wholly blank rows skipped")
    return out


def _true(value) -> bool:
    return str(value).strip().lower() == "true"


def astrocyte_counts() -> dict:
    """The counts of the astrocyte screen, from labels, coordinates and the published power flags."""
    table3 = ASTRO["files"]["supplementary_table_3"]
    pairs = read_sheet(
        table3["path"], table3["sheet"], ["Pair", "Enhancer", "GeneSymbol", "Hit", "EnhancerCoord"]
    )
    power = read_sheet(
        table3["path"],
        table3["power_sheet"],
        ["Pair", "Hit", "WellPowered_at_FC_0.15", "WellPowered_at_FC_0.25"],
    )
    table1 = ASTRO["files"]["supplementary_table_1"]
    elements = read_sheet(table1["path"], table1["sheet"], ["Enh", "Coord", "Tested", "Hit"])

    positives = [r for r in pairs if _true(r["Hit"])]
    keys = fresh.locus_keys("astrocyte", positives, "EnhancerCoord", "GeneSymbol")
    loci = cell2.count_loci(keys)
    by_chrom: dict[str, int] = {}
    for key in keys:
        by_chrom[key.chrom] = by_chrom.get(key.chrom, 0) + 1

    # the two tables must agree on the number of positives, or one of them is being misread
    if sum(1 for r in power if _true(r["Hit"])) != len(positives):
        raise ValueError("the pair table and the power table disagree on the number of positives")

    return {
        "candidate_elements": len(elements),
        "elements_tested": sum(1 for r in elements if _true(r["Tested"])),
        "elements_called_hit": sum(1 for r in elements if _true(r["Hit"])),
        "pairs_tested": len(pairs),
        "measured_positives": len(positives),
        "distinct_measured_genes": len({k.gene for k in keys}),
        "chromosomes": len(by_chrom),
        "positives_per_chromosome": dict(sorted(by_chrom.items())),
        "independent_loci": loci,
        "positives_per_locus_max": max(
            sum(1 for g in cell2.group(keys) if g == gid) for gid in set(cell2.group(keys))
        ),
        "well_powered_pairs_at_fc_0_15": sum(1 for r in power if _true(r["WellPowered_at_FC_0.15"])),
        "well_powered_pairs_at_fc_0_25": sum(1 for r in power if _true(r["WellPowered_at_FC_0.25"])),
        "power_is_published": (
            "the authors publish a per-pair power flag (Supplementary Table 3F, WellPowered at fold "
            "change 0.15 and 0.25). This is the quantity the second-cell-type lane could not compute "
            "for the benchmark's held-out cell types; it is read here as a label, not as a figure of "
            "merit, and no power calculation is performed by this lane"
        ),
    }


def candidates(astro: dict) -> fresh.Verdict:
    """Every candidate examined, each with how its freshness was established and which floor it misses.

    The counts for the astrocyte screen are computed above from its published tables. The counts for
    the other candidates are what their publications state; where a publication does not state a
    count, it is left as None and the floor reads as unchecked rather than cleared.
    """
    return fresh.Verdict(
        candidates=[
            fresh.Candidate(
                name="Green 2025 astrocyte CRISPRi (AstroREG)",
                cell_context="cultured primary human astrocytes, fetal-derived (Lonza CC-2565), "
                "two donor lines",
                assay="CRISPRi (dCas9-KRAB) with single-cell RNA-seq readout (CROP-seq)",
                accession=f"doi {ASTRO['doi']}, PMID {ASTRO['pmid']}, {ASTRO['pmcid']}; "
                "Supplementary Tables 1 and 3",
                pairs_tested=astro["pairs_tested"],
                positives=astro["measured_positives"],
                independent_loci=astro["independent_loci"],
                distinct_genes=astro["distinct_measured_genes"],
                chromosomes=astro["chromosomes"],
                read_here="no outcome from this set has been read, scored or quoted here. "
                "Established by grep over docs/, genomeos/, scripts/, tests/ and the knowledge "
                "READMEs for Voineagu, AstroREG, PsychENCODE, EGrf, CROP-seq, NHA, s41593, 02154 "
                "and PMID 41413662: no hit for any of them, and the study is absent from the "
                "project's CRISPRi exposure ledger. data/knowledge/crispri/ holds only the two "
                "EPCrisprBenchmark files",
                access="open access, CC BY 4.0; supplementary tables served by the publisher over "
                "https with no login and no payment",
                notes="the project has read ENCODE astrocyte chromatin peaks in its node reader "
                "lanes, so astrocyte *chromatin features* are development-exposed even though no "
                "outcome of this CRISPRi set is. Supplementary Table 1B annotates its elements "
                "against K562 validated enhancers (Yao 2022), a column this lane did not read",
            ),
            fresh.Candidate(
                name="Wang 2026 primary human T cell CRISPR enhancer annotation",
                cell_context="primary human T cells, tumour-infiltrating CD8 T cells, CAR T cells",
                assay="CRISPRi (dCas9-KRAB) and Cas9-indel mutagenesis",
                accession="bioRxiv 2026, licence forbids PMC archiving; no accession established",
                pairs_tested=None,
                positives=None,
                independent_loci=3,
                read_here="no outcome read here; 'primary T cells' appears once in the record, in a "
                "related-work overlap audit of Pacalin 2024, with no outcome read",
                access="preprint; the licence does not permit PMC archiving, so the full text and "
                "its supplementary tables were not retrieved at 0 cost",
                notes="the screen is built around three genes (PDCD1, HAVCR2, TBX21), so under the "
                "registered convention its positives cannot span more than 3 loci however many "
                "elements were tiled. Pair and positive counts were not established",
            ),
            fresh.Candidate(
                name="Caragine and Le 2025, 2.8 Mb TAD dissection in six cancer lines",
                cell_context="six human cancer cell lines",
                assay="CRISPRi tiling of one topologically associating domain",
                accession="Nature Communications 2025, s41467-025-56568-5",
                pairs_tested=None,
                positives=None,
                independent_loci=1,
                read_here="no outcome read here; the study is absent from the exposure ledger",
                access="open access",
                notes="a contiguous 2.8 Mb tiling on one chromosome chains into a single locus under "
                "cell2.INDEPENDENT_LOCUS_RULE, because each element lies within 1 Mb of the next. "
                "Dense tiling of one neighbourhood cannot clear a locus floor by construction",
            ),
            fresh.Candidate(
                name="EPCrisprBenchmark held-out five cell types (Gschwind 2025)",
                cell_context="GM12878, HCT116, Jurkat, WTC11 (and K562 training)",
                assay="CRISPRi, pooled from nine upstream screens",
                accession="EngreitzLab/CRISPR_comparison, heldout_5_cell_types.GRCh38",
                pairs_tested=4378,
                positives=190,
                independent_loci=22,
                read_here="yes: effect sizes, significance calls, per-screen prevalences, AUPRC and "
                "AlphaGenome deletion values have all been read and scored, in all five cell types. "
                "It is this project's primary registered endpoint",
                access="open",
                notes="development evidence, not fresh. Recorded here so the table states why the "
                "set this project already holds cannot answer the question",
            ),
            fresh.Candidate(
                name="IGVF K562 MHC CRISPRi Perturb-seq",
                cell_context="K562",
                assay="CRISPRi (dCas9-KRAB) Perturb-seq, targeted panel",
                accession="IGVFDS7132YVKO",
                pairs_tested=None,
                positives=None,
                independent_loci=None,
                read_here="yes: scored as this project's independent benchmark, "
                "data/results/indep_mhc_crispri.json",
                access="open, CC BY 4.0",
                notes="K562, and already scored. Listed for completeness",
            ),
        ]
    )


def manifest(astro: dict) -> dict:
    sources = []
    inputs = []
    for label, spec in ASTRO["files"].items():
        sources.append(
            {
                "accession": f"{ASTRO['citation']} -- {label.replace('_', ' ')}",
                "version": "publisher open-access supplementary store; pinned here by sha256",
                "url": spec["url"],
            }
        )
        inputs.append(
            mf.input_entry(
                spec["path"],
                partition=None,
                role="labels, element coordinates, measured gene symbols and the published "
                "per-pair power flags only; no effect size, p-value, FDR or model score was read",
            )
        )
    return {
        "sources": sources,
        "inputs": inputs,
        "assembly": "GRCh38",
        "coordinates": (
            "n/a: the publication states no coordinate convention for its chrom:start-end element "
            "strings, so this lane does not assert one. The only use made of them is the 1 Mb "
            "locus grouping, whose outcome a one-base shift cannot change; no coordinate from "
            "this set is joined to a project interval here, and any later join must settle the "
            "convention first"
        ),
        "parameters": {
            "independent_locus_span_bp": cell2.INDEPENDENT_LOCUS_SPAN,
            "locus_floor": fresh.LOCUS_FLOOR,
            "positive_floor": fresh.POSITIVE_FLOOR,
            "excluded_cell": fresh.EXCLUDED_CELL,
            "permitted_columns": sorted(fresh.PERMITTED_COLUMNS),
            "alphagenome_requests": 0,
            "model_requests": 0,
            "performance_metrics_computed": 0,
            "gains_computed": 0,
            "money_spent": 0,
        },
        "exclusions": [
            "no effect size, log fold change, p-value, FDR, Z statistic, gene-expression level or "
            "sensitivity figure was read from any candidate table: fresh.permit refuses the column "
            "before the workbook is opened",
            "K562 is not a candidate cell context: the committed result already covers it",
            "a candidate whose outcome has been read here is development evidence and is recorded as "
            "not fresh rather than counted",
            "mouse screens are not candidates: they share no hg38 element with this project",
            "no candidate requiring a login, a payment or controlled access was retrieved",
        ],
        "partitions": {
            "astrocyte_pairs": f"all {astro['pairs_tested']} element-gene pairs of the published "
            f"screen; {astro['measured_positives']} carry the authors' hit call"
        },
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


def main() -> int:
    print("Reading candidate tables (labels, coordinates and power flags only)")
    astro = astrocyte_counts()
    verdict = candidates(astro)
    payload = verdict.to_dict()
    payload["astrocyte_counts"] = astro
    payload["authorisation"] = (
        "Albert, 2026-10-02: use the freed slot for a 0-request search for a fresh non-K562 CRISPRi "
        "set (>=20 loci, >=30 positives, never read here)"
    )
    payload["what_a_test_here_could_claim"] = [
        "a first reading of the frozen model's enhancer-gene separation in a cell context this "
        "project has never scored, on a set whose outcomes it has never read",
        f"an interval over {astro['independent_loci']} independent loci, which is at or above the "
        "ten-cluster minimum, so a resampling interval is reportable rather than refused as "
        "'interval unreliable'",
        "a power statement taken from the authors' own published per-pair flag rather than computed "
        "here, which is what the second-cell-type lane could not obtain",
    ]
    payload["what_a_test_here_could_not_claim"] = [
        "not a replication of the K562 result: a different assay readout (single-cell CROP-seq "
        "against the benchmark's mixed FlowFISH and scRNA-seq screens), a different element "
        "selection (PsychENCODE astrocyte enhancers), a different hit threshold and a different "
        "cell lineage, so a difference cannot be attributed to the cell type",
        "not independence of the genomic regions: the locus grouping is cell2's operational "
        "convention, and AlphaGenome was trained on ENCODE tracks genome-wide, which include "
        "astrocyte chromatin; the labels are unseen, the regions are not",
        "not a cell-type-specific claim about astrocytes: two donor lines of cultured fetal-derived "
        "astrocytes are not primary brain tissue",
        "no claim at all until the comparison is registered before any score is computed",
    ]
    payload["left_undone"] = [
        "no comparison registered and no score computed: this lane is eligibility only",
        "the astrocyte element set has not been intersected with the project's hg38 element "
        "inventory, so the overlap with already-scored elements is not quantified",
        "the genome build of the astrocyte coordinates is inferred from the publication's use of "
        "GENCODE v32 and PsychENCODE enhancers and has not been confirmed against a stated build",
        "Supplementary Table 3's 158 pairs and Table 1B's 145 hit elements differ because 1B carries "
        "one representative gene per element; the pair table is authoritative and is what was counted",
    ]
    path = save_result(RESULT, payload, manifest=manifest(astro))
    print(json.dumps({k: payload[k] for k in ("verdict", "eligible_sets", "candidates_examined")}))
    print(
        f"astrocyte: {astro['measured_positives']} positives, {astro['independent_loci']} loci, "
        f"{astro['distinct_measured_genes']} genes, {astro['chromosomes']} chromosomes "
        f"(floors: {fresh.POSITIVE_FLOOR} positives, {fresh.LOCUS_FLOOR} loci)"
    )
    print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
