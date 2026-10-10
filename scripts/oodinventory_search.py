# SPDX-License-Identifier: AGPL-3.0-or-later
"""The search the rule in `data/results/oodinventory_rule.json` was written for: which public
cell types have BOTH open accessibility and open RNA-seq from one source, and are ABSENT from
the AlphaGenome RNA track metadata we hold.

    uv run python scripts/oodinventory_search.py

ADDS NO BIOSAMPLE, runs no out-of-domain test, makes 0 model requests and spends nothing
(2026-09-27 stop list). It reads catalogue metadata over HTTP and downloads no experimental data.

WHERE IT LOOKS, and why that is the whole of it. The rule requires accessibility AND RNA-seq
from the SAME source, open, with no application. That requirement is what narrows the field, not
a shortage of effort:

  - ENCODE (www.encodeproject.org) is searched exhaustively through its public JSON API: every
    released human experiment of each accessibility assay (DNase-seq, ATAC-seq, snATAC-seq) and
    each RNA-seq assay (total, polyA plus, polyA minus, long read, snRNA-seq, scRNA-seq),
    grouped by `biosample_ontology.term_id`. Both modalities then come from one archive under one
    pipeline, which is what "same source" means. The perturbation RNA-seq assays (shRNA, siRNA,
    CRISPR, CRISPRi) are excluded: their biosamples are genetically modified, so they are not a
    cell type's own transcriptome. CAGE, RAMPAGE, PRO-cap, GRO-cap, Bru-seq, microRNA-seq, small
    RNA-seq, RNA-PET, PAS-seq and RNA microarray are excluded because they are not RNA-seq. GM
    DNase-seq, FAIRE-seq and MNase-seq are excluded because the rule names DNase and ATAC.
  - GTEx supplies 54 of the metadata's own RNA rows and publishes no matched accessibility, so it
    cannot supply a pair. Its individual-level data is dbGaP-controlled, which the rule's licence
    clause excludes outright.
  - HCA and Tabula Sapiens publish single-cell RNA openly and no matched accessibility of their
    own; where a project does pair scRNA with scATAC the pair is usually in a controlled tier.
  - Roadmap Epigenomics is folded into ENCODE's portal and is reached by the same search.
  These three are recorded as searched-and-excluded in the result, with the clause each fails,
  rather than left unmentioned.

HOW THE ONTOLOGY CLAUSE IS APPLIED. The rule says a candidate is absent only when no covered
CURIE is equal to it, an ancestor of it, or a descendant of it, and that the relation is
DECLARED with evidence rather than assumed. Here the declaration comes from EBI OLS4
(www.ebi.ac.uk/ols4), which serves CL, UBERON, EFO and CLO: for each candidate we fetch
`hierarchicalAncestors` and, unless the term is a donor cell line with no cell-type children,
`hierarchicalDescendants`, and any covered CURIE in either list is declared with OLS named as
the evidence. Two kinds of candidate OLS cannot answer are declared by hand, each with its
reason, in `DECLARED_BY_HAND`: ENCODE's own NTR terms, which exist in no public ontology, and
one term whose id migrated (`CL:4006001` is the covered `NTR:0000524` renamed, same name, same
cells).

WHAT THE RULE'S CLAUSE DOES NOT CATCH, stated because the result has to say it. The clause is
is_a/part_of: it excludes a parent and a child, and it does NOT exclude a SIBLING. So a
lymphoblastoid line from a different donor than the covered GM12878 reads as absent although it
is the same cell type, and a cancer cell line reads as absent although a sibling line of the
same tissue is covered. Those candidates are flagged `sibling_of_covered` and counted
separately. The rule is applied as committed and NOT amended: the flag is a disclosure, so that
whoever reads the list can see which rows the rule admits for a reason it should not have.
"""

from __future__ import annotations

import collections
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import oodinventory as ood  # noqa: E402
from genomeos.results import save_result  # noqa: E402

NAME = "oodinventory_candidates"
ENTRY = "scripts/oodinventory_search.py"
OWN_CODE = (
    "genomeos/attribution/oodinventory.py",
    "scripts/oodinventory_curies.py",
    "scripts/oodinventory_search.py",
    "tests/test_oodinventory.py",
)

ENCODE = "https://www.encodeproject.org"
OLS = "https://www.ebi.ac.uk/ols4/api"
ENCODE_LICENCE = (
    "ENCODE open access: released data are freely available with no application and no "
    "data-use agreement (encodeproject.org/help/citing-encode/)"
)

ACCESSIBILITY_ASSAYS = ("DNase-seq", "ATAC-seq", "snATAC-seq")
RNA_ASSAYS = (
    "total RNA-seq",
    "polyA plus RNA-seq",
    "polyA minus RNA-seq",
    "long read RNA-seq",
    "snRNA-seq",
    "scRNA-seq",
)
EXCLUDED_ASSAYS = {
    "shRNA RNA-seq": "genetically modified biosample: a knockdown, not the cell type's own RNA",
    "siRNA RNA-seq": "genetically modified biosample: a knockdown, not the cell type's own RNA",
    "CRISPR RNA-seq": "genetically modified biosample",
    "CRISPRi RNA-seq": "genetically modified biosample",
    "microRNA-seq": "not RNA-seq: a different RNA population",
    "small RNA-seq": "not RNA-seq: a different RNA population",
    "microRNA counts": "not RNA-seq",
    "CAGE": "not RNA-seq: 5' ends",
    "RAMPAGE": "not RNA-seq: 5' ends",
    "PRO-cap": "not RNA-seq: nascent 5' ends",
    "GRO-cap": "not RNA-seq: nascent 5' ends",
    "PRO-seq": "not RNA-seq: nascent transcription",
    "Bru-seq": "not RNA-seq: metabolic labelling",
    "BruChase-seq": "not RNA-seq: metabolic labelling",
    "BruUV-seq": "not RNA-seq: metabolic labelling",
    "RNA-PET": "not RNA-seq",
    "PAS-seq": "not RNA-seq: 3' ends",
    "RNA microarray": "not sequencing",
    "GM DNase-seq": "genetically modified biosample",
    "FAIRE-seq": "the rule names DNase and ATAC",
    "MNase-seq": "the rule names DNase and ATAC",
}

#: Sources searched and excluded by a named clause of the rule, recorded rather than left out.
OTHER_SOURCES = [
    {
        "source": "GTEx",
        "clause_failed": "same source, and licence",
        "detail": "GTEx supplies 54 of the metadata's own rna_seq rows and publishes no matched "
        "chromatin accessibility of its own, so it cannot supply a pair at all. Its "
        "individual-level data is dbGaP-controlled, which the rule's licence clause excludes.",
    },
    {
        "source": "Human Cell Atlas (HCA) Data Portal",
        "clause_failed": "accessibility, and same source",
        "detail": "HCA's open releases are single-cell RNA. Where a project pairs scRNA with "
        "scATAC the matched pair is in a managed-access tier, which the licence clause excludes; "
        "pairing an open HCA RNA project with an accessibility dataset from another archive fails "
        "the same-source clause.",
    },
    {
        "source": "Tabula Sapiens",
        "clause_failed": "accessibility",
        "detail": "Tabula Sapiens publishes single-cell RNA across tissues under CC BY 4.0 and no "
        "chromatin accessibility, so no pair exists within the source.",
    },
    {
        "source": "Roadmap Epigenomics",
        "clause_failed": "none: it is inside the ENCODE search",
        "detail": "Roadmap's DNase and RNA datasets are served from the ENCODE portal and are "
        "already covered by the exhaustive ENCODE query below, so they are not a separate source.",
    },
    {
        "source": "GEO / SRA",
        "clause_failed": "same source",
        "detail": "GEO is an archive of independent submissions with no common pipeline, so a "
        "DNase series and an RNA series in GEO are not one source under the rule unless one "
        "submission carries both; any such submission that also appears on the ENCODE portal is "
        "already in the ENCODE search.",
    },
]

#: Relations OLS cannot serve, declared by hand with the evidence, as the rule's clause requires.
DECLARED_BY_HAND: dict[str, tuple[tuple[str, ...], str]] = {
    # ENCODE's NTR terms are portal-local and exist in no public ontology, so OLS has no entry
    # for them. Each NTR term name states its parent cell type, and that parent is covered.
    "NTR:0000497": (
        ("CL:0000895",),
        "ENCODE NTR term 'activated naive CD4-positive, alpha-beta T cell': an activated "
        "state of the covered CL:0000895 naive thymus-derived CD4-positive alpha-beta T cell",
    ),
    "NTR:0000499": (
        ("CL:0000897",),
        "'activated CD4-positive, alpha-beta memory T cell': an activated state of the covered CL:0000897",
    ),
    "NTR:0000502": (
        ("CL:0000900",),
        "'activated naive CD8-positive, alpha-beta T cell': an activated state of the covered CL:0000900",
    ),
    "NTR:0000504": (
        ("CL:0000909",),
        "'activated CD8-positive, alpha-beta memory T cell': an activated state of the covered CL:0000909",
    ),
    "NTR:0000506": (("CL:0000236",), "'activated B cell': a state of the covered CL:0000236 B cell"),
    "NTR:0000537": (("CL:0000084",), "'activated T-cell': a state of the covered CL:0000084 T-cell"),
    "NTR:0000633": (
        ("CL:0000084",),
        "'activated T-helper 1 cell': a CD4 T-cell subset, inside the covered "
        "CL:0000084 T-cell and CL:0000624",
    ),
    "NTR:0000634": (
        ("CL:0000084",),
        "'activated T-helper 2 cell': a CD4 T-cell subset, inside the covered "
        "CL:0000084 T-cell and CL:0000624",
    ),
    "NTR:0000635": (
        ("CL:0000084",),
        "'activated T-helper 17 cell': a CD4 T-cell subset, inside the covered "
        "CL:0000084 T-cell and CL:0000624",
    ),
    "NTR:0000660": (
        ("CL:0000897",),
        "'activated effector memory CD4-positive, alpha-beta T cell': inside the "
        "covered CL:0000897 CD4-positive alpha-beta memory T cell",
    ),
    "NTR:0000663": (
        ("CL:0000624",),
        "'stimulated activated CD4-positive, alpha-beta T cell': a state of the covered CL:0000624",
    ),
    "NTR:0000672": (("CL:0000084",), "'activated gamma-delta T cell': inside the covered CL:0000084 T-cell"),
    # An id migration, not a new cell type: the covered NTR:0000524 carries the identical name.
    "CL:4006001": (
        ("NTR:0000524",),
        "the same term under a migrated id: the covered set carries NTR:0000524 "
        "'fibroblast of skin of scalp', name-identical. CURIE matching alone would have "
        "called this absent, which is why the rule declares relations.",
    ),
}

#: A candidate whose rule-absence rests on the clause not excluding siblings, with the sibling.
SIBLING_NOTE = "sibling_of_covered"


def _get(url: str, cache: Path | None, key: str, pause: float = 0.0) -> dict[str, Any]:
    """One JSON GET, optionally cached on disk so a re-run is cheap and reproducible."""
    if cache is not None:
        hit = cache / f"{key}.json"
        if hit.exists():
            return json.loads(hit.read_text())
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=120) as handle:
            payload = json.load(handle)
    except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as exc:
        payload = {"_error": f"{type(exc).__name__}: {exc}"}
    if cache is not None:
        cache.mkdir(parents=True, exist_ok=True)
        (cache / f"{key}.json").write_text(json.dumps(payload))
    if pause:
        time.sleep(pause)
    return payload


def encode_query(assay: str) -> str:
    return (
        f"{ENCODE}/search/?type=Experiment&status=released"
        "&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens"
        f"&assay_title={urllib.parse.quote_plus(assay)}"
        "&field=accession&field=assay_title&field=biosample_ontology.term_id"
        "&field=biosample_ontology.term_name&field=biosample_ontology.classification"
        "&format=json&limit=all"
    )


def ols_url(curie: str, relation: str) -> str | None:
    prefix = curie.split(":")[0]
    onto = {"CL": "cl", "UBERON": "uberon", "EFO": "efo", "CLO": "clo"}.get(prefix)
    if onto is None:
        return None
    iri = (
        f"http://www.ebi.ac.uk/efo/{curie.replace(':', '_')}"
        if prefix == "EFO"
        else f"http://purl.obolibrary.org/obo/{curie.replace(':', '_')}"
    )
    quoted = urllib.parse.quote(urllib.parse.quote(iri, safe=""), safe="")
    return f"{OLS}/ontologies/{onto}/terms/{quoted}/{relation}?size=1000"


def ols_terms(curie: str, relation: str, cache: Path | None) -> tuple[list[str], str]:
    url = ols_url(curie, relation)
    if url is None:
        return [], f"no OLS ontology for the {curie.split(':')[0]} prefix"
    payload = _get(url, cache, f"ols_{relation}_{curie.replace(':', '_')}", pause=0.05)
    if "_error" in payload:
        return [], payload["_error"]
    terms = [t["obo_id"] for t in payload.get("_embedded", {}).get("terms", []) if t.get("obo_id")]
    return terms, "ok"


def _donor_line(name: str | None, classification: str | None) -> bool:
    """A named donor cell line (GM12878-style): a leaf term with no cell-type children."""
    return bool(classification == "cell line" and name and re.fullmatch(r"(GM|HG)\d+", name))


def gather_encode(cache: Path | None) -> dict[str, Any]:
    acc: dict[str, list[str]] = collections.defaultdict(list)
    rna: dict[str, list[str]] = collections.defaultdict(list)
    names: dict[str, str] = {}
    classes: dict[str, str] = {}
    totals: dict[str, int] = {}
    for assay in ACCESSIBILITY_ASSAYS + RNA_ASSAYS:
        payload = _get(encode_query(assay), cache, f"encode_{assay.replace(' ', '_')}")
        graph = payload.get("@graph", [])
        totals[assay] = len(graph)
        for row in graph:
            onto = row.get("biosample_ontology") or {}
            term = onto.get("term_id")
            if not term:
                continue
            names[term] = onto.get("term_name") or ""
            classes[term] = onto.get("classification") or ""
            bucket = acc if assay in ACCESSIBILITY_ASSAYS else rna
            bucket[term].append(f"{row['accession']} ({assay})")
    return {
        "accessibility": dict(acc),
        "rna": dict(rna),
        "names": names,
        "classes": classes,
        "experiment_counts": totals,
    }


def build_candidates(
    harvest: dict[str, Any], covered: set[str], cache: Path | None
) -> tuple[list[ood.Candidate], dict[str, Any]]:
    acc, rna = harvest["accessibility"], harvest["rna"]
    names, classes = harvest["names"], harvest["classes"]
    both = sorted(set(acc) & set(rna))
    candidates: list[ood.Candidate] = []
    ols_problems: list[str] = []
    for term in both:
        hand = DECLARED_BY_HAND.get(term)
        related: list[str] = []
        evidence: list[str] = []
        if term in covered:
            related.append(term)
            evidence.append("the CURIE is itself in the covered set")
        if hand is not None:
            related.extend(c for c in hand[0] if c in covered)
            evidence.append(hand[1])
        ancestors, up_ok = ols_terms(term, "hierarchicalAncestors", cache)
        if up_ok != "ok":
            ols_problems.append(f"{term} ancestors: {up_ok}")
        hits = [c for c in ancestors if c in covered]
        if hits:
            related.extend(hits)
            evidence.append(f"EBI OLS4 hierarchicalAncestors of {term} include {sorted(hits)}")
        down_ok = "not queried"
        if not _donor_line(names.get(term), classes.get(term)):
            descendants, down_ok = ols_terms(term, "hierarchicalDescendants", cache)
            if down_ok != "ok":
                ols_problems.append(f"{term} descendants: {down_ok}")
            down = [c for c in descendants if c in covered]
            if down:
                related.extend(down)
                evidence.append(f"EBI OLS4 hierarchicalDescendants of {term} include {sorted(down)}")
        else:
            down_ok = "ok"
            evidence.append(
                "descendants not queried: a named donor cell line is a leaf term with no cell-type children"
            )
        # The relation counts as CHECKED only when both directions were answered, or a hand
        # declaration covers the term. An OLS failure leaves it relation_unchecked, which the
        # rule refuses to read as absence.
        checked = bool(hand) or term in covered or (up_ok == "ok" and down_ok == "ok")
        sibling = _sibling_flag(term, names.get(term, ""), classes.get(term, ""), covered)
        candidates.append(
            ood.Candidate(
                curie=term,
                label=names.get(term, ""),
                source="ENCODE",
                accessibility=", ".join(sorted(acc[term])[:6])
                + (f" (+{len(acc[term]) - 6} more)" if len(acc[term]) > 6 else ""),
                rna=", ".join(sorted(rna[term])[:6])
                + (f" (+{len(rna[term]) - 6} more)" if len(rna[term]) > 6 else ""),
                licence=ENCODE_LICENCE,
                size="processed files not downloaded; "
                f"{len(acc[term])} accessibility and {len(rna[term])} RNA experiments on the "
                "portal, each with released GRCh38 signal/peaks or gene quantifications",
                related_covered=tuple(dict.fromkeys(related)),
                relation_checked=checked,
                notes="; ".join(evidence),
                extra={
                    "classification": classes.get(term, ""),
                    "accessibility_experiment_count": str(len(acc[term])),
                    "rna_experiment_count": str(len(rna[term])),
                    "rule_clause_not_excluding": sibling or "",
                    "rule_clause_not_excluding_kind": _sibling_kind(
                        term, names.get(term, ""), classes.get(term, ""), covered
                    ),
                },
            )
        )
    return candidates, {"terms_with_both": len(both), "ols_problems": ols_problems}


#: The covered neighbour a named cell line's rule-absence rests on, so the flag is specific
#: rather than a generic sentence. Read off the covered set's own 52 cell-line terms.
LINE_NEIGHBOURS = {
    "EFO:0001247": "T47D is a breast ductal carcinoma line; the covered set holds MCF-7 "
    "(EFO:0001203), MCF 10A (EFO:0001200) and breast epithelium (UBERON:0008367)",
    "EFO:0002793": "HL-60 is an acute myeloid leukaemia line; the covered set holds K562 "
    "(EFO:0002067, chronic myeloid leukaemia) and CL:0001059 common myeloid progenitor",
    "EFO:0005718": "Ishikawa is an endometrial adenocarcinoma line; the covered set holds the "
    "uterus itself (UBERON:0000995) and uterine smooth muscle cell (CL:0002601) as tissue and "
    "primary-cell terms, but no endometrial line",
}


def _sibling_kind(term: str, name: str, classification: str, covered: set[str]) -> str:
    """Which category of not-excluded-by-the-clause this row falls in, '' when none."""
    if term in covered:
        return ""
    if re.fullmatch(r"(GM|HG)\d+", name or ""):
        return "donor_variant_of_a_covered_lymphoblastoid_line"
    if classification == "cell line":
        return "immortalised_or_cancer_line_with_covered_neighbours"
    return ""


def _sibling_flag(term: str, name: str, classification: str, covered: set[str]) -> str | None:
    """Whether this row's rule-absence rests on the clause not excluding siblings."""
    kind = _sibling_kind(term, name, classification, covered)
    if kind == "donor_variant_of_a_covered_lymphoblastoid_line":
        return (
            "the covered set holds GM12878 (EFO:0002784), GM12891 and GM12892, which are the "
            "SAME cell type from other donors; the rule's is_a clause does not exclude a sibling, "
            "so this row reads absent for a reason it should not have"
        )
    if kind == "immortalised_or_cancer_line_with_covered_neighbours":
        return LINE_NEIGHBOURS.get(
            term,
            "an immortalised or cancer line; the covered set holds 52 cell-line terms and the "
            "rule's is_a clause does not exclude a sibling line",
        )
    return None


def _shortlist_caveats(rows: list[dict[str, Any]]) -> list[dict[str, str]]:
    """What a reader of the shortlist would otherwise have to work out from the accessions.

    A row only passes the rule's modality clause; it does not promise the modality is the one a
    measured-RNA comparison can use. A term whose only RNA is single-nucleus needs a pseudobulk
    before any per-gene value exists, and that is a step, not a download.
    """
    out: list[dict[str, str]] = []
    for row in rows:
        if not row["absent_by_rule"]:
            continue
        rna = row["rna"] or ""
        if "RNA-seq" in rna and "snRNA-seq" not in rna and "scRNA-seq" not in rna:
            continue
        if "snRNA-seq" in rna or "scRNA-seq" in rna:
            out.append(
                {
                    "curie": row["curie"],
                    "label": row["label"],
                    "caveat": "every RNA experiment on this term is single-nucleus or "
                    "single-cell, so a per-gene value needs a pseudobulk first. The rule admits "
                    "it ('bulk or single-cell with a usable pseudobulk'); it is a step before any "
                    "comparison, not a file to download.",
                }
            )
    return out


def result_manifest(read_now: list[str], harvest: dict[str, Any]) -> dict[str, Any]:
    declared = [p for p in read_now if p.endswith("alphagenome_track_metadata_copy.csv")] or [
        str(ood.METADATA_PATH)
    ]
    inputs = [
        mf.input_entry(p, partition="not an evaluation partition: the catalogue being counted")
        for p in declared
    ]
    return {
        "sources": [
            {
                "accession": "AlphaGenome published output-track metadata, local copy at "
                "data/cache/entex/alphagenome_track_metadata_copy.csv",
                "version": "as held in this checkout. A copy of a published TRACK LIST, not of "
                "the training set.",
            },
            {
                "accession": "ENCODE portal search API, released human experiments of "
                + ", ".join(ACCESSIBILITY_ASSAYS + RNA_ASSAYS)
                + "; metadata only, no experimental file fetched",
                "version": f"queried {time.strftime('%Y-%m-%d')}; experiment counts per assay in "
                f"encode_experiment_counts; {ENCODE_LICENCE}",
            },
            {
                "accession": "EBI OLS4 hierarchicalAncestors and hierarchicalDescendants for CL, "
                "UBERON, EFO and CLO, for the rule's ontology clause",
                "version": f"queried {time.strftime('%Y-%m-%d')}; OLS4 serves the current release "
                "of each ontology and is not pinned here, which is recorded as a limitation",
            },
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "n/a: no coordinate is read. This result counts ontology terms and names "
        "experiment accessions.",
        "coordinates": "n/a: no interval is read or written",
        "parameters": {
            "accessibility_assays": list(ACCESSIBILITY_ASSAYS),
            "rna_assays": list(RNA_ASSAYS),
            "match_field": ood.MATCH_FIELD,
            "selection_rule_version": ood.SELECTION_RULE["version"],
            "encode_experiment_counts": harvest["experiment_counts"],
        },
        "exclusions": [f"{k}: {v}" for k, v in sorted(EXCLUDED_ASSAYS.items())]
        + [
            "every non-human experiment, and every experiment not status=released",
            "every experimental data file: this search reads catalogue metadata only",
            "no biosample is added and no out-of-domain test is run (2026-09-27 stop list)",
        ],
        "partitions": "n/a: nothing is fitted and nothing is held out; the output is a list",
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


def main() -> int:
    mf.trace_begin()
    started = time.time()
    cache_dir = os.environ.get("OODINVENTORY_CACHE")
    cache = Path(cache_dir) if cache_dir else None
    covered = ood.covered_curies(ROOT / ood.METADATA_PATH)
    harvest = gather_encode(cache)
    candidates, search = build_candidates(harvest, covered, cache)
    record: dict[str, Any] = dict(ood.inventory(candidates, covered))
    rows = record["candidates"]
    absent = [r for r in rows if r["absent_by_rule"]]
    siblings = [r for r in absent if r.get("rule_clause_not_excluding")]
    record.update(
        {
            "lane": "lane-oodinventory",
            "stage": "the search the committed rule was written for",
            "rule_result": "data/results/oodinventory_rule.json",
            "encode_terms_with_accessibility": len(harvest["accessibility"]),
            "encode_terms_with_rna": len(harvest["rna"]),
            "encode_terms_with_both": search["terms_with_both"],
            "absent_by_rule_count": len(absent),
            "absent_and_sibling_of_covered_count": len(siblings),
            "absent_breakdown": dict(
                collections.Counter(
                    r.get("rule_clause_not_excluding_kind") or "no_covered_sibling" for r in absent
                )
            ),
            "absent_after_the_sibling_disclosure": [
                r["curie"] for r in absent if not r.get("rule_clause_not_excluding")
            ],
            "shortlist": [
                {
                    "curie": r["curie"],
                    "label": r["label"],
                    "classification": r["classification"],
                    "source": r["source"],
                    "accessibility": r["accessibility"],
                    "accessibility_experiment_count": r["accessibility_experiment_count"],
                    "rna": r["rna"],
                    "rna_experiment_count": r["rna_experiment_count"],
                    "licence": r["licence"],
                    "size": r["size"],
                    "absent_by_rule": r["absent_by_rule"],
                    "rule_clause_not_excluding": r["rule_clause_not_excluding"],
                    "rule_clause_not_excluding_kind": r["rule_clause_not_excluding_kind"],
                }
                for r in absent
                if r.get("rule_clause_not_excluding_kind") != "donor_variant_of_a_covered_lymphoblastoid_line"
            ],
            "shortlist_note": (
                "The absent rows that are NOT a different donor of a covered lymphoblastoid "
                "line. This is the list to read: the donor variants are absent only because the "
                "rule's clause does not exclude a sibling, and they are the same cell type as "
                "the covered GM12878."
            ),
            "shortlist_caveats": _shortlist_caveats(rows),
            "other_sources_searched_and_excluded": OTHER_SOURCES,
            "excluded_assays": EXCLUDED_ASSAYS,
            "ols_problems": search["ols_problems"],
            "rule_clause_limitation": (
                "The rule's ontology clause is is_a/part_of: it excludes a parent and a child and "
                "NOT a sibling. Rows flagged rule_clause_not_excluding are absent under the rule "
                "as committed although a covered term is the same kind of biosample. The rule is "
                "applied as committed and not amended; the flag is the disclosure."
            ),
            "wording_limit": ood.SELECTION_RULE["wording_limit"],
            "no_biosample_added": ood.SELECTION_RULE["no_biosample"],
            "alphagenome_requests": 0,
            "money": "none: catalogue metadata over HTTP",
            "seconds": round(time.time() - started, 1),
        }
    )
    record["result_manifest"] = result_manifest(mf.traced_reads(), harvest)
    path = save_result(NAME, record)
    print(
        f"wrote {path}: {search['terms_with_both']} ENCODE terms with both modalities, "
        f"{len(absent)} absent by the rule, {len(siblings)} of those only because the clause "
        f"does not exclude siblings; answer {record['answer']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
