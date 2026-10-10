# SPDX-License-Identifier: AGPL-3.0-or-later
"""An inventory of public cell types whose ontology term is ABSENT from the AlphaGenome RNA
track metadata we hold — and of whether a paired accessibility + RNA dataset exists for any of
them. It adds no biosample, reads no RNA value, fits nothing, and its output is a list.

STOP-LIST CHECK. The 2026-09-27 stop list says "No more biosamples". This module does not add
one and does not run the out-of-domain test. It establishes whether that test is even possible,
so Albert's decision about lifting the clause is not taken blind.

WHY THIS EXISTS. M2h found the compiled rules fit measured RNA better when scaled by one fitted
constant (+0.2146; held-out Delta MSE +0.04339 [98.33% CI +0.0379, +0.04894]), and M2i
replicated it on RNA fresh to this project (+0.03727 [+0.03202, +0.04299], nothing fitted).
lane-trainexp then measured that independence from AlphaGenome's training data is UNKNOWN for
13 of 13: `data/cache/entex/alphagenome_track_metadata_copy.csv` carries no ENCODE experiment
accession on any `rna_seq` row, so only biosample+assay matching is possible, and a probable
track exists for all thirteen. The deeper point is not fixable by accession: AlphaGenome was
trained on RNA-seq of those same 13 cell types, so an accession-level "not this track" would
still not make the comparison independent. Only a cell type OUTSIDE the training domain would.

WHAT "ABSENT" CAN AND CANNOT MEAN. **The metadata file is what we hold, not what AlphaGenome
saw.** It is a copy of the published track list. A cell type absent from it is absent from the
track list; that is NOT the same as absent from the training set, which is not published at the
granularity that would settle it. Every reading this module produces is a statement about a
published catalogue, and it licenses no claim about the model. Said again in the result files.

THE MATCH FIELD. The file's columns are `name`, `strand`, `Assay title`, `ontology_curie`,
`biosample_name`, `biosample_type`, `biosample_life_stage`, `gtex_tissue`, `data_source`,
`endedness`, `genetically_modified`, `nonzero_mean`, `output`, `histone_mark`. We match on
`ontology_curie`, not `biosample_name`: the CURIE is a controlled identifier that ENCODE, HCA,
GEO's richer submissions and Tabula Sapiens also carry, so it can be compared across sources,
whereas `biosample_name` is free text ("T-cell", "endothelial cell of umbilical vein") and a
source that writes it differently would read as absent for a spelling reason. Every `rna_seq`
row in the file has a non-empty `ontology_curie`, so nothing is lost by using it.

THE ONTOLOGY RELATION, which is the rule's one real subtlety. CURIE equality alone would call
a cell type absent when the metadata covers its parent. A child of a covered parent is NOT
absent: AlphaGenome's `UBERON:0002048` (lung) track was trained on lung RNA, so "type II
pneumocyte" is inside that domain even though its CL term is not in the list. The converse also
disqualifies: a candidate whose term is a PARENT of a covered term has the covered cells inside
its own territory. So a candidate is absent only when no covered CURIE is equal to it, an
ancestor of it, or a descendant of it. We hold no ontology graph in this repository, so the
relation cannot be computed here: each candidate must DECLARE the covered CURIEs it is related
to, with the evidence, and `absent_by_rule` is False as soon as one is declared. An
undeclared-and-unchecked candidate is recorded as `relation_unchecked`, never as absent.
"""

from __future__ import annotations

import csv
import hashlib
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from pathlib import Path

METADATA_PATH = Path("data/cache/entex/alphagenome_track_metadata_copy.csv")

#: The one `output` value whose rows define the RNA training domain we compare against.
RNA_OUTPUT = "rna_seq"

#: The field matched on, and the reason, restated for the result file.
MATCH_FIELD = "ontology_curie"

SELECTION_RULE = {
    "version": "1",
    "question": (
        "Does any public cell type exist that is ABSENT from AlphaGenome's RNA track metadata "
        "as we hold it, AND has public accessibility (DNase or ATAC) and RNA-seq from the same "
        "source?"
    ),
    "match_field": MATCH_FIELD,
    "match_field_reason": (
        "ontology_curie is a controlled identifier shared with ENCODE, HCA and Tabula Sapiens; "
        "biosample_name is free text and a different spelling would read as absent."
    ),
    "covered_set": (
        "the distinct ontology_curie values of every row whose output == 'rna_seq' in "
        "data/cache/entex/alphagenome_track_metadata_copy.csv"
    ),
    "absent": (
        "A candidate CURIE is ABSENT only when no covered CURIE is (a) equal to it, (b) an "
        "ancestor of it, or (c) a descendant of it. A child term of a covered parent is NOT "
        "absent: the parent's RNA track was trained on tissue that contains the child. A parent "
        "term of a covered child is NOT absent either: the covered cells lie inside it."
    ),
    "relation_evidence": (
        "No ontology graph is held in this repository, so the ancestor/descendant relation is "
        "not computed. Each candidate declares the covered CURIEs it is related to, with the "
        "evidence. One declared relation makes absent_by_rule False. A candidate whose relation "
        "was not checked is recorded relation_unchecked and is NOT counted absent."
    ),
    "accessibility": (
        "A DNase-seq or ATAC-seq (including scATAC-seq) experiment on the same biosample term, "
        "with processed per-region signal or peaks obtainable without an application."
    ),
    "rna": (
        "An RNA-seq experiment (bulk or single-cell with a usable pseudobulk) on the same biosample term."
    ),
    "same_source": (
        "Both experiments published by one archive or consortium under one accession family and "
        "one processing pipeline, so the pair needs no cross-archive harmonisation: e.g. two "
        "ENCODE experiments on one ENCODE biosample, or one HCA/Tabula Sapiens project carrying "
        "both modalities."
    ),
    "licence": (
        "Open access, no application, no data-use agreement, no controlled-access tier. dbGaP, "
        "EGA and any 'apply for access' resource do not qualify."
    ),
    "minimum": "One cell type is enough to answer YES. Zero is NO_SOURCE_QUALIFIES.",
    "acceptable_answer": (
        "NO_SOURCE_QUALIFIES is a complete answer: Albert's biosample question then closes without a lift."
    ),
    "wording_limit": (
        "The metadata file is what we hold, not what AlphaGenome saw. 'Absent from the metadata' "
        "is not 'absent from the training set'. This inventory licenses no claim about the model."
    ),
    "no_biosample": ("This lane adds no biosample and runs no out-of-domain test. 2026-09-27 stop list."),
}


def rna_rows(path: Path | str = METADATA_PATH) -> list[dict[str, str]]:
    """Every row of the metadata whose ``output`` is :data:`RNA_OUTPUT`."""
    with Path(path).open(newline="") as handle:
        return [row for row in csv.DictReader(handle) if row.get("output") == RNA_OUTPUT]


def covered_curies(path: Path | str = METADATA_PATH) -> set[str]:
    """The set the rule calls covered: distinct ``ontology_curie`` over the RNA rows."""
    return {row[MATCH_FIELD].strip() for row in rna_rows(path) if row[MATCH_FIELD].strip()}


def curie_digest(curies: Iterable[str]) -> str:
    """A sha256 over the sorted covered set, so a result names the set it was judged against."""
    joined = "\n".join(sorted(curies))
    return hashlib.sha256(joined.encode()).hexdigest()


@dataclass(frozen=True)
class Candidate:
    """One public cell type considered, and what the rule says about it.

    ``related_covered`` holds covered CURIEs declared equal to, ancestor of, or descendant of
    ``curie``, each with its evidence. Empty means the relation was NOT checked, which the rule
    records as ``relation_unchecked`` rather than as absence.
    """

    curie: str
    label: str
    source: str
    accessibility: str | None = None
    rna: str | None = None
    licence: str | None = None
    size: str | None = None
    related_covered: tuple[str, ...] = ()
    relation_checked: bool = False
    notes: str = ""
    disqualified: str | None = None
    extra: dict[str, str] = field(default_factory=dict)


def verdict(candidate: Candidate, covered: Sequence[str] | set[str]) -> dict[str, object]:
    """Apply the committed rule to one candidate. Returns the row that goes in the list."""
    covered = set(covered)
    declared = tuple(c for c in candidate.related_covered if c in covered)
    in_metadata = candidate.curie in covered
    if not candidate.relation_checked:
        status = "relation_unchecked"
        absent = False
    elif in_metadata:
        status = "curie_present"
        absent = False
    elif declared:
        status = "related_to_covered"
        absent = False
    else:
        status = "absent"
        absent = True
    qualifies = bool(
        absent
        and candidate.accessibility
        and candidate.rna
        and candidate.licence
        and not candidate.disqualified
    )
    return {
        "curie": candidate.curie,
        "label": candidate.label,
        "source": candidate.source,
        "accessibility": candidate.accessibility,
        "rna": candidate.rna,
        "licence": candidate.licence,
        "size": candidate.size,
        "curie_in_metadata": in_metadata,
        "related_covered": list(declared),
        "relation_checked": candidate.relation_checked,
        "absent_status": status,
        "absent_by_rule": absent,
        "disqualified": candidate.disqualified,
        "qualifies": qualifies,
        "notes": candidate.notes,
        **candidate.extra,
    }


def inventory(candidates: Iterable[Candidate], covered: Sequence[str] | set[str]) -> dict[str, object]:
    """The whole list, plus the one answer the lane was asked for."""
    covered = set(covered)
    rows = [verdict(c, covered) for c in candidates]
    winners = [r for r in rows if r["qualifies"]]
    return {
        "covered_curie_count": len(covered),
        "covered_curie_digest": curie_digest(covered),
        "match_field": MATCH_FIELD,
        "candidates": rows,
        "qualifying_count": len(winners),
        "answer": "QUALIFIES" if winners else "NO_SOURCE_QUALIFIES",
        "wording_limit": SELECTION_RULE["wording_limit"],
    }
