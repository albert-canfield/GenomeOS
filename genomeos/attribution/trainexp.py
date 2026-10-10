# SPDX-License-Identifier: AGPL-3.0-or-later
"""POST HOC, DECIDING NOTHING: whether M2i's fresh RNA is independent of AlphaGenome's training data.

This module QUALIFIES M2i (`f6e71c6`). It does not revise it, and nothing here can change M2i's
registered reading, its outcome or its numbers. M2i's own files are read-only to this lane.

THE GAP IT ADDRESSES. M2i showed that M2h's one fitted constant (+0.2145911381525326, `ac5580b`)
holds on thirteen ENCODE RNA-seq experiments THIS PROJECT had never read, with non-overlap proved at
experiment, file and replicate-library level against a 954-accession read set. But the rules whose
magnitude is being tested originate from AlphaGenome, and AlphaGenome was trained on ENCODE RNA-seq
of these very cell types. A measurement fresh to this project may still sit inside the MODEL's
training set, in which case the replication is not independent of the predictor's origin. The words
"trained" and "training exposure" appear nowhere in M2i's registration except for its gene split.

WHAT CAN AND CANNOT BE MEASURED FROM THE RECORD. AlphaGenome's recorded output-track metadata is on
disk (`METADATA`). It carries, per track, `name`, `strand`, `Assay title`, `ontology_curie`,
`biosample_name`, `biosample_type`, `biosample_life_stage`, `gtex_tissue`, `data_source`,
`endedness`, `genetically_modified`, `nonzero_mean`, `output` and `histone_mark`. It carries NO
ENCODE experiment, file, library or biosample accession for any RNA-seq track: of its 667 `rna_seq`
rows (371 distinct track names, each on two strands) exactly 0 carry one, and the only accessions
anywhere in the file are the six PRO-cap track names (`ACCESSION_BEARING_OUTPUTS`). So accession
matching -- the only method that could settle whether a specific experiment is behind a specific
track -- IS NOT AVAILABLE for any RNA-seq track, and this lane reports that rather than inferring a
training set from a publication's prose.

THREE CLASSES, and the one that must not be collapsed. `CLASS_MATCHED` is a same-accession match and
is strong. `CLASS_PROBABLE` is a same-biosample-and-assay match with no accession available and is
SUGGESTIVE ONLY: it establishes that AlphaGenome has an output track of that biosample and that
assay, hence was trained on SOME ENCODE RNA-seq of it, but NOT that the particular experiment M2i
read is the one behind that track. `CLASS_UNKNOWN` is no basis to decide. **`CLASS_UNKNOWN` IS NOT
"NOT A TRAINING TRACK"**, and neither is `CLASS_PROBABLE` a proof of exposure; conflating either
with a verdict is the single error this module exists to prevent.
"""

from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Any

from genomeos.attribution import context_evidence as ce

#: This lane's name, as the board and the result carry it.
LANE = "lane-trainexp"

#: AlphaGenome's own recorded output-track metadata, already on disk and GIT-IGNORED (.gitignore:20).
#: The same authority `context_evidence.TRACK_METADATA` and `cellcover.TRACK_METADATA` read.
METADATA = Path("data/cache/entex/alphagenome_track_metadata_copy.csv")

#: M2i's thirteen fresh experiments, and the read set they were proved disjoint from.
FRESH_SOURCES = Path("data/results/freshrna_sources.json")

#: Items 7/7d/7h's twenty-seven already-read quantifications (thirteen experiments).
READ_SOURCES = Path("data/results/readergate_rna_sources.json")

#: The `output` value of every track whose `name` carries an ENCODE accession, measured not assumed.
ACCESSION_BEARING_OUTPUTS = ("procap",)

#: The `output` value of the tracks this lane matches against.
RNA_OUTPUT = "rna_seq"

#: An ENCODE accession, in any of the four forms a source file names.
ACCESSION = re.compile(r"^ENC(?:SR|FF|BS|LB)[0-9A-Z]{6}$")

CLASS_MATCHED = "matched"
CLASS_PROBABLE = "probable"
CLASS_UNKNOWN = "unknown"

#: What each class means and how strong it is. Carried into the result verbatim.
CLASSES: dict[str, str] = {
    CLASS_MATCHED: (
        "MATCHED: the recorded track metadata names this experiment's own ENCODE accession, so the"
        " experiment IS one AlphaGenome records as an output track of its own. Accession matching is"
        " STRONG: it identifies the experiment, not merely its kind"
    ),
    CLASS_PROBABLE: (
        "PROBABLE: the recorded track metadata carries NO accession for any RNA-seq track, so the"
        " strong method is unavailable; a track of the SAME biosample (ontology CURIE) and the SAME"
        " assay title exists. Biosample+assay matching is SUGGESTIVE ONLY: it establishes that"
        " AlphaGenome has an output track of this biosample and assay, hence was trained on SOME"
        " ENCODE RNA-seq of it, but NOT that this experiment is the one behind that track"
    ),
    CLASS_UNKNOWN: (
        "UNKNOWN: no basis in the record to decide. The metadata has no accession to match and no"
        " RNA-seq track of this biosample under this assay. THIS IS NOT 'NOT A TRAINING TRACK': it"
        " is the absence of a basis, and it must never be read, reported or summarised as evidence"
        " of independence"
    ),
}

#: The method a match was made by, named per experiment so a reader never has to infer it.
METHOD_ACCESSION = "ENCODE accession, as the track metadata carries it (STRONG)"
METHOD_BIOSAMPLE_ASSAY = (
    "biosample ontology CURIE + assay title, the metadata carrying no accession (SUGGESTIVE ONLY)"
)
METHOD_NONE = "none available: no accession in the metadata and no track of this biosample and assay"

#: What is knowable about AlphaGenome's training-set composition from this project's record alone.
TRAINING_SET_RECORD = (
    "AlphaGenome's training-set composition is NOT documented at track level in anything this"
    " project holds. What IS knowable from the record is which OUTPUT tracks the model carries"
    " metadata for, by assay title and biosample ontology CURIE, and that the only ENCODE accessions"
    " the metadata names at all are the six PRO-cap track names"
    " (ENCSR072YCM, ENCSR182QNJ, ENCSR740IPL, ENCSR797DEF, ENCSR801ECP, ENCSR860TYZ) -- none of them"
    " RNA-seq. An output track is a training TARGET, so a biosample+assay match says the model saw"
    " some ENCODE RNA-seq of that biosample and assay; it does not say which experiment. This lane"
    " therefore reports 'unknown' for the experiment-level question and does NOT infer a training"
    " set from a publication's prose."
)

#: The two cells M2i itself disclosed as sharing a biosample with an already-read experiment.
SHARED_BIOSAMPLE_CELLS = ("cardiac muscle cell", "hepatocyte")

#: The label every quantity this module produces carries.
POST_HOC = (
    "POST HOC, computed after M2i's registered result (f6e71c6) and its post-hoc stage (0119387)"
    " were committed. It QUALIFIES M2i and DECIDES NOTHING: it does not revise M2i's reading, its"
    " outcome, its interval or its scale, and nothing is fitted anywhere in it."
)


class MetadataAbsentError(RuntimeError):
    """The recorded track metadata is not on this machine, so nothing is matched and nothing claimed."""


def cell_terms() -> dict[str, tuple[str, str]]:
    """freshrna/readergate cell id -> (the ENCODE biosample term name, its ontology CURIE).

    Derived mechanically from `context_evidence.READER_TERMS` through the same `ident` the compiler
    uses, so no hand-written table can drift from the authority.
    """
    return {ce.ident(name): (name, term) for name, term in ce.READER_TERMS.items()}


def read_tracks(path: Path = METADATA) -> list[dict[str, str]]:
    """Every row of AlphaGenome's recorded output-track metadata, as it is on disk."""
    if not path.exists():
        raise MetadataAbsentError(
            f"{path} is not on this machine (it is git-ignored machine-local data). Nothing is"
            " matched and no independence is claimed; this is reported, not fetched."
        )
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh))


def rna_tracks(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """The `rna_seq` rows alone: the tracks a measured RNA quantification could correspond to."""
    return [r for r in rows if r.get("output") == RNA_OUTPUT]


def accession_index(rows: list[dict[str, str]]) -> dict[str, list[str]]:
    """Every ENCODE accession the metadata names anywhere -> the `output` values it appears under.

    The whole point of measuring this is that for `rna_seq` it comes back EMPTY, which is why the
    strong matching method is unavailable and why `CLASS_MATCHED` is 0.
    """
    out: dict[str, list[str]] = {}
    for r in rows:
        for field in (r.get("name") or "").split():
            if ACCESSION.match(field):
                outs = out.setdefault(field, [])
                if r.get("output") and r["output"] not in outs:
                    outs.append(r["output"])
    return out


def biosample_assay_index(rna: list[dict[str, str]]) -> dict[tuple[str, str], list[str]]:
    """(ontology CURIE, assay title) -> the distinct RNA-seq track names recorded for it."""
    out: dict[tuple[str, str], list[str]] = {}
    for r in rna:
        names = out.setdefault((r.get("ontology_curie", ""), r.get("Assay title", "")), [])
        if r.get("name") and r["name"] not in names:
            names.append(r["name"])
    return out


def curie_index(rna: list[dict[str, str]]) -> dict[str, list[str]]:
    """ontology CURIE -> the assay titles the metadata records RNA-seq tracks of it under.

    Reported beside an `unknown` so a reader can see WHY there was no basis: the biosample may be
    absent from the metadata altogether, or present under another assay.
    """
    out: dict[str, list[str]] = {}
    for r in rna:
        titles = out.setdefault(r.get("ontology_curie", ""), [])
        if r.get("Assay title") and r["Assay title"] not in titles:
            titles.append(r["Assay title"])
    return out


def classify_one(
    cell: str,
    experiment: str,
    biosample: str,
    curie: str,
    assay: str,
    accessions: dict[str, list[str]],
    by_biosample_assay: dict[tuple[str, str], list[str]],
    by_curie: dict[str, list[str]],
    files: int = 0,
) -> dict[str, Any]:
    """One experiment's class, the method used, and the evidence it rests on. Decides nothing."""
    if experiment in accessions:
        return {
            "cell": cell,
            "experiment": experiment,
            "biosample": biosample,
            "ontology_curie": curie,
            "assay": assay,
            "files": files,
            "class": CLASS_MATCHED,
            "method": METHOD_ACCESSION,
            "evidence": f"the metadata names {experiment} under output(s) "
            + ", ".join(accessions[experiment]),
        }
    tracks = by_biosample_assay.get((curie, assay), [])
    if tracks:
        return {
            "cell": cell,
            "experiment": experiment,
            "biosample": biosample,
            "ontology_curie": curie,
            "assay": assay,
            "files": files,
            "class": CLASS_PROBABLE,
            "method": METHOD_BIOSAMPLE_ASSAY,
            "evidence": f"{len(tracks)} recorded RNA-seq track name(s) of {curie} under '{assay}': "
            + ", ".join(sorted(tracks)),
        }
    return {
        "cell": cell,
        "experiment": experiment,
        "biosample": biosample,
        "ontology_curie": curie,
        "assay": assay,
        "files": files,
        "class": CLASS_UNKNOWN,
        "method": METHOD_NONE,
        "evidence": (
            f"no RNA-seq track of {curie} under '{assay}'; the assay titles the metadata does record"
            f" RNA-seq tracks of {curie} under are: "
            + (", ".join(sorted(by_curie.get(curie, []))) or "none at all")
        ),
    }


def tally(entries: list[dict[str, Any]]) -> dict[str, dict[str, int]]:
    """Per class, the experiments and the quantification files in it."""
    out = {k: {"experiments": 0, "files": 0} for k in CLASSES}
    for e in entries:
        out[e["class"]]["experiments"] += 1
        out[e["class"]]["files"] += int(e.get("files") or 0)
    return out


def independence_sentence(entries: list[dict[str, Any]]) -> str:
    """The sentence to be carried BESIDE every M2i 'better', filled in from what was measured.

    The template M2i's qualification takes is "fresh to this project; its independence from
    AlphaGenome's training data is <measured: k of 13 / unknown>". `k of 13` would be the count
    settled by accession matching. Because the recorded metadata carries no accession for any
    RNA-seq track, nothing is settled either way and the honest filling is `unknown for k of 13`,
    with the biosample+assay finding stated as the suggestive thing it is.
    """
    n = len(entries)
    matched = sum(1 for e in entries if e["class"] == CLASS_MATCHED)
    probable = sum(1 for e in entries if e["class"] == CLASS_PROBABLE)
    if matched == n and n:
        return (
            f"fresh to this project; its independence from AlphaGenome's training data is {matched}"
            f" of {n} inside the model's recorded tracks by ENCODE accession"
        )
    settled = matched
    return (
        "fresh to this project; its independence from AlphaGenome's training data is unknown for"
        f" {n - settled} of {n}"
        " -- AlphaGenome's recorded track metadata carries no ENCODE accession for any RNA-seq"
        f" track, and a track of the same biosample and assay exists for {probable} of {n}, which is"
        " suggestive of training exposure and settles nothing either way"
    )


def subset_indices(cells: tuple[str, ...] | list[str], exclude: tuple[str, ...] | list[str]) -> list[int]:
    """The positions of the cells to KEEP, in the order M2i's per-cell vectors carry them."""
    drop = {ce.ident(c) for c in exclude} | set(exclude)
    return [i for i, c in enumerate(cells) if c not in drop]


def excluded_cells(cells: tuple[str, ...] | list[str], exclude: tuple[str, ...] | list[str]) -> list[str]:
    """The cells actually dropped, so a result never reports a subset it did not take."""
    keep = set(subset_indices(cells, exclude))
    return [c for i, c in enumerate(cells) if i not in keep]


#: What the eleven-cell subset is and is not, carried into the result verbatim.
ELEVEN_CELL_NOTE = (
    "POST HOC: M2i's comparison recomputed on eleven of the thirteen cells, dropping cardiac muscle"
    " cell and hepatocyte -- the two whose biosamples (ENCBS077RUJ, ENCBS343AKO, ENCBS520VFV,"
    " ENCBS806KLM) an already-read experiment also used, which is the one non-independence M2i"
    " itself disclosed (`non_overlap.cells_sharing_a_read_biosample`). It isolates that disclosed"
    " overlap and NOTHING ELSE: it is not a test of independence from AlphaGenome's training data,"
    " which the record cannot settle. Everything about it is post hoc. NOTHING IS FITTED: the scale"
    " is M2i's constant, D' is recentred within each gene across the eleven kept cells exactly as"
    " M2i's own nine-polyA-plus subset recentres, and the interval is M2i's gene-level percentile"
    " bootstrap at its registered resamples, seed and percentiles."
)


# ---- what was MEASURED, 2026-10-10, and why it is recorded here ---------------------------------
#
# These are the values `scripts/trainexp_match.py` computed on 2026-10-10 from the metadata and the
# quantifications named above. They are recorded in committed code because the result file they were
# written into, `data/results/trainexp_independence.json`, COULD NOT BE COMMITTED: `commit_own.sh`
# refused it (exit 2) because another lane's uncommitted `genomeos/attribution/onetarget2.py` -- a
# fix held outside a commit, which this lane is not permitted to touch -- sits on the static import
# closure of anything that reaches `context_evidence` (the authority for the thirteen biosample
# CURIEs), so the revision the result would stamp names a tree that does not contain all the code on
# its counting path. The refusal is correct and was not routed around: `--force` was not used, the
# result file stays written locally and uncommitted, and these constants plus the commit message are
# the record. Re-running the script once `onetarget2.py` is committed reproduces every value here.

#: The recorded track metadata as it was read, measured not assumed.
METADATA_AS_READ = {
    "rows": 2841,
    "rna_seq_rows": 667,
    "rna_seq_distinct_track_names": 371,
    "rna_seq_rows_carrying_an_encode_accession": 0,
    "accessions_anywhere_in_the_file": (
        "ENCSR072YCM",
        "ENCSR182QNJ",
        "ENCSR740IPL",
        "ENCSR797DEF",
        "ENCSR801ECP",
        "ENCSR860TYZ",
    ),
    "those_accessions_output": "procap, none of them rna_seq",
    "accession_matching_available_for_rna_seq": False,
}

#: M2i's thirteen fresh experiments: experiments and files per class, by biosample CURIE + assay.
FRESH_TALLY = {
    CLASS_MATCHED: {"experiments": 0, "files": 0},
    CLASS_PROBABLE: {"experiments": 13, "files": 22},
    CLASS_UNKNOWN: {"experiments": 0, "files": 0},
}

#: Items 7/7d/7h's twenty-seven already-read quantifications, thirteen experiments, all polyA plus.
#: The one `unknown` is astrocyte (ENCSR920UAO, CL:0000127), which the metadata records RNA-seq
#: tracks of under total RNA-seq but not under polyA plus RNA-seq.
ALREADY_READ_TALLY = {
    CLASS_MATCHED: {"experiments": 0, "files": 0},
    CLASS_PROBABLE: {"experiments": 12, "files": 26},
    CLASS_UNKNOWN: {"experiments": 1, "files": 1},
}

#: THE SENTENCE to be carried beside every M2i "better", word for word. Produced by
#: `independence_sentence` on the thirteen fresh entries as they classified, and asserted equal to it.
INDEPENDENCE_SENTENCE = (
    "fresh to this project; its independence from AlphaGenome's training data is unknown for 13 of 13"
    " -- AlphaGenome's recorded track metadata carries no ENCODE accession for any RNA-seq track, and"
    " a track of the same biosample and assay exists for 13 of 13, which is suggestive of training"
    " exposure and settles nothing either way"
)

#: The eleven-cell subset, POST HOC and deciding nothing: M2i's constant scale, M2i's interval method.
ELEVEN_CELL_MEASURED = {
    "delta_mse": 0.0459557845802319,
    "ci": (0.040272388859847665, 0.05209648171116512),
    "genes": 8328,
    "pairs": 91608,
    "thirteen_cell_delta_mse": 0.03727241565491676,
    "thirteen_cell_ci": (0.03202429380417028, 0.04299157342045463),
    "thirteen_reproduces_the_committed_result": True,
    "cells_dropped": ("hepatocyte", "cardiac_muscle_cell"),
}
