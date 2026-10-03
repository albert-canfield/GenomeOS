# SPDX-License-Identifier: AGPL-3.0-or-later
"""A documentation census of the knockdown fields of the Replogle 2022 pseudobulk file.

    uv run --frozen python scripts/n1_kd_semantics_census.py

The question this answers, asked before any study is authorised: can the semantics of the CRISPRi
knockdown-efficacy fields be established FROM THE PRODUCERS' OWN DOCUMENTATION, so that an
eligibility threshold could be fixed in advance on a known scale? N1's committed registration
refused those fields as undocumented (data/results/n1_registration_amendment_2.json,
`run_v3.never_read`: "any knockdown field"), and a threshold on an unknown scale constrains nothing.

THIS CENSUS CONTAINS NO MEASUREMENT OF ANY STUDY VALUE AND OPENS NO STUDY FILE. Not one byte of
`K562_gwps_normalized_bulk_01.h5ad`, `K562_gwps_raw_bulk_01.h5ad` or any published p-value table is
read -- not a header, not a row, not a distribution. What is read is documentation: two Figshare+
deposit records (the public metadata API, not the data files), the paper's full text from Europe PMC,
and the two code repositories the paper names as the authors' own. Every NOT ESTABLISHED verdict is a
mechanical count over those documents (`absent`), not an assertion, and every ESTABLISHED verdict
carries the producers' sentence verbatim.

What it deliberately does NOT do: choose a threshold. That is a registration decision.

Writes: data/results/n1_knockdown_field_semantics_census.json
Caches: data/cache/kdsemantics/ (the fetched documents, so the counts are reproducible)
"""

from __future__ import annotations

import json
import re
import sys
import tarfile
import urllib.request
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

NAME = "n1_knockdown_field_semantics_census"
CACHE = Path("data/cache/kdsemantics")
ORIGINAL = RESULTS_DIR / "n1_registration.json"
AMENDMENT_1 = RESULTS_DIR / "n1_registration_amendment_1.json"
AMENDMENT_2 = RESULTS_DIR / "n1_registration_amendment_2.json"

#: This lane's own files. Everything else uncommitted in this shared checkout belongs to a peer.
OWN_CODE = (
    "scripts/n1_kd_semantics_census.py",
    "tests/test_kd_semantics_census.py",
)

#: The documents, each one the producers' own: their deposit records, their paper, their code.
DOCS: dict[str, dict[str, str]] = {
    "figshare_20029387": {
        "url": "https://api.figshare.com/v2/articles/20029387",
        "file": "figshare_20029387.json",
        "what": "the deposit that hosts the pseudobulk files, description text (public metadata API)",
    },
    "figshare_21632564": {
        "url": "https://api.figshare.com/v2/articles/21632564",
        "file": "figshare_21632564.json",
        "what": "the producers' supplemental deposit, description text (public metadata API)",
    },
    "pmc9380471": {
        "url": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9380471/fullTextXML",
        "file": "pmc9380471.xml",
        "what": "Replogle et al. 2022, Cell 185:2559, full text with STAR Methods and figure legends",
    },
    "code_perturbseq_gi": {
        "url": "https://codeload.github.com/thomasmaxwellnorman/Perturbseq_GI/tar.gz/refs/heads/master",
        "file": "Perturbseq_GI-master.tar.gz",
        "what": "the 'Perturb-seq analysis codebase' the paper's key resources table names (Norman 2019)",
    },
    "code_guide_calling": {
        "url": "https://codeload.github.com/josephreplogle/guide_calling/tar.gz/refs/heads/master",
        "file": "guide_calling-master.tar.gz",
        "what": "the 'sgRNA assignment scripts' the paper's key resources table names",
    },
}

#: The field names as they are spelled in the file N1 refused to read them from.
FIELDS_2022 = ("control_expr", "fold_expr", "pct_expr", "num_cells_filtered")

#: Producer sentences the census requires to be present, verbatim. A missing one is a hard failure:
#: the census would otherwise record an ESTABLISHED verdict against a quote it never found.
QUOTES: dict[str, str] = {
    "knockdown_definition": (
        "Knockdown was computed as the ratio of mean (unnormalized) expression of the target gene "
        "within perturbed cells vs. that in cells with non-targeting sgRNAs."
    ),
    "fractional_change_definition": (
        "The fractional change in expression is defined as the expression in the targeted cells minus "
        "the expression in non-targeting cells, relative to the expression in the non-targeting cell "
        "population"
    ),
    "producer_threshold_30": (
        "an on-target knockdown, if measured, of at least 30% (i.e. the target of perturbation was "
        "either knocked down by at least 30% or was not detected, a broad attempt to remove "
        "non-functional perturbations)"
    ),
    "producer_threshold_60": "an on-target knockdown of at least 60%",
    "median_knockdown": (
        "We observed a median target knockdown of 85.5% in K562 cells and 91.6% in RPE1 cells"
    ),
    "obs_is_all_the_deposit_says": ("the .obs annotation details single-cells/pseudobulk populations"),
    "code_availability": (
        "Our codebase for Perturb-seq analysis is available at "
        "https://github.com/thomasmaxwellnorman/Perturbseq_GI and "
        "https://github.com/josephreplogle/guide_calling"
    ),
}

#: The 2019 counterpart the committed record says does not exist (amendment 1, item 2: "pct_expr has
#: no 2019 counterpart"). It does, in the codebase the paper names, and the census records the line.
PCT_2019_LINE = (
    "mean_pop.cells['pct_first_expr'] = "
    "mean_pop.cells['diff_first_expr']/mean_pop.cells['control_first_expr']"
)


def fetch(key: str) -> Path:
    """One document, cached. Documentation only: no study data file is ever a target here."""
    spec = DOCS[key]
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / spec["file"]
    if not path.exists():
        req = urllib.request.Request(spec["url"], headers={"User-Agent": "GenomeOS/kdsemantics"})
        with urllib.request.urlopen(req, timeout=180) as r:  # noqa: S310 - fixed https documents
            path.write_bytes(r.read())
    return path


def plain_text(raw: bytes) -> str:
    """Tags out, entities for the two that matter, runs of whitespace to one space."""
    t = raw.decode("utf-8", errors="replace")
    t = re.sub(r"<[^>]+>", " ", t)
    t = t.replace("&gt;", ">").replace("&lt;", "<").replace("&amp;", "&")
    t = t.replace("—", "-").replace("–", "-").replace("’", "'")
    t = t.replace("“", '"').replace("”", '"').replace("‘", "'")
    return re.sub(r"\s+", " ", t)


def archive_text(path: Path) -> str:
    """Every text member of a source tarball, concatenated. Notebooks are JSON, so this sees them."""
    out: list[str] = []
    with tarfile.open(path, "r:gz") as tf:
        for m in tf.getmembers():
            if not m.isfile() or m.size > 8_000_000:
                continue
            f = tf.extractfile(m)
            if f is None:
                continue
            chunk = f.read()
            if b"\x00" in chunk[:2048]:
                continue  # a binary member (the PDFs of guide calls)
            out.append(f"\n### {m.name}\n")
            out.append(chunk.decode("utf-8", errors="replace"))
    return "".join(out)


def occurrences(text: str, needles: tuple[str, ...]) -> dict[str, int]:
    """How many times each spelling appears. Zero is the evidence a NOT ESTABLISHED verdict rests on."""
    return {n: len(re.findall(re.escape(n), text)) for n in needles}


def find_quote(text: str, quote: str) -> bool:
    """The producer's sentence, present or not, insensitive to the whitespace a converter introduces."""
    pattern = r"\s*".join(re.escape(w) for w in quote.split())
    return re.search(pattern, text, re.I) is not None


def census(documents: dict[str, str]) -> dict[str, Any]:
    """The verdict per field, from the counts and the quotes alone.

    The distinction the whole census exists to keep: the producers document a QUANTITY ("on-target
    knockdown"), its units, its direction and two thresholds they themselves used. They document no
    COLUMN of the distributed file. A quantity whose carrier is unnamed cannot be gated in advance,
    because the number would be fixed against a scale no document ties to the bytes being read.
    """
    paper = documents["pmc9380471"]
    deposits = documents["figshare_20029387"] + documents["figshare_21632564"]
    code = documents["code_perturbseq_gi"] + documents["code_guide_calling"]
    everywhere = paper + deposits + code

    counts = {
        "paper": occurrences(paper, FIELDS_2022),
        "deposit_descriptions": occurrences(deposits, FIELDS_2022),
        "named_producer_code": occurrences(code, FIELDS_2022),
    }
    quotes = {k: find_quote(everywhere, v) for k, v in QUOTES.items()}
    absent_everywhere = [f for f in FIELDS_2022 if occurrences(everywhere, (f,))[f] == 0]

    quantity = {
        "item": "on-target knockdown, the quantity",
        "verdict": "ESTABLISHED",
        "definition": QUOTES["knockdown_definition"],
        "definition_where": (
            "Replogle et al. 2022, Cell 185:2559 (PMC9380471), STAR Methods, 'Leverage scores for "
            "quantifying perturbation penetrance and variability'"
        ),
        "units": (
            "dimensionless as a ratio (perturbed mean over non-targeting mean of the target gene), "
            "reported in the paper as a percent reduction; 100% knockdown is a ratio of 0 and a "
            "fractional change of -1"
        ),
        "units_where": "the same section, with Figure S3's legend for the fractional change",
        "direction": (
            "as a RATIO, smaller is better knockdown; as a PERCENT KNOCKDOWN or a fractional change, "
            "larger magnitude of reduction is better knockdown. The two run opposite ways, which is "
            "exactly why a threshold cannot be fixed against an unnamed carrier"
        ),
        "direction_where": QUOTES["fractional_change_definition"] + " (-1 implies 100% knockdown).",
        "producers_own_thresholds": [
            {
                "level": "at least 30% on-target knockdown, or the target not detected",
                "used_for": "the paper's 'strong perturbations', its broad functional filter",
                "quote": QUOTES["producer_threshold_30"],
                "where": "PMC9380471, STAR Methods, 'Functional analyses of strong perturbations'",
            },
            {
                "level": "at least 60% on-target knockdown",
                "used_for": "the mitochondrial clustering of Figures 6A, 6C, S12A and S12E",
                "quote": QUOTES["producer_threshold_60"],
                "where": "PMC9380471, STAR Methods, the mitochondrial clustering criteria",
            },
        ],
        "scale_reference": {
            "quote": QUOTES["median_knockdown"],
            "where": "PMC9380471, Results, with Figure 1B",
        },
    }

    fields = [
        {
            "field": f,
            "where_it_lives": (
                "`.obs` of K562_gwps_normalized_bulk_01.h5ad (Figshare+ 20029387), one row per "
                "pseudobulk perturbation population -- NOT OPENED BY THIS CENSUS"
            ),
            "verdict": "NOT ESTABLISHED",
            "looked_at": [DOCS[k]["what"] for k in DOCS],
            "what_was_missing": (
                f"the spelling `{f}` occurs 0 times in the paper's full text (including STAR Methods, "
                "all figure legends and the key resources table), 0 times in either Figshare+ deposit "
                "description, and 0 times in both code repositories the paper names as the authors' "
                "own. No producer document defines it, gives its units, or states its direction"
            ),
            "counts": {k: counts[k][f] for k in counts},
        }
        for f in FIELDS_2022
    ]

    return {
        "quantity": quantity,
        "fields": fields,
        "fields_absent_from_every_producer_document": absent_everywhere,
        "counts": counts,
        "quotes_found": quotes,
        "deposit_says_only_this_about_obs": {
            "quote": QUOTES["obs_is_all_the_deposit_says"],
            "where": "Figshare+ 20029387 description (public metadata API), the whole of it on .obs",
            "verdict": "NOT ESTABLISHED",
            "reason": "it names no column",
        },
        "the_2022_pipeline_code_is_not_released": {
            "verdict": "ESTABLISHED (as an absence)",
            "quote": QUOTES["code_availability"],
            "where": "PMC9380471, 'Data and code availability'",
            "reading": (
                "the two repositories the authors release are Norman 2019's analysis codebase and the "
                "sgRNA caller. Neither contains any of the four spellings. The 2022 pseudobulk "
                "pipeline that wrote these columns is not among them, and the section's remedy is "
                "'available from the lead contact upon request'"
            ),
        },
        "correction_to_the_committed_record": {
            "what_the_record_says": (
                "data/results/n1_registration_amendment_1.json, changes[1].evidence[4]."
                "does_not_establish: 'pct_expr has no 2019 counterpart'"
            ),
            "what_this_census_finds": (
                "a counterpart does exist in the codebase the paper names, under the 2019 naming "
                "scheme, computed as (perturbed - control) / control"
            ),
            "line": PCT_2019_LINE,
            "where": (
                "thomasmaxwellnorman/Perturbseq_GI, GI_generate_populations.ipynb and "
                "GI_activation_of_neighboring_genes.ipynb"
            ),
            "verdict": "NOT ESTABLISHED for the 2022 file",
            "why_it_changes_nothing": (
                "the 2019 names carry a slot (`pct_first_expr`, `pct_second_expr`) because that "
                "experiment pairs two targets; the 2022 file's columns are spelled without one and "
                "appear in no document. Reading one family as the other is an inference from naming "
                "analogy across a different experiment and a changed pipeline -- the 2022 pipeline "
                "added per-gemgroup depth adjustment, so 'unnormalized' is ambiguous between raw and "
                "adjusted counts. Amendment 1's verdict therefore stands on a corrected premise, and "
                "amending the recorded sentence is a registration decision, not this census's"
            ),
        },
        "answer": {
            "question": (
                "can an eligibility gate with a knockdown threshold fixed in advance be built on this "
                "dataset's producer documentation?"
            ),
            "verdict": "NO, not on the file's knockdown columns",
            "because": (
                "the quantity is documented and the producers state their own thresholds, but no "
                "producer document names a single column of the distributed file. A threshold fixed "
                "against `fold_expr` would be a number on a scale no document ties to those bytes, "
                "and the ratio and the percent run in opposite directions, so even the sign of the "
                "comparison is an inference"
            ),
            "what_would_resolve_it": [
                "the 2022 pipeline code, or a statement by the authors defining the columns and units",
                "deriving knockdown from documented quantities instead -- the raw pseudobulk's target "
                "column against its non-targeting rows is the producers' own definition, computable "
                "without reading any undocumented column -- which is a separate registration, with "
                "its own outcome exposure, and is NOT authorised by this census",
            ],
            "not_chosen_here": "no threshold is chosen; that is a registration decision",
        },
    }


def main() -> int:
    paths = {k: fetch(k) for k in DOCS}
    documents = {
        k: (archive_text(p) if p.name.endswith(".tar.gz") else plain_text(p.read_bytes()))
        for k, p in paths.items()
    }
    body = census(documents)
    missing = [k for k, ok in body["quotes_found"].items() if not ok]
    if missing:
        print(f"REFUSED: producer sentences not found in the documents: {missing}", file=sys.stderr)
        return 1
    if body["fields_absent_from_every_producer_document"] != list(FIELDS_2022):
        print("a field name was found in a producer document; the census must be rewritten by hand")
        return 1

    payload: dict[str, Any] = {
        "result": NAME,
        "question": (
            "can the semantics of the CRISPRi knockdown-efficacy fields be established from the "
            "producers' own documentation, so that an eligibility threshold could be fixed in advance?"
        ),
        "contains_no_study_measurement": True,
        "no_study_file_was_opened": (
            "no byte of K562_gwps_normalized_bulk_01.h5ad, K562_gwps_raw_bulk_01.h5ad or any "
            "published p-value table was read: not a header, not a row, not a value, not a "
            "distribution. Only documentation was read"
        ),
        "documents_read": {k: {**DOCS[k], "bytes": paths[k].stat().st_size} for k in DOCS},
        "registration_this_answers_for": {
            "refusal_that_stands": (
                "data/results/n1_registration_amendment_2.json, run_v3.never_read: 'any knockdown "
                "field'; question_changes[0]: 'the undocumented knockdown fields are not read'"
            ),
            "nothing_here_amends_it": True,
        },
        **body,
    }
    payload["result_manifest"] = {
        "sources": [
            {"accession": "Replogle et al. 2022, Cell, full text (documentation)", "version": "PMC9380471"},
            {"accession": "Figshare+ 20029387 description (metadata API)", "version": "article record"},
            {"accession": "Figshare+ 21632564 description (metadata API)", "version": "article record"},
            {"accession": "thomasmaxwellnorman/Perturbseq_GI (producer code)", "version": "master"},
            {"accession": "josephreplogle/guide_calling (producer code)", "version": "master"},
            {
                "accession": "GenomeOS n1_registration, amendment 1 and amendment 2",
                "version": "the committed record",
            },
        ],
        "inputs": [
            mf.input_entry(ORIGINAL, partition=None),
            mf.input_entry(AMENDMENT_1, partition=None),
            mf.input_entry(AMENDMENT_2, partition=None),
            *[mf.input_entry(paths[k], partition=None, url=DOCS[k]["url"]) for k in DOCS],
        ],
        "assembly": "n/a: a documentation census; no coordinate is read",
        "coordinates": "n/a: a documentation census; no coordinate is read",
        "parameters": {
            "fields_censused": list(FIELDS_2022),
            "producer_documents": list(DOCS),
            "threshold_chosen": "none: choosing one is a registration decision",
        },
        "exclusions": [
            "every study data file: the two pseudobulk h5ads and the published p-value table are "
            "never opened, so no knockdown value of any kind is read",
        ],
        "partitions": "n/a: a census of documents, not of rows",
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }
    for p in (ORIGINAL, AMENDMENT_1, AMENDMENT_2):
        json.loads(p.read_text())  # the record this answers for must parse and be present
    out = save_result(NAME, payload)
    counts = body["counts"]
    print(f"wrote {out}")
    print(f"quantity: {body['quantity']['verdict']}; answer: {body['answer']['verdict']}")
    for f in FIELDS_2022:
        print(f"  {f:<20} NOT ESTABLISHED  " + json.dumps({k: counts[k][f] for k in counts}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
