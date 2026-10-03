# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""The first independent check of `silencer_like`: assay-validated ReSE silencers against the
model's signed direction.

`silencer_like` (genomeos/predict/enhancer_target.py) counts, per element, that on deleting it in
one deterministic run the single largest-moving gene in the 1 Mb window moved UP by at least
MIN_EFFECT log2. It is winner-take-all over one gene in one direction, ties go to activation, and
by review R4 no probability attaches to it. Nothing outside the model has ever been consulted about
a single one of the 230,839 elements that carry it.

THE TWO QUANTITIES ARE NOT THE SAME THING, and this module exists to compare them anyway, with
that difference registered first. The derived call is a *direction over a whole 1 Mb window* in a
sequence model. ReSE measures *one fragment's repressive ability in a lentiviral reporter assay*
that scores cell survival under transcriptional repression of an apoptosis-inducing fusion. A
fragment can repress a reporter it is cloned beside while the model's largest window mover in its
native locus is some other gene going the other way, and neither instrument is wrong when that
happens.

The module is split in two ON PURPOSE:

  * `scope()` states what will be measured and never sees a direction. It takes its element data
    only through `element_index()`, which returns `(start, end, id, has_predicted_target)` tuples --
    there is physically no action field in what `scope()` is handed, so a registration built from it
    cannot carry an outcome by accident.
  * `readings()` is the only function that reads a direction, and it is not called until the
    registration is committed.
"""

from __future__ import annotations

import bisect
import gzip
import json
import random
import re
from collections import Counter
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent.parent
RESE_CACHE = ROOT / "data/cache/silenceragree/rese_records.json"
ELEMENTS_DIR = ROOT / "data/knowledge/alphagenome/all_elements"
CCRE_DIR = ROOT / "data/results"

#: The ReSE interval is the record's own DESCRIPTION, and it is stated in GRCh37 on every record.
#: The GRCh38 coordinates come from `genomicinfo`, which RefSeq lifted; the description's own
#: numbers are GRCh37 and MUST NOT be joined against this project's GRCh38 elements.
RESE_DESCRIPTION = re.compile(r"^ReSE screen-validated silencer (GRCh\d+)_chr([0-9XYMT]+):(\d+)-(\d+)$")
#: Records that mention ReSE only in `otherdesignations` or in the summary prose are EXCLUDED: for
#: those the record is some larger functional element (a promoter, an enhancer) that merely contains
#: the assayed fragment, so `genomicinfo` is the parent's span and not the silencer's.
RESE_EXCLUDED_REASON = (
    "ReSE is named somewhere in the record but the record is not itself the assayed fragment, so its "
    "genomicinfo span is a parent region (measured: median 3.55x the fragment, up to 138.9x) and "
    "joining on it would attribute the assay to sequence the assay never carried"
)
#: ReSE validated in exactly these lines. The project scores four; the source covers two of them.
RESE_CELL_LINES = ("HepG2", "K562")
PROJECT_CELL_LINES = ("GM12878", "HepG2", "IMR-90", "K562")

CELL_CLAUSE = re.compile(r"identified as (?:a )?functional silencers? in ([^.]*)", re.IGNORECASE)

SOURCE = {
    "screen": "ReSE (repressive ability of silencer elements), a lentiviral screen scoring cell "
    "survival under transcriptional repression of an apoptosis-inducing fusion protein",
    "paper": "Pang B, Snyder MP. Systematic identification of silencers in human cells. "
    "Nat Genet 2020;52(3):254-263",
    "doi": "10.1038/s41588-020-0578-5",
    "curated_into": "RefSeq Functional Elements, reached through NCBI Gene esummary",
    "query": "ReSE[All Fields] AND human[orgn]",
    "licence": "NCBI places no restrictions on the use or distribution of the data it holds; "
    "RefSeq content created by the US government is not subject to copyright",
    "cost": "none. Free NCBI eutils, no API key. 0 AlphaGenome requests, 0 paid requests, no money.",
}

#: The matched control. 37.7% of named elements are already silencer_like, so a random set of scored
#: elements agrees at roughly that rate BY CONSTRUCTION and the whole check is worthless without a
#: control drawn the same way the ReSE set lands.
CONTROL = {
    "matched_on": "the (chromosome, cCRE class) pair, cell by cell: for each stratum the ReSE join "
    "occupies, the control draws the same NUMBER of elements from that same chromosome and that same "
    "cCRE class",
    "drawn_from": "scored elements carrying a named predicted target, which is the only population "
    "in which a direction exists to compare",
    "excludes": "every element the ReSE join itself used, so the control cannot re-draw the test set",
    "draws": 1000,
    "seed": 20261003,
    "why_class_matters": "the join is dELS and pELS only, and those two classes need not carry the "
    "genome-wide 37.7% between them; a control matched on chromosome alone would import whatever "
    "class composition the genome has",
}

#: Registered BEFORE any agreement rate exists. The reading is the DIFFERENCE in percentage points
#: between the ReSE set's silencer_like rate and the matched control's mean rate.
BANDS = {
    "CONSISTENT": "the ReSE set's rate exceeds the matched control's mean by >= 10 percentage points "
    "AND the 95% interval on the difference excludes 0",
    "DISPLACED_BUT_SMALL": "the 95% interval on the difference excludes 0 and the difference is "
    "under 10 percentage points in absolute value",
    "NO_SIGNAL": "the 95% interval on the difference includes 0",
    "INVERTED": "the ReSE set's rate falls below the matched control's mean by >= 10 percentage "
    "points AND the 95% interval on the difference excludes 0",
}
BAND_THRESHOLD_POINTS = 10.0
CONTROL_INTERVAL = "the 2.5th and 97.5th percentiles of the 1000 matched control rates"

WHAT_AGREEMENT_LICENSES = (
    "That the model's signed direction carries SOME information about assayed repressive ability, on "
    "dELS and pELS elements, in K562 and HepG2, at the rate measured and no better. It would license "
    "reporting silencer_like as a direction with one external point of contact where it previously "
    "had none."
)
WHAT_AGREEMENT_DOES_NOT_LICENSE = (
    "It is NOT a validation of the 230,839. The joined set is 0.5% of them, two cell lines of the "
    "four scored, and two cCRE classes. It does NOT make any element a silencer: molecular_role "
    "silencer is defined in genomeos/lang/grammar.py as 'a repressive element; never inferred from a "
    "direction of effect alone', and this check reads exactly a direction of effect. It does not "
    "attach a probability to any element: review R4 stands and nothing here converts a rate over a "
    "set into a certainty for a member of it. A LOW rate is NOT a refutation of the model either: a "
    "reporter fragment and a 1 Mb window argmax can disagree with neither being wrong."
)

TIE_BIAS = (
    "The model's own rule breaks a tie toward ACTIVATION ('ties go to activation, the usual case'), "
    "which biases this check AGAINST agreement on a silencer set, because every exact tie is counted "
    "as a disagreement. The bias is quantifiable and is quantified: the per-element committed results "
    "keep only the winning direction and its signed magnitude, but the AlphaGenome element cache "
    "keeps max_drop_log2fc and max_rise_log2fc for every gene in the window, so the exact-tie count "
    "is countable on the joined set rather than argued about."
)

FOUR_LINE_NEGATIVE = (
    "No reachable curated source covers the project's four scored cell lines. Every ReSE record "
    "validates in K562 or HepG2 and NOTHING validates in GM12878 or IMR-90 -- measured on the "
    "records themselves, not inferred. So any work planned as 'silencers, genome-wide, all four "
    "lines' is not fundable by any source that exists today, and that is the more durable finding of "
    "this lane whatever the agreement rate turns out to be."
)

FALSIFIERS = {
    "control_disagrees_with_itself": "the control procedure is run against a NULL test set -- 1000 "
    "further matched draws standing in for the ReSE set -- and if the difference between a drawn set "
    "and the control mean has a 95% interval excluding 0, the instrument manufactures signal from "
    "nothing and NO reading is published",
    "join_not_reproducible": "the joined element count recomputed from the committed per-element "
    "results differs from the count this registration states, in which case the join is not "
    "reproducible and the reading is withheld",
}

#: Keys a run would fill. A registration holding any of them, at any depth, is refused.
FORBIDDEN_IN_REGISTRATION = (
    "result",
    "verdict",
    "passed",
    "figures",
    "reading",
    "overall",
    "agreement_rate",
    "agreement_count",
    "band",
    "control_rate",
    "difference_points",
    "ties",
)


def rese_records() -> dict[str, dict[str, Any]]:
    """Every record the one NCBI query returned, from the local cache."""
    if not RESE_CACHE.exists():
        raise FileNotFoundError(f"{RESE_CACHE} is not fetched; run scripts/silenceragree_fetch.py")
    return json.loads(RESE_CACHE.read_text())["records"]


def rese_intervals(records: dict[str, dict[str, Any]]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """The joinable GRCh38 ReSE intervals, and an accounting of every record that is not one.

    The GRCh38 coordinates are `genomicinfo`'s. The description's are GRCh37 and are carried only so
    the lift is auditable; they are never joined against.
    """
    kept: list[dict[str, Any]] = []
    rejected = {"not_the_assayed_fragment": 0, "no_genomicinfo": 0, "unplaced_on_grch38": 0}
    builds: Counter[str] = Counter()
    length_equal = length_differs = relocated = 0
    for uid, rec in sorted(records.items()):
        hit = RESE_DESCRIPTION.match((rec.get("description") or "").strip())
        if not hit:
            rejected["not_the_assayed_fragment"] += 1
            continue
        build, chrom37, s37, e37 = hit.group(1), hit.group(2), int(hit.group(3)), int(hit.group(4))
        builds[build] += 1
        info = rec.get("genomicinfo") or []
        if not info:
            rejected["no_genomicinfo"] += 1
            continue
        loc = str(info[0].get("chrloc") or "")
        if not loc:
            rejected["unplaced_on_grch38"] += 1
            continue
        start = min(info[0]["chrstart"], info[0]["chrstop"])
        end = max(info[0]["chrstart"], info[0]["chrstop"])
        if loc != chrom37:
            relocated += 1
        if (end - start + 1) == (e37 - s37 + 1):
            length_equal += 1
        else:
            length_differs += 1
        clause = CELL_CLAUSE.search(rec.get("summary") or "")
        cells = sorted({c for c in RESE_CELL_LINES if clause and c in clause.group(1)})
        kept.append(
            {
                "uid": uid,
                "chrom": f"chr{loc}",
                "start": start,
                "end": end,
                "grch37": f"chr{chrom37}:{s37}-{e37}",
                "cells": cells,
            }
        )
    accounting = {
        "records_returned": len(records),
        "joinable_intervals": len(kept),
        "rejected": rejected,
        "rejected_reason": RESE_EXCLUDED_REASON,
        "build_tag_on_every_description": dict(builds),
        "grch38_span_equals_the_fragment_length": length_equal,
        "grch38_span_differs_in_length": length_differs,
        "relocated_to_another_chromosome_by_the_lift": relocated,
        "coordinates_joined_on": "genomicinfo, which is GRCh38 (NC_*.1x accessions)",
    }
    return kept, accounting


def element_index(chrom: str) -> list[tuple[int, int, str, bool]]:
    """`(start, end, id, has_predicted_target)` for one chromosome, sorted by start.

    NO DIRECTION IS RETURNED. This is the only door `scope()` has onto the element results, so a
    registration built through it cannot carry a direction even by mistake.
    """
    path = ELEMENTS_DIR / f"{chrom}.json"
    if not path.exists():
        return []
    rows = json.loads(path.read_text())
    return sorted((r["start"], r["end"], r["id"], r.get("predicted") is not None) for r in rows)


def ccre_classes(chrom: str) -> dict[str, str]:
    """Element id -> cCRE class, from the committed per-chromosome bed."""
    path = CCRE_DIR / f"ccres_{chrom}.bed.gz"
    if not path.exists():
        return {}
    out = {}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            field = line.rstrip("\n").split("\t")
            out[field[3]] = field[4]
    return out


def join(intervals: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Each ReSE fragment paired with the scored element it overlaps most, where one exists.

    Carries no direction: the pairing is by coordinate and the element is named by id. One row per
    ReSE fragment; a fragment overlapping several elements keeps the greatest overlap, and how often
    that happened is reported rather than hidden.
    """
    by_chrom: dict[str, list[dict[str, Any]]] = {}
    for item in intervals:
        by_chrom.setdefault(item["chrom"], []).append(item)
    rows: list[dict[str, Any]] = []
    counts = {
        "overlapping_a_scored_element": 0,
        "overlapping_one_with_a_named_target": 0,
        "no_scored_element": 0,
        "scored_but_no_gene_moved_enough": 0,
        "on_a_chromosome_with_no_scored_file": 0,
        "fragment_spanning_more_than_one_named_element": 0,
    }
    for chrom in sorted(by_chrom):
        index = element_index(chrom)
        if not index:
            counts["on_a_chromosome_with_no_scored_file"] += len(by_chrom[chrom])
            continue
        starts = [x[0] for x in index]
        classes = ccre_classes(chrom)
        for item in by_chrom[chrom]:
            s, e = item["start"], item["end"]
            hits = []
            for j in range(bisect.bisect_left(starts, s - 10_000), len(index)):
                est, een, eid, named = index[j]
                if est > e:
                    break
                if est < e and een > s:
                    hits.append((est, een, eid, named))
            if not hits:
                counts["no_scored_element"] += 1
                continue
            counts["overlapping_a_scored_element"] += 1
            named_hits = [h for h in hits if h[3]]
            if not named_hits:
                counts["scored_but_no_gene_moved_enough"] += 1
                continue
            counts["overlapping_one_with_a_named_target"] += 1
            if len(named_hits) > 1:
                counts["fragment_spanning_more_than_one_named_element"] += 1
            best = max(named_hits, key=lambda h: min(h[1], e) - max(h[0], s))
            rows.append(
                {
                    "uid": item["uid"],
                    "chrom": chrom,
                    "cells": item["cells"],
                    "element": best[2],
                    "el_class": classes.get(best[2], "unknown"),
                    "overlap_bp": min(best[1], e) - max(best[0], s),
                }
            )
    return rows, counts


def base_rate() -> dict[str, int | float]:
    """The genome-wide silencer_like level the control has to beat, from the committed results.

    This reads the direction field and is therefore NOT part of `scope()`. It is the one quantity the
    registration needs that is a direction count, it is a figure the earlier probe already published,
    and it is the test set's base rate rather than its agreement rate.
    """
    total = named = repress = 0
    for path in sorted(ELEMENTS_DIR.glob("chr*.json")):
        for row in json.loads(path.read_text()):
            total += 1
            pred = row.get("predicted")
            if pred:
                named += 1
                repress += pred["action"] == "represses"
    return {
        "elements_scored": total,
        "with_a_named_predicted_target": named,
        "silencer_like": repress,
        "rate_of_named": round(repress / named, 5) if named else 0.0,
    }


def scope() -> dict[str, Any]:
    """What this study will measure, with no direction read anywhere in it.

    Everything here comes from the ReSE records and from `element_index()`, which returns no action
    field. The function's whole output is a denominator and its composition.
    """
    kept, accounting = rese_intervals(rese_records())
    rows, counts = join(kept)
    strata = Counter((r["chrom"], r["el_class"]) for r in rows)
    return {
        "rese": accounting,
        "join": counts,
        "denominator": len(rows),
        "by_element_class": dict(Counter(r["el_class"] for r in rows).most_common()),
        "by_chromosome": dict(sorted(Counter(r["chrom"] for r in rows).items())),
        "by_validation_cell_line": dict(
            Counter("+".join(r["cells"]) or "(cell line not parsed)" for r in rows).most_common()
        ),
        "strata": {f"{c}|{k}": n for (c, k), n in sorted(strata.items())},
        "strata_count": len(strata),
        "distinct_elements": len({r["element"] for r in rows}),
    }


def control_draws(rows: list[dict[str, Any]], draws: int, seed: int) -> list[list[str]]:
    """`draws` control sets, each matched to `rows` stratum by stratum on (chromosome, cCRE class)."""
    need = Counter((r["chrom"], r["el_class"]) for r in rows)
    used = {r["element"] for r in rows}
    pools: dict[tuple[str, str], list[str]] = {}
    for chrom in sorted({c for c, _ in need}):
        classes = ccre_classes(chrom)
        for start, end, eid, named in element_index(chrom):  # noqa: B007
            if not named:
                continue
            key = (chrom, classes.get(eid, "unknown"))
            if key in need and eid not in used:
                pools.setdefault(key, []).append(eid)
    rng = random.Random(seed)
    out = []
    for _ in range(draws):
        picked: list[str] = []
        for key, n in sorted(need.items()):
            pool = pools.get(key, [])
            if len(pool) < n:
                raise ValueError(f"stratum {key} needs {n} control elements and the pool holds {len(pool)}")
            picked.extend(rng.sample(pool, n))
        out.append(picked)
    return out


#: The second reading, registered beside the first. R1 is the call AS SHIPPED: `predicted.action`,
#: the argmax over all 371 tracks, most of which are neither K562 nor HepG2. R2 holds the same gene
#: and asks that ONE cell line's own track, which is the comparison ReSE actually licenses, since a
#: fragment validated in K562 says nothing about a direction argmaxed in occipital lobe.
READING_1 = (
    "AS SHIPPED: predicted.action == 'represses' on the committed per-element result, the argmax over "
    "every track. This is the quantity 230,839 elements carry and the one a reader of the project "
    "would meet."
)
READING_2 = (
    "CELL-MATCHED: the sign of predicted_by_cell[<the ReSE record's own validation cell line>] for "
    "that same predicted gene -- positive means the gene rose on deletion, i.e. silencer-like in "
    "that cell. No magnitude threshold is applied, because MIN_EFFECT was already spent selecting "
    "the gene; a value of exactly 0.0 counts as NOT silencer-like, which is the model's own tie rule "
    "carried over unchanged. A fragment validated in both lines is counted silencer-like if either "
    "line's value is positive, and how many such fragments there are is stated. The count whose "
    "|value| < 0.01 -- effectively no movement in the matched cell at all -- is reported beside the "
    "rate rather than removed from it."
)
READING_3 = (
    "DESCRIPTIVE, no band and no control: the exact-tie count on the joined set, read from the "
    "AlphaGenome element cache's per-gene max_drop_log2fc and max_rise_log2fc. An exact tie is "
    "-max_drop == max_rise on the gene the rule selected, and every one of them was awarded to "
    "activation. This quantifies the tie bias instead of arguing about it."
)
NEAR_ZERO = 0.01
CACHE_DIR = ROOT / "data/knowledge/alphagenome/elements"


#: The overlap predicate and the tie-break, written out so the population is reproducible from the
#: registration alone rather than from reading the code.
JOIN_RULE = {
    "overlap": "a ReSE fragment [rs, re) and a scored element [es, ee) overlap when es < re AND "
    "ee > rs. Touching ends do not overlap. No padding and no minimum overlap: one shared base is "
    "an overlap, and the overlap in bases is recorded on every row so a stricter rule can be applied "
    "afterwards without re-joining.",
    "one_row_per_fragment": "the unit is the ReSE FRAGMENT, not the element, so a fragment never "
    "contributes more than one row and a fragment cannot be counted twice by overlapping two "
    "elements.",
    "fragment_overlapping_several_elements": "the element with the GREATEST overlap in bases wins, "
    "ties by that measure broken by the order the elements come off the sorted index, which is "
    "deterministic. How many fragments were in this position is reported as "
    "`fragment_spanning_more_than_one_named_element` and not hidden.",
    "element_overlapping_several_fragments": "an element CAN be named by more than one fragment, in "
    "which case it contributes a row per fragment. `distinct_elements` is reported beside the "
    "denominator so the reader can see whether this happened.",
    "eligibility": "the element must carry a NAMED predicted target. An element on which no gene "
    "moved by MIN_EFFECT has no direction at all and is counted in the join block, never scored as a "
    "disagreement.",
}

#: The two base rates the coordinator required stated in advance.
BASE_RATES = {
    "curated_label_share_among_the_joined_elements": "100% BY CONSTRUCTION, and saying so is the "
    "point. Every joined element overlaps an assay-validated ReSE silencer, because that is what the "
    "join selects on; this base rate is degenerate and carries no information, and quoting it as "
    "though it did would be the error. It is stated so that no reader mistakes the agreement rate "
    "for a comparison against it.",
    "silencer_like_share_measured_on_the_joined_set_strata": "the informative one. NOT the global "
    "37.7%: the exact share of named-target elements that are silencer_like within the (chromosome, "
    "cCRE class) strata the join occupies, computed over the WHOLE stratum population and not a "
    "sample, with the joined elements themselves excluded so the base rate cannot contain the "
    "observation. This is the number the agreement rate is compared against, and the global figure "
    "is carried only to show how far the strata sit from it.",
}

#: Prior exposure, checked against the committed record rather than assumed.
PRIOR_EXPOSURE = {
    "what_the_probe_lane_committed": "commit b15bbae, a docs/ROADMAP.md row and nothing else. It "
    "wrote no file into data/, its board entry data/work/lane-silencer.json carries note: null and "
    "files: [], and no result, script or artefact of it exists in the tree.",
    "what_that_row_records_it_saw": [
        "the global silencer_like totals, 230,839 of 612,323 named elements across 75 results",
        "4,998 ReSE records, and 1,864 of 4,992 overlapping a project cCRE",
        "the cCRE class breakdown of that overlap, dELS 1,237",
        "1,585 overlapping an already-scored element",
        "SilencerDB's validated tier at offset exactly (0,0) over all 4,998",
        "GM12878 0 validated silencers of 23,894, and IMR90 1",
    ],
    "what_it_states_about_the_agreement_rate": "in its own row: it 'deliberately did not compute an "
    "agreement rate', leaving that to a registered lane.",
    "so_is_this_test_blind": "BLIND BY ATTESTATION, NOT BLIND BY AUDIT, and the distinction is "
    "registered rather than smoothed over. The committed row shows the join count, the class "
    "breakdown and the global silencer_like total -- no cross-tabulation of the curated label "
    "against the direction appears anywhere in it. But because that lane committed no artefact, "
    "there is nothing in the tree that could CONFIRM OR REFUTE whether it computed one and did not "
    "report it. The absence rests on the lane's own statement. A reader who wants a stronger claim "
    "than attestation cannot have it from this record, and this registration does not offer one.",
    "what_this_lane_had_seen_before_this_registration_was_written": "the ReSE records and their "
    "coordinates; the join counts and strata, computed through scope(), which is handed no action "
    "field; the genome-wide base rate 230,839/612,323, which the probe had already published. NO "
    "agreement count, NO cell-matched count and NO control rate existed in this lane before the "
    "registration was committed, and the code is split so that scope() could not have produced one.",
}

#: Independence, which is the scientifically decisive question and which is NOT established.
INDEPENDENCE = {
    "status": "NOT ESTABLISHED. This check is registered as a comparison against an external "
    "SOURCE, and NOT as a test against data the model has provably never seen.",
    "what_can_be_established_offline": [
        "The model predicts on K562 and HepG2 tracks directly: every scored element carries "
        "predicted_by_cell values for HepG2, IMR-90, K562 and GM12878, so functional-genomics data "
        "from the two ReSE cell lines is inside the corpus the model was fitted to.",
        "ReSE's candidate regions were chosen as ACCESSIBLE CHROMATIN in those same cell lines, so "
        "the regions the assay tested are regions whose accessibility the model was trained on. The "
        "element coordinates are therefore NOT independent of the training data in any useful sense.",
        "The assay's own readout -- cell survival under reporter repression, log-scale enrichment of "
        "a lentiviral library -- is not among the model's output modalities (RNA-seq, ATAC, DNase, "
        "ChIP, contact maps, splicing), so the ReSE LABEL is not a quantity the model was fitted to "
        "predict.",
    ],
    "what_cannot_be_established": "whether the ReSE screen's data, or any derivative of it, entered "
    "the model's training corpus. The training set is not published at the resolution that would "
    "answer this, and nothing in this checkout can answer it. It is registered as unknown rather "
    "than assumed favourable.",
    "so_the_word_independent": "'independent' is used here of the SOURCE -- a curated external assay "
    "the project had never consulted -- and NOT of the training data. Any agreement found is an "
    "agreement between a model and a label that may or may not be downstream of something it saw. "
    "The project already registers this class of caveat and does not soften it: lane-labelgene's "
    "GTEx arm carries 'AlphaGenome's GTEx tracks were trained on this same RNA-seq, so this arm may "
    "be recall of seen data rather than biology and cannot distinguish them.' The same sentence "
    "applies here with the training status unknown instead of known.",
    "why_this_matters_more_than_the_rate": "overstating independence would be worse than finding "
    "nothing, because the number would be quoted. A check whose independence is assumed is not a "
    "check.",
}

#: The nearest prior art in the project, so this reading is not presented as the first of its kind.
PRIOR_ART = (
    "The project has already compared this model to a reporter assay, on the ACTIVATION side: the "
    "ENCODE4 joint lentiMPRA library, 51,376 non-overlapping elements in K562, HepG2 and WTC11 "
    "(docs/ALPHAGENOME.md, 2026-09-12). Predicted DNase correlates with measured reporter activity "
    "at rank 0.443 in K562 and 0.323 in HepG2 in the matched cell, and 0.030 and 0.119 in the other "
    "cell. That result's own conclusion names the limit this study inherits: a distal enhancer-like "
    "element is barely more active than a sequence in no class at all, which is 'the known limit of "
    "a reporter that measures a sequence out of its node'. The repressive side has never been looked "
    "at, which is what this study adds."
)

READING_1 = (
    "AS SHIPPED: predicted.action == 'represses' on the committed per-element result, the argmax over "
    "every track. This is the quantity 230,839 elements carry and the one a reader of the project "
    "would meet."
)
READING_2 = (
    "CELL-MATCHED: the sign of predicted_by_cell[<the ReSE record's own validation cell line>] for "
    "that same predicted gene -- positive means the gene rose on deletion, i.e. silencer-like in "
    "that cell. No magnitude threshold is applied, because MIN_EFFECT was already spent selecting "
    "the gene; a value of exactly 0.0 counts as NOT silencer-like, which is the model's own tie rule "
    "carried over unchanged. A fragment validated in both lines is counted silencer-like if either "
    "line's value is positive, and how many such fragments there are is stated. The count whose "
    "|value| < 0.01 -- effectively no movement in the matched cell at all -- is reported beside the "
    "rate rather than removed from it."
)
READING_3 = (
    "DESCRIPTIVE, no band and no control: the exact-tie count on the joined set, read from the "
    "AlphaGenome element cache's per-gene max_drop_log2fc and max_rise_log2fc. An exact tie is "
    "-max_drop == max_rise on the gene the rule selected, and every one of them was awarded to "
    "activation. This quantifies the tie bias instead of arguing about it."
)
NEAR_ZERO = 0.01
CACHE_DIR = ROOT / "data/knowledge/alphagenome/elements"

#: The interval is CLUSTERED, because elements near one another are not independent draws and a
#: binomial interval over neighbouring elements would be confidently wrong.
INTERVAL = {
    "method": "a cluster bootstrap over CHROMOSOMES: chromosomes are resampled with replacement, "
    "and within each draw the ReSE rate and the matched control rate are recomputed over only the "
    "chromosomes drawn, so the difference carries the between-chromosome variation.",
    "clusters": "chromosomes carrying at least one joined element",
    "minimum_clusters": 10,
    "if_fewer": "NO INTERVAL IS REPORTED and the reading is INCONCLUSIVE_TOO_FEW_CLUSTERS. A naive "
    "interval is not offered in its place.",
    "draws": 2000,
    "seed": 20261003,
    "why_not_binomial": "1,181 joined elements are not 1,181 independent draws. Many sit in the same "
    "locus, several are named by fragments from the same ReSE region, and silencer_like is itself "
    "spatially correlated because neighbouring elements share the genes in their windows. A "
    "binomial interval would be narrow and wrong.",
}

BANDS = {
    "INCONCLUSIVE_TOO_FEW_CLUSTERS": "fewer than 10 chromosomes carry a joined element, so no "
    "interval is computed and no band is read",
    "INCONCLUSIVE_UNDERPOWERED": "the interval includes 0 AND reaches beyond the band threshold, so "
    "'no signal' and 'consistent' are both inside it and this population cannot tell them apart. "
    "THIS IS A COMPLETE READING: it says the check is underpowered at this denominator.",
    "NO_SIGNAL": "the interval includes 0 and stays inside the band threshold, so the difference is "
    "both unsigned and small",
    "CONSISTENT": "the difference exceeds the stratum-matched base rate by >= 10 percentage points "
    "and the interval excludes 0",
    "INVERTED": "the difference falls below the stratum-matched base rate by >= 10 percentage points "
    "and the interval excludes 0",
    "DISPLACED_BUT_SMALL": "the interval excludes 0 and the difference is under 10 percentage points "
    "in absolute value",
}
BAND_THRESHOLD_POINTS = 10.0
BANDS_ARE_RELATIVE = (
    "Every band is on the DIFFERENCE in percentage points between the observed rate and the "
    "stratum-matched base rate, never on the observed rate itself. An absolute band at a base rate "
    "near 37.7% would be near-vacuous."
)
POWER = {
    "stated_in_advance": True,
    "how": "before the real agreement rate is computed, the whole band machinery is run on SYNTHETIC "
    "test sets built from the control: a matched control draw is taken as the test set and a chosen "
    "number of its non-silencer_like members is flipped to silencer_like to raise its rate by a "
    "known shift. The band returned is recorded over many repeats, at shifts of 0, 5, 10, 15 and 20 "
    "percentage points.",
    "repeats_per_shift": 200,
    "what_it_answers": "whether 1,181 elements clustered in 23 chromosomes can distinguish the "
    "registered bands at all. If the 0-point shift does not read NO_SIGNAL most of the time, or the "
    "15-point shift does not read CONSISTENT most of the time, then the check is underpowered and "
    "THAT IS THE FINDING.",
    "reported_whatever_it_says": "the power table is published beside the reading, not consulted and "
    "discarded.",
}
CONTROL_INTERVAL = "the 2.5th and 97.5th percentiles of the clustered bootstrap on the difference"

FALSIFIERS = {
    "control_disagrees_with_itself": "the control procedure is run against a NULL test set -- a "
    "further matched draw standing in for the ReSE set -- and if the clustered interval on that "
    "difference excludes 0, the instrument manufactures signal from nothing and NO reading is "
    "published",
    "join_not_reproducible": "the joined element count recomputed from the committed per-element "
    "results differs from the count this registration states, in which case the join is not "
    "reproducible and the reading is withheld",
    "power_absent": "if the registered power run cannot reach CONSISTENT at a 15-point shift in a "
    "majority of repeats, the reading is published as INCONCLUSIVE_UNDERPOWERED whatever the "
    "observed difference is",
}


def classify(
    difference_points: float,
    low: float,
    high: float,
    clusters: int,
    min_clusters: int = 10,
    threshold: float = BAND_THRESHOLD_POINTS,
) -> str:
    """The registered band. Relative to the base rate, clustered, and with power inside the rule.

    The order of the tests is the rule and nothing else decides. `low` and `high` are the clustered
    interval's ends on the DIFFERENCE, in percentage points.
    """
    if clusters < min_clusters:
        return "INCONCLUSIVE_TOO_FEW_CLUSTERS"
    if low <= 0.0 <= high:
        if high >= threshold or low <= -threshold:
            return "INCONCLUSIVE_UNDERPOWERED"
        return "NO_SIGNAL"
    if difference_points >= threshold:
        return "CONSISTENT"
    if difference_points <= -threshold:
        return "INVERTED"
    return "DISPLACED_BUT_SMALL"


def stratum_matched_base_rate(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """silencer_like's EXACT share within the strata the join occupies, joined elements excluded.

    Over the whole stratum population, not a sample, so it carries no sampling error of its own.
    This is the base rate the bands are relative to.
    """
    need = Counter((r["chrom"], r["el_class"]) for r in rows)
    used = {r["element"] for r in rows}
    hits = total = 0
    per: dict[str, dict[str, int]] = {}
    for chrom in sorted({c for c, _ in need}):
        classes = ccre_classes(chrom)
        table = directions(chrom)
        for eid, d in table.items():
            if eid in used:
                continue
            if (chrom, classes.get(eid, "unknown")) not in need:
                continue
            total += 1
            hit = d["action"] == "represses"
            hits += hit
            slot = per.setdefault(f"{chrom}|{classes.get(eid, 'unknown')}", {"n": 0, "silencer_like": 0})
            slot["n"] += 1
            slot["silencer_like"] += hit
    return {
        "population": "every named-target element in the join's (chromosome, cCRE class) strata, "
        "with the joined elements themselves excluded",
        "elements": total,
        "silencer_like": hits,
        "percent": _rate(hits, total),
        "per_stratum": per,
    }


def directions(chrom: str) -> dict[str, dict[str, Any]]:
    """Element id -> the direction fields. The ONLY reader of an action in this module's run path."""
    path = ELEMENTS_DIR / f"{chrom}.json"
    if not path.exists():
        return {}
    out = {}
    for row in json.loads(path.read_text()):
        pred = row.get("predicted")
        if pred:
            out[row["id"]] = {
                "action": pred["action"],
                "gene": pred["gene"],
                "log2fc": pred["log2_fold_change"],
                "by_cell": row.get("predicted_by_cell") or {},
            }
    return out


def _rate(hits: int, n: int) -> float:
    return round(100.0 * hits / n, 3) if n else 0.0


def clustered_interval(
    per_cluster: dict[str, tuple[int, int]],
    base_percent: float,
    draws: int,
    seed: int,
) -> dict[str, Any]:
    """A cluster bootstrap over chromosomes on the DIFFERENCE from the base rate.

    `per_cluster` maps a chromosome to `(silencer_like, n)` for the set being measured. Each draw
    resamples chromosomes with replacement and recomputes the pooled rate over the drawn
    chromosomes; the statistic is that rate minus `base_percent`.
    """
    names = sorted(per_cluster)
    if len(names) < INTERVAL["minimum_clusters"]:
        return {"clusters": len(names), "interval_percent_points": None, "reason": INTERVAL["if_fewer"]}
    rng = random.Random(seed)
    stats = []
    for _ in range(draws):
        drawn = [names[rng.randrange(len(names))] for _ in names]
        hits = sum(per_cluster[c][0] for c in drawn)
        n = sum(per_cluster[c][1] for c in drawn)
        if n:
            stats.append(_rate(hits, n) - base_percent)
    stats.sort()
    return {
        "clusters": len(names),
        "draws": len(stats),
        "interval_percent_points": [
            round(stats[int(0.025 * len(stats))], 3),
            round(stats[int(0.975 * len(stats))], 3),
        ],
    }


def power_table(
    rows: list[dict[str, Any]], base: float, shifts: tuple[int, ...] = (0, 5, 10, 15, 20)
) -> dict:
    """What bands this denominator can reach, computed BEFORE the real rate is read.

    A matched control draw stands in for the ReSE set, so it carries the base rate; a chosen number
    of its non-silencer_like members is flipped to raise its rate by a known shift; and the full band
    machinery is run on it. The table says whether 1,181 elements in 23 clusters can tell the
    registered bands apart at all.
    """
    chroms = sorted({r["chrom"] for r in rows})
    by_id: dict[str, dict[str, Any]] = {}
    for c in chroms:
        by_id.update(directions(c))
    order = [(r["chrom"], r["element"]) for r in rows]
    sets = control_draws(rows, POWER["repeats_per_shift"], CONTROL["seed"] + 7)
    out: dict[str, Any] = {}
    for shift in shifts:
        bands: Counter[str] = Counter()
        for repeat, drawn in enumerate(sets):
            per: dict[str, list[int]] = {}
            flagged = []
            for i, eid in enumerate(drawn):
                chrom = order[i][0]
                hit = (by_id.get(eid) or {}).get("action") == "represses"
                per.setdefault(chrom, [0, 0])
                per[chrom][1] += 1
                if hit:
                    per[chrom][0] += 1
                else:
                    flagged.append(chrom)
            want = int(round(len(drawn) * shift / 100.0))
            rng = random.Random(CONTROL["seed"] + 7 + repeat)
            for chrom in rng.sample(flagged, min(want, len(flagged))):
                per[chrom][0] += 1
            cl = {c: (v[0], v[1]) for c, v in per.items()}
            hits = sum(v[0] for v in cl.values())
            n = sum(v[1] for v in cl.values())
            iv = clustered_interval(cl, base, 400, CONTROL["seed"] + 7 + repeat)
            ends = iv["interval_percent_points"]
            if ends is None:
                bands["INCONCLUSIVE_TOO_FEW_CLUSTERS"] += 1
                continue
            bands[classify(_rate(hits, n) - base, ends[0], ends[1], iv["clusters"])] += 1
        out[f"shift_{shift}_points"] = dict(bands.most_common())
    return out


def exact_ties(
    rows: list[dict[str, Any]], by_chrom_dirs: dict[str, dict[str, dict[str, Any]]]
) -> dict[str, Any]:
    """How many of the joined elements' selected genes were an EXACT tie awarded to activation."""
    seen = ties = missing = 0
    margins: list[float] = []
    for chrom in sorted({r["chrom"] for r in rows}):
        path = CACHE_DIR / f"{chrom}.json.gz"
        if not path.exists():
            missing += sum(1 for r in rows if r["chrom"] == chrom)
            continue
        with gzip.open(path, "rt") as handle:
            cache = json.load(handle)
        for row in rows:
            if row["chrom"] != chrom:
                continue
            entry = cache.get(row["element"])
            want = by_chrom_dirs[chrom].get(row["element"], {}).get("gene")
            gene = next((g for g in (entry or {}).get("genes", []) if g["gene"] == want), None)
            if gene is None:
                missing += 1
                continue
            seen += 1
            drop, rise = -gene["max_drop_log2fc"], gene["max_rise_log2fc"]
            if drop == rise:
                ties += 1
            margins.append(round(rise - drop, 6))
    return {
        "what_a_tie_is": "-max_drop_log2fc == max_rise_log2fc on the gene the rule selected, which "
        "the model's own rule awards to activation",
        "selected_genes_found_in_the_cache": seen,
        "not_found_in_the_cache": missing,
        "exact_ties_awarded_to_activation": ties,
        "tie_share_of_those_checked_percent": _rate(ties, seen),
        "margins_within_0.001_of_a_tie": sum(1 for m in margins if abs(m) <= 0.001),
        "direction_of_the_bias": "every tie is counted as a DISAGREEMENT by this study, so the tie "
        "rule biases the measured agreement DOWNWARD on a silencer set",
    }


def readings() -> dict[str, Any]:
    """The measurement. NOT called before the registration is committed."""
    kept, accounting = rese_intervals(rese_records())
    rows, counts = join(kept)
    chroms = sorted({r["chrom"] for r in rows})
    dirs = {c: directions(c) for c in chroms}
    base = stratum_matched_base_rate(rows)
    base_pc = base["percent"]

    r1 = {c: [0, 0] for c in chroms}
    r2 = {c: [0, 0] for c in chroms}
    r2_zero = r2_near = r2_either = 0
    for row in rows:
        d = dirs[row["chrom"]][row["element"]]
        r1[row["chrom"]][1] += 1
        if d["action"] == "represses":
            r1[row["chrom"]][0] += 1
        cells = row["cells"] or list(RESE_CELL_LINES)
        values = [d["by_cell"][c] for c in cells if c in d["by_cell"]]
        if values:
            r2[row["chrom"]][1] += 1
            if any(v > 0 for v in values):
                r2[row["chrom"]][0] += 1
                if len(cells) > 1:
                    r2_either += 1
            if all(v == 0 for v in values):
                r2_zero += 1
            if all(abs(v) < NEAR_ZERO for v in values):
                r2_near += 1

    # The cell-matched base rate has to be matched too: the same cell line per stratum member.
    cell_of = {r["element"]: (r["cells"] or list(RESE_CELL_LINES)) for r in rows}
    order = [r["element"] for r in rows]
    by_id: dict[str, dict[str, Any]] = {}
    for c in chroms:
        by_id.update(dirs[c])
    sets = control_draws(rows, CONTROL["draws"], CONTROL["seed"])
    c2: list[float] = []
    for drawn in sets:
        hit = 0
        for i, eid in enumerate(drawn):
            d = by_id.get(eid)
            if d is None:
                continue
            values = [d["by_cell"][c] for c in cell_of[order[i]] if c in d["by_cell"]]
            if values and any(v > 0 for v in values):
                hit += 1
        c2.append(_rate(hit, len(drawn)))
    base_pc_2 = round(sum(c2) / len(c2), 3)

    def read(per: dict[str, list[int]], basis: float, label: str) -> dict[str, Any]:
        cl = {c: (v[0], v[1]) for c, v in per.items() if v[1]}
        hits = sum(v[0] for v in cl.values())
        n = sum(v[1] for v in cl.values())
        obs = _rate(hits, n)
        iv = clustered_interval(cl, basis, INTERVAL["draws"], INTERVAL["seed"])
        ends = iv["interval_percent_points"]
        band = (
            classify(obs - basis, ends[0], ends[1], iv["clusters"])
            if ends
            else "INCONCLUSIVE_TOO_FEW_CLUSTERS"
        )
        return {
            "what": label,
            "denominator": n,
            "silencer_like": hits,
            "observed_percent": obs,
            "base_rate_percent": basis,
            "difference_points": round(obs - basis, 3),
            "clustered_interval_on_the_difference": iv,
            "band": band,
            "per_chromosome": {c: {"silencer_like": v[0], "n": v[1]} for c, v in sorted(cl.items())},
        }

    # The falsifier: the whole machinery against a NULL test set drawn the same way.
    null = control_draws(rows, 1, CONTROL["seed"] + 1)[0]
    null_per: dict[str, list[int]] = {}
    for i, eid in enumerate(null):
        chrom = rows[i]["chrom"]
        null_per.setdefault(chrom, [0, 0])
        null_per[chrom][1] += 1
        if (by_id.get(eid) or {}).get("action") == "represses":
            null_per[chrom][0] += 1
    null_read = read(null_per, base_pc, "a matched control draw standing in for the ReSE set")

    return {
        "denominator": len(rows),
        "rese": accounting,
        "join": counts,
        "join_rule": JOIN_RULE,
        "base_rates": {
            "curated_label_share_among_the_joined_elements_percent": 100.0,
            "curated_label_share_is_degenerate": BASE_RATES[
                "curated_label_share_among_the_joined_elements"
            ],
            "stratum_matched_as_shipped": base,
            "stratum_matched_cell_matched_percent": base_pc_2,
            "global_for_comparison": base_rate(),
        },
        "reading_1_as_shipped": read(r1, base_pc, READING_1),
        "reading_2_cell_matched": {
            **read(r2, base_pc_2, READING_2),
            "exactly_zero_in_the_matched_cell": r2_zero,
            "absolute_value_under_0.01_in_the_matched_cell": r2_near,
            "counted_on_either_line_because_validated_in_both": r2_either,
        },
        "reading_3_tie_bias": exact_ties(rows, dirs),
        "power": power_table(rows, base_pc),
        "falsifier_control_disagrees_with_itself": {
            "fired": null_read["band"] not in ("NO_SIGNAL", "INCONCLUSIVE_UNDERPOWERED"),
            "null_reading": null_read,
        },
    }
