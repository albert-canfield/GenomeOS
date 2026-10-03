# SPDX-License-Identifier: AGPL-3.0-or-later
"""The remaining one-target consumers, moved and held against their registered falsifiers.

    uv run --frozen python scripts/onetarget2_run.py [--chrom chr21] [--modules closure,...]

Writes data/results/onetarget2.json, the result the registration at 24adf33 / aed8ae9 named. Each
arm is ONE module of that registration's census, read twice -- as it reads today through the compact
one-target head, and through `attribution.targets.window_reading`, the one shared window reader --
and judged against THAT MODULE'S OWN falsifier, never against the set's.

The file never implies completeness: `modules_moved` and `modules_outstanding` are both written, with
the outstanding ones carrying the reason they are outstanding (a live peer lane holds the path, or the
file is inside the paid study's frozen 52-file import closure).

Before any arm is read, `onetarget2.check_all()` runs and `check_head_invariant()` is applied to every
arm's census: ONE any-gene head disagreement refutes the invariant the whole wave stands on and stops
the run, rather than being written as a result.

0 AlphaGenome requests, no network, no money, no key read: the per-element response cache is already
on this machine, and one chromosome's archive is held at a time.
"""

from __future__ import annotations

import argparse
import sys
import time
from datetime import date
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import closure as cl  # noqa: E402
from genomeos.attribution import onetarget2 as ot  # noqa: E402
from genomeos.attribution.targets import ELEMENT_CACHE, ElementResponses  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "onetarget2"
REFERENCE = Path("data/reference")

OWN_CODE = (
    "genomeos/attribution/onetarget2.py",
    "genomeos/attribution/targets.py",
    "genomeos/attribution/closure.py",
    "genomeos/attribution/measured.py",
    "genomeos/attribution/not_open_profile.py",
    "scripts/onetarget2_register.py",
    "scripts/onetarget2_run.py",
    "tests/test_onetarget2.py",
    "tests/test_closure_window.py",
    "tests/test_measured_window.py",
    "tests/test_not_open_profile_window.py",
)


#: The stores this run reads, enumerated DELIBERATELY and not taken from the tracer. Declaring what
#: `manifest.traced_reads()` reports would make the contract circular -- the writer would declare
#: whatever it happened to open and the reconciliation could never fail. So the groups below are
#: written out from what the arms are known to read, and if the enumeration is wrong the write is
#: refused, which is the check doing its job. It refused this run once, naming 155 undeclared files.
RESPONSE_CACHE = Path("data/knowledge/alphagenome/elements")
COMPACT_TABLES = Path("data/knowledge/alphagenome/all_elements")
RESULTS = Path("data/results")
RUN_SUMMARIES = ("enhancer_targets_all", "constrained_targets", "enhancer_targets")
#: The measured layer's own assay sources, read by `measured.Layer.load`.
ASSAY_SOURCES = (
    "data/knowledge/crispri/EPCrisprBenchmark_combined_data.training_K562.GRCh38.tsv.gz",
    "data/knowledge/crispri/EPCrisprBenchmark_combined_data.heldout_5_cell_types.GRCh38.tsv.gz",
    "data/knowledge/mpra/ENCFF475FKV.bed.gz",
    "data/knowledge/mpra/ENCFF769REH.bed.gz",
    "data/knowledge/mpra/ENCFF802FUV.bed.gz",
    "data/knowledge/satmut/elements.tsv.gz",
    "data/knowledge/vista/locus.tsv.gz",
)
#: Read by `not_open_profile` through `compile._ccres` and `compile._Interspersed`, chr21 only
#: because that is the chromosome its arm runs on.
PROFILE_AXES = ("data/results/ccres_chr21.bed.gz", "data/results/rmsk_chr21.bed.gz")
#: `pilot_bio.gene_tss` reads a content-addressed cache whose file NAME is a digest, so it cannot be
#: predicted; the directory is declared as a tree instead.
TSS_CACHE = Path("data/cache/holdout/genes")


def inputs(chroms: tuple[str, ...], profile_chrom: str) -> list[dict[str, Any]]:
    """Every file this run reads, declared by group and hashed member by member."""
    here = Path(__file__).resolve().parents[1]
    out: list[dict[str, Any]] = []

    def group(label: str, paths: list[Path], why: str) -> None:
        present = [p for p in paths if (here / p).exists()]
        if present:
            out.append(mf.files_entry(label, present, partition=why))

    group(
        "the sweep's per-element response cache, one archive per chromosome",
        [RESPONSE_CACHE / f"{c}.json.gz" for c in chroms],
        "the window: every gene the sweep scored at an element, read at threshold=0.0. One archive "
        "is held at a time, so the memory bill is the largest chromosome and not the tree",
    )
    group(
        "the compact one-target tables, one per chromosome",
        [COMPACT_TABLES / f"{c}.json" for c in chroms],
        "the head each consumer reads today, reached through the run summaries' elements_where field",
    )
    group(
        "the deletion runs' committed summaries",
        [RESULTS / f"{n}_{c}.json" for n in RUN_SUMMARIES for c in chroms],
        "targets.attributed reads all three in order and each one's elements_where points at the "
        "compact table above",
    )
    group(
        "GENCODE v50, per chromosome",
        [Path(f"data/reference/gencode_v50_{c}.gff3.gz") for c in chroms],
        "the protein-coding symbol set that forms the coding window, and the TSS of each gene a "
        "window rule's distance is measured to",
    )
    group(
        "the measured layer's assay sources",
        [Path(x) for x in ASSAY_SOURCES],
        "CRISPRi, lentiMPRA, saturation mutagenesis and VISTA as measured.Layer.load reads them; "
        "the CRISPRi genes_regulated list is what measured.py's falsifier is asked about",
    )
    group(
        f"the element axes for {profile_chrom}",
        [Path(x) for x in PROFILE_AXES],
        "cCREs and interspersed repeats, read by not_open_profile through compile's own loaders to "
        "type each rule's element",
    )
    if (here / TSS_CACHE).exists():
        out.append(
            mf.input_entry(
                here / TSS_CACHE,
                partition="pilot_bio.gene_tss's content-addressed cache. Declared as a TREE because "
                "each file is named by its own digest, so the name cannot be written down in advance",
            )
        )
    inv = here / ot.PRIOR_EVIDENCE
    if inv.exists():
        out.append(
            mf.input_entry(
                inv,
                partition="the committed 2026-09-27 result, read by check_the_quoted_invariant to "
                "verify the invariant's quoted figures against the file rather than recall them",
            )
        )
    return out


def coding_of(chrom: str) -> set[str]:
    """The chromosome's protein-coding symbols, from the per-chromosome GENCODE v50 file.

    The same source `scripts/onetarget_consumers.py` used on 2026-09-27, so the coding set this
    wave reads is the one the earlier one-target lanes read.
    """
    from genomeos.genome import Annotation

    ann = Annotation.from_gff3(REFERENCE / f"gencode_v50_{chrom}.gff3.gz", {chrom})
    return {g.symbol for g in ann.genes.values() if g.type == "protein_coding"}


def arm_closure(chrom: str, responses: ElementResponses, coding: set[str]) -> dict[str, Any]:
    """closure.py. Falsifier, as registered: a gene whose closure gains or loses an element when the
    elements are grouped by every gene at the bar instead of by the head alone."""
    t = time.time()
    old = cl.attributed_elements(chrom)
    new, census = cl.window_elements(chrom, responses, coding)
    ot.check_head_invariant(census, f"closure.window_elements {chrom}")
    gained = {g: sorted({x["id"] for x in new[g]} - {x["id"] for x in old.get(g, [])}) for g in new}
    gained = {g: v for g, v in gained.items() if v}
    lost = {g: sorted({x["id"] for x in old[g]} - {x["id"] for x in new.get(g, [])}) for g in old}
    lost = {g: v for g, v in lost.items() if v}
    new_genes = sorted(g for g in new if g not in old)
    return {
        "module": "genomeos/attribution/closure.py",
        "falsifier": next(
            c["falsifier"] for c in ot.CENSUS if c["module"] == "genomeos/attribution/closure.py"
        ),
        "chrom": chrom,
        "committed_reading_unchanged": {
            "function": "attributed_elements",
            "genes": len(old),
            "memberships": sum(len(v) for v in old.values()),
            "note": "judge() still calls this with no reader, so no committed closure figure moves",
        },
        "window_reading": {
            "function": "window_elements",
            "genes": len(new),
            "memberships": census["memberships"],
        },
        "census": census,
        "genes_that_gain_an_element": len(gained),
        "element_memberships_gained": sum(len(v) for v in gained.values()),
        "genes_new_to_the_closure": len(new_genes),
        "genes_new_to_the_closure_named": new_genes,
        "genes_that_lose_an_element": len(lost),
        "genes_that_lose_an_element_named": sorted(lost),
        "verdict_against_its_own_falsifier": (
            "FIRES on the gain half and does NOT fire on the loss half"
            if gained and not lost
            else "does not fire"
            if not gained and not lost
            else "FIRES on the loss half, which must not happen: the head is always at its own bar"
        ),
        "the_closure_test_can_now_judge_genes_it_could_not": (
            "a gene with no attributed element cannot be judged by this module at all, so the genes "
            "new to the closure are the part of this arm that changes what the test can reach. "
            "Whether their closure HOLDS is not asked here and is not answered by this run"
        ),
        "seconds": round(time.time() - t, 2),
    }


CHROMS = tuple(f"chr{c}" for c in [*range(1, 23), "X", "Y"])


def arm_measured_genome_wide(responses: ElementResponses) -> dict[str, Any]:
    """measured.py over EVERY chromosome, because its falsifier lives in the CRISPRi population.

    chr21 carries 18 measured CRISPRi rows and the falsifier does not fire on them. Eighteen is not
    an answer about this module, so the arm is run genome-wide: the screen is K562 and its pairs are
    spread across the genome. One chromosome's archive is held at a time, so the memory bill is the
    largest chromosome and not the tree.
    """
    t = time.time()
    per: dict[str, Any] = {}
    arms = []
    for ch in CHROMS:
        try:
            coding = coding_of(ch)
        except Exception as ex:  # a chromosome with no GENCODE file on this machine
            per[ch] = {"skipped": f"no GENCODE v50 file: {type(ex).__name__}"}
            continue
        a = arm_measured(ch, responses, coding)
        arms.append(a)
        per[ch] = {
            "rows": a["rows"],
            "crispri_rows": a["rows_with_a_crispri_measurement"],
            "fires": a["falsifier_fires_on"],
            "head_disagrees": a["census"]["head_disagrees"],
            "coding_head_disagrees": a["census"]["coding_head_disagrees"],
            "not_cached": a["census"]["not_cached"],
        }
    fired = [{**f, "chrom": a["chrom"]} for a in arms for f in a["fired"]]
    total = {
        k: sum(p.get(k, 0) for p in per.values() if "skipped" not in p)
        for k in ("rows", "crispri_rows", "fires", "head_disagrees", "coding_head_disagrees", "not_cached")
    }
    ot.check_head_invariant({"head_disagrees": total["head_disagrees"]}, "measured.rows genome-wide")
    return {
        "module": "genomeos/attribution/measured.py",
        "scope": "genome-wide",
        "falsifier": arms[0]["falsifier"] if arms else "",
        "chromosomes_read": sorted(k for k, v in per.items() if "skipped" not in v),
        "chromosomes_skipped": {k: v["skipped"] for k, v in per.items() if "skipped" in v},
        "committed_fields_identical_on_every_chromosome": all(a["committed_fields_identical"] for a in arms),
        "total": total,
        "rows_whose_head_already_agrees": sum(a["rows_whose_head_already_agrees"] for a in arms),
        "rows_where_a_regulated_gene_is_anywhere_in_the_window": sum(
            a["rows_where_a_regulated_gene_is_anywhere_in_the_window"] for a in arms
        ),
        "per_chromosome": per,
        "fired": fired,
        "verdict_against_its_own_falsifier": (
            f"FIRES on {total['fires']} element(s) across {len({f['chrom'] for f in fired})} chromosome(s)"
            if fired
            else "does not fire genome-wide, which would retire this module"
        ),
        "what_the_firing_rows_are_and_are_NOT": arms[0]["what_the_firing_rows_are_and_are_NOT"]
        if arms
        else "",
        "seconds": round(time.time() - t, 2),
    }


def arm_measured(chrom: str, responses: ElementResponses, coding: set[str]) -> dict[str, Any]:
    """measured.py. Falsifier, as registered: an element whose measured regulated gene is at the bar
    in the window but is not the head."""
    from genomeos.attribution import measured as me

    t = time.time()
    plain = me.rows(chrom)
    windowed = me.rows(chrom, responses=responses, coding=coding)
    control = all(
        {k: v for k, v in b.items() if k != "window"} == a for a, b in zip(plain, windowed, strict=True)
    )
    census = {
        "elements": len(windowed),
        "not_cached": sum(1 for r in windowed if r["window"]["not_cached"]),
        "head_disagrees": sum(1 for r in windowed if r["window"]["head_agrees"] is False),
        "coding_head_disagrees": sum(1 for r in windowed if r["window"]["coding_head_agrees"] is False),
    }
    ot.check_head_invariant(census, f"measured.rows {chrom}")
    cris = [r for r in windowed if "crispri" in r["assays"]]
    fired = [
        r
        for r in cris
        if r["agreement"].get("crispri") != me.AGREES and r["window"]["regulated_gene_in_window"]
    ]
    return {
        "module": "genomeos/attribution/measured.py",
        "falsifier": next(
            c["falsifier"] for c in ot.CENSUS if c["module"] == "genomeos/attribution/measured.py"
        ),
        "chrom": chrom,
        "committed_fields_identical": control,
        "rows": len(windowed),
        "rows_with_a_crispri_measurement": len(cris),
        "census": census,
        "rows_whose_head_already_agrees": sum(1 for r in cris if r["agreement"].get("crispri") == me.AGREES),
        "rows_where_a_regulated_gene_is_anywhere_in_the_window": sum(
            1 for r in cris if r["window"]["best_rank_of_a_regulated_gene"]
        ),
        "falsifier_fires_on": len(fired),
        "fired": [
            {
                "id": r["id"],
                "compiled_claim_names": r["predicted_gene"],
                "crispri_verdict_on_that_claim": r["agreement"].get("crispri"),
                "regulated_genes_in_the_same_window": r["window"]["regulated_genes_in_window"],
                "rank_in_the_window": r["window"]["best_rank_of_a_regulated_gene"],
            }
            for r in fired
        ],
        "verdict_against_its_own_falsifier": (
            f"FIRES on {len(fired)} element(s)"
            if fired
            else "does not fire on this population: no element whose measured regulated gene is at "
            "the bar and is not the head. On a population this small that is a statement about the "
            "population, not a retirement of the module"
        ),
        "what_the_firing_rows_are_and_are_NOT": (
            "every verdict among them must be read as it stands. A row fires when the compiled "
            "claim's own gene did not get an `agrees` AND a gene the screen found regulated is at "
            "the bar in the same element's window. If that verdict is "
            "`predicted_gene_not_tested`, the screen never tested the model's gene, so the model "
            "was NOT contradicted at that element -- it was not asked, and the window shows the "
            "screen did regulate a gene the model also named. That is a claim made "
            "unfalsifiable by the projection, which is not the same as a claim shown wrong, and "
            "neither is it evidence that the element regulates that gene"
        ),
        "seconds": round(time.time() - t, 2),
    }


def arm_not_open_profile(chrom: str, responses: ElementResponses, coding: set[str]) -> dict[str, Any]:
    """not_open_profile.py. Falsifier, as registered: the number of predicted rules changing."""
    from genomeos.attribution import not_open_profile as nop

    t = time.time()
    plain = nop.rules(chrom)
    windowed = nop.window_rules(chrom, responses=responses, coding=coding)
    base = [r for r in plain if r.source == nop.SOURCE_PREDICTED]
    kept = [r for r in windowed if r.source == nop.SOURCE_PREDICTED]
    extra = [r for r in windowed if r.source == nop.SOURCE_PREDICTED_WINDOW]
    return {
        "module": "genomeos/attribution/not_open_profile.py",
        "falsifier": next(
            c["falsifier"] for c in ot.CENSUS if c["module"] == "genomeos/attribution/not_open_profile.py"
        ),
        "chrom": chrom,
        "compiled_rules_unchanged": kept == base,
        "rules_is_not_changed_it_is_called": (
            "window_rules() calls rules() and appends; rules() itself is byte-identical to its "
            "committed form, so tests/test_not_open_profile.py's element-for-element pin to "
            "context_evidence.rule_loci holds by construction rather than by re-checking"
        ),
        "compiled_predicted_rules": len(base),
        "window_rules_added": len(extra),
        "distinct_genes_compiled": len({r.gene for r in base}),
        "distinct_genes_with_the_window": len({r.gene for r in kept + extra}),
        "elements_naming_more_than_one_gene": len({r.element for r in extra}),
        "window_rules_by_effect_band": {
            b: sum(1 for r in extra if r.effect_band == b) for b in nop.EFFECT_BANDS
        },
        "verdict_against_its_own_falsifier": (
            "FIRES: the predicted rule count rises"
            if extra
            else "does not fire: no element in "
            "this population names a second coding gene, and the module is withdrawn on its own "
            "evidence"
        ),
        "window_rules_are_in_no_committed_program": (
            "they carry SOURCE_PREDICTED_WINDOW, which is not in SOURCES, so a caller counting "
            "compiled rules is unchanged. Their `cell` is the ELEMENT's compiled context and not "
            "the gene's own strongest track, which the function's docstring states as a limit"
        ),
        "seconds": round(time.time() - t, 2),
    }


ARMS = {
    "closure": arm_closure,
    "measured": arm_measured,
    "not_open_profile": arm_not_open_profile,
}


def outstanding() -> list[dict[str, str]]:
    """The census's identity modules this run does not move, each with the reason it cannot."""
    out = []
    for c in ot.CENSUS:
        if c["class"] != "identity" or c["module"] in {a["module"] for a in []}:
            continue
        if c["status"] == "peer_held":
            out.append(
                {
                    "module": c["module"],
                    "reason": "a live peer lane holds this path on the work board; the commit gate "
                    "refuses a staged path a peer holds and no hold was released",
                    "class": "peer",
                }
            )
        elif c["status"] == "frozen_closure":
            out.append(
                {
                    "module": c["module"],
                    "reason": "inside the paid study's frozen 52-file import closure; "
                    "sender_closure() hashes the working tree, so even an uncommitted edit refuses "
                    "the frozen send. It clears when the paid run resolves or is declined",
                    "class": "closure",
                }
            )
    return out


def payload(
    chrom: str, arms: list[dict[str, Any]], seconds: float, entries: list[dict[str, Any]]
) -> dict[str, Any]:
    moved = [a["module"] for a in arms]
    blocked = outstanding()
    free_not_yet = [m for m in ot.free_to_move() if m not in moved]
    return {
        "result": RESULT,
        "date": date.today().isoformat(),
        "lane": ot.LANE,
        "registered_at": "24adf33 (code), aed8ae9 (the registration file)",
        "question": ot.QUESTION,
        "the_invariant_the_wave_stands_on": ot.THE_INVARIANT,
        "the_invariants_own_limit_is_load_bearing": (
            "the CODING head is not invariant the way the any-gene head is: the 2026-09-27 result "
            "found 10 elements of 4,794 with no `predicted_coding` field although the cache named a "
            "coding gene. Every arm therefore counts coding_head_disagrees BY NAME and only an "
            "ANY-GENE disagreement stops the run"
        ),
        "a_falsifier_per_module": (
            "each arm carries the falsifier its own census entry registered, quoted from the "
            "registration, and is judged against that one alone"
        ),
        "modules_moved": moved,
        "modules_free_but_not_yet_moved": free_not_yet,
        "modules_outstanding": blocked,
        "this_file_does_not_imply_completeness": (
            f"{len(moved)} of the census's {len(ot.by_class('identity'))} identity modules are "
            "moved here. The rest are listed above with the reason each is not"
        ),
        "arms": arms,
        "requests": {
            "alphagenome_requests": 0,
            "money": "none",
            "network": "none: the response cache is already on this machine",
        },
        "seconds": round(seconds, 2),
        "result_manifest": {
            "sources": [
                {
                    "accession": "the AlphaGenome deletion sweep's per-element response cache",
                    "version": f"{ELEMENT_CACHE}, written by scripts/enhancer_targets_all.py at "
                    "threshold=0.0; machine-local and git-ignored by the data boundary",
                },
                {
                    "accession": "GENCODE v50, per chromosome",
                    "version": f"{REFERENCE}/gencode_v50_<chrom>.gff3.gz, the same file the "
                    "2026-09-27 one-target lane read its coding set from",
                },
            ],
            "inputs": entries,
            "input_count": len(entries),
            "assembly": "GRCh38: the assembly the deletion sweep and the compact tables were written on",
            "coordinates": "n/a: this run compares gene memberships and element ids; no interval of "
            "its own is formed",
            "parameters": {
                "chrom": chrom,
                "min_effect": 0.1,
                "bar_rule": "predict_target's own size rule applied to every gene of the window "
                "instead of to its head only (targets.genes_at_bar)",
                "coding_set": "protein_coding symbols of the chromosome's GENCODE v50 file",
            },
            "exclusions": [
                "every module of the census whose class is not identity: two are registered as "
                "invariant and NOT to be moved, five read a head without being a live consumer",
                "every module listed in modules_outstanding, with its reason",
            ],
            "partitions": {
                "committed_reading": "attributed_elements, one coding gene per element, unchanged",
                "window_reading": "window_elements, every coding gene at the bar",
            },
            "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
        },
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", default="chr21")
    ap.add_argument("--modules", default="closure,measured,not_open_profile")
    ap.add_argument(
        "--measured-genome-wide",
        action="store_true",
        help="run the measured.py arm over every chromosome: its falsifier lives in the CRISPRi "
        "population, and chr21 holds 18 of its 1,505 rows",
    )
    ap.add_argument("--no-save", action="store_true")
    args = ap.parse_args()

    ot.check_all()
    t = time.time()
    coding = coding_of(args.chrom)
    responses = ElementResponses()
    arms = []
    for name in args.modules.split(","):
        if not name:
            continue
        if name == "measured" and args.measured_genome_wide:
            arms.append(arm_measured_genome_wide(responses))
        else:
            arms.append(ARMS[name](args.chrom, responses, coding))
    read_chroms = CHROMS if args.measured_genome_wide else (args.chrom,)
    p = payload(args.chrom, arms, time.time() - t, inputs(read_chroms, args.chrom))
    for a in arms:
        print(f"{a['module']}: {a['verdict_against_its_own_falsifier']}")
    if args.no_save:
        print("not saved (--no-save)")
        return 0
    print(f"wrote: {save_result(RESULT, p)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
