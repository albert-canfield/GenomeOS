# SPDX-License-Identifier: AGPL-3.0-or-later
"""The per-assertion direction-v2 emitter: the class comes from the rule, and the rule is unchanged.

Every test on the classifier is on the classifier and on the wording it must carry, and reads no
result, so it passes or fails on the code and not on a number.

Four tests do read trees, and the docstring that said none of them did is amended rather than kept:
the historical premise is pinned to `ad07bf5^` by sha, a second pin shows the probe can see a
carrier, the live claim is narrowed to "every carrier of a v2 reason descends from the v2 result",
and the shapes that rule must tell apart are planted in a temporary directory. A committed result's
absence fails here; it never skips.
"""

from __future__ import annotations

import json
import subprocess
from collections.abc import Callable
from pathlib import Path

import pytest

from genomeos.attribution import direction_v2 as dv
from genomeos.attribution import respmap_v2 as rv
from genomeos.predict import enhancer_target as et
from tests.committed_data import must_be_committed

ROOT = Path(__file__).resolve().parents[1]
SOURCE_RELATIVE = "data/results/response_map_increment2.json"
INCREMENT_3_RELATIVE = "data/results/response_map_increment3.json"
SOURCE = ROOT / SOURCE_RELATIVE


def assertion(
    *,
    element: str = "EH38E0000001",
    chrom: str = "chr21",
    gene: str = "GENE1",
    cell: str = "K562",
    action: str = "activates",
    value: float = -0.8,
) -> dict:
    """One response-map assertion, shaped as `response_map2._compiled_rule` writes it."""
    pc = {
        "gene": gene,
        "action": action,
        "log2_fold_change": value,
        "tissue": cell,
        "basis": "predicted: expression change on deleting the element",
    }
    return {
        "id": f"r2|{chrom}:100-200|{element}|{gene}",
        "status": "predicted",
        "measurement": {"direction": action, "value": value, "is_a_measurement": False},
        "quoted": {"id": element, "predicted_coding": pc},
    }


def row(**fields) -> dict:
    """A cached gene row, with only the fields `direction_v2.retained_values` reads."""
    base = {
        "gene": "GENE1",
        "max_drop_tissue": "",
        "max_drop_log2fc": 0.0,
        "max_rise_tissue": "",
        "max_rise_log2fc": 0.0,
        "by_cell": {},
    }
    base.update(fields)
    return base


# ---- the premise this lane had to check before counting anything ------------------------------


def test_direction_v2_script_writes_no_file() -> None:
    """The premise check, held as a test: the counting script records no per-assertion class."""
    src = (ROOT / "scripts/direction_v2.py").read_text()
    for forbidden in ("json.dump", "write_text", "save_result", "open("):
        assert forbidden not in src, f"{forbidden} now appears; the emitter may be redundant"


# ---- the premise this lane had to check, pinned to the tree it is about ------------------------

#: The commit that built the per-assertion emitter. Its message IS the claim: no committed artefact
#: carried the v2 class of a single published assertion, so the class is emitted and not inferred.
EMITTER_COMMIT = "ad07bf5"
#: The tree the claim is about: the one the emitter was written against, named by sha and not by
#: "now", so the claim cannot be falsified by work done after it was true.
TREE_THE_CLAIM_IS_ABOUT = f"{EMITTER_COMMIT}^"
#: The commit that first committed this lane's own v2 result. Pinned for the opposite reason: it is
#: a tree that HAS carriers, so the probe below is shown to see one rather than asserted to.
FIRST_V2_RESULT_COMMIT = "230509a"
#: Every result whose name starts with this is the registration itself or its own result: the root
#: of the descent relation, not a descendant of it.
V2_RESULT_STEM = "respmap_direction_v2"

PINNED_NOT_LIVE = (
    "amended additively, 2026-10-02. 'No committed result carried a v2 reason token' is a claim "
    "about the tree the emitter was built against, and it is now pinned to that tree by sha. "
    "Asserted against the live tree instead, it was a verdict that moved while nothing it covers "
    "moved: cellcover.json, clause1.json and both of their registrations have since been committed, "
    "every one of them digesting data/results/respmap_direction_v2.json among its manifest inputs, "
    "and the assertion read four registered descendants as a defect and failed a push. The live "
    "claim is narrowed, not dropped: a carrier must have a registered parent. docs/LESSONS.md, "
    "'The verdict moved twice without the mechanism moving once'"
)

ABSENCE_IS_A_DEFECT = (
    "a committed result's absence FAILS and does not skip. Three guards here skipped on it -- the "
    "data/results directory, response_map_increment2.json and response_map_increment3.json, all "
    "three in the commit -- and a stub may substitute a dependency's behaviour, never its "
    "existence. tests/local_data.py's needs_local_data marker is for the git-IGNORED stores and "
    "raises on a path git does not ignore, which is this same distinction asked the other way round"
)


def _git(*args: str) -> str:
    """git's stdout. Exit 1 is a grep's 'no match'; anything above it is a refusal, never an answer.

    Decoded with `errors="replace"` because the descent rule follows a manifest input wherever it
    points and some committed inputs are gzip: measured, a parent named `*.json.gz` made this raise
    UnicodeDecodeError past the JSON guard below, so a planted orphan failed the test with the wrong
    reason. Undecodable bytes are not a JSON object, which is the answer `payload_in_tree` needs.
    """
    r = subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, errors="replace")
    if r.returncode not in (0, 1):
        raise AssertionError(f"git {' '.join(args)} refused ({r.returncode}): {r.stderr.strip()[-300:]}")
    return r.stdout


def carriers_in_tree(revision: str) -> list[str]:
    """The results in `revision`'s tree holding any v2 reason token, by path, read from the tuple.

    All eight tokens, not the one the first version of this test grepped for: the tokens come from
    `dv.UNRESOLVED_REASONS` and are never typed here.
    """
    patterns: list[str] = []
    for token in dv.UNRESOLVED_REASONS:
        patterns += ["-e", token]
    out = _git("grep", "-l", "--fixed-strings", *patterns, revision, "--", "data/results")
    return sorted(line.split(":", 1)[1] for line in out.splitlines() if ":" in line)


def results_in_tree(revision: str) -> frozenset[str]:
    """Every path under data/results that `revision` commits."""
    return frozenset(
        p for p in _git("ls-tree", "-r", "--name-only", revision, "--", "data/results").splitlines() if p
    )


def payload_in_tree(revision: str, path: str) -> dict | None:
    """`path`'s JSON object in `revision`, or None when it is not an object this rule can read."""
    try:
        got = json.loads(_git("show", f"{revision}:{path}"))
    except (json.JSONDecodeError, AssertionError):
        return None
    return got if isinstance(got, dict) else None


def is_the_v2_result(path: str) -> bool:
    """Whether `path` is the v2 result or its registration, which are the root and not descendants."""
    return Path(path).name.startswith(V2_RESULT_STEM)


def manifest_input_paths(payload: dict) -> list[str]:
    """Every path `result_manifest.inputs` names, a group input's members included.

    `inputs` and not the prose in `sources`: an input is digested by sha256, and that is what makes
    naming a parent a registration rather than a mention.
    """
    out: list[str] = []
    for entry in (payload.get("result_manifest") or {}).get("inputs") or []:
        if not isinstance(entry, dict):
            continue
        if isinstance(entry.get("path"), str):
            out.append(entry["path"])
        for member in entry.get("members") or []:
            if isinstance(member, dict) and isinstance(member.get("path"), str):
                out.append(member["path"])
    return out


def carriers_from_nowhere(
    carriers: list[str],
    holds: Callable[[str], bool],
    payload: Callable[[str], dict | None],
) -> list[str]:
    """Those `carriers` that are not registered descendants of the v2 result.

    A descendant digests `respmap_direction_v2.json` or its registration among its manifest inputs,
    directly or through another result that does. `holds` and `payload` are the only two questions
    the rule asks of a tree, so the same rule runs against a git revision and against a plain
    directory -- the second is how the test below is shown to be able to fail at all.
    """

    def descends(name: str, chain: frozenset[str]) -> bool:
        if name in chain:
            return False  # a cycle of results naming each other registers nothing
        got = payload(name)
        if got is None:
            return False
        for parent in manifest_input_paths(got):
            if is_the_v2_result(parent):
                return True
            if holds(parent) and descends(parent, chain | {name}):
                return True
        return False

    return sorted(c for c in carriers if not is_the_v2_result(c) and not descends(c, frozenset()))


# `must_be_committed` was written here at 4f44dbf and is now `tests/committed_data.py`, MOVED rather
# than copied: 55 skip guards across 30 test files stood on paths git tracks, so a second copy of the
# rule would have been the defect the rule is about. The local `ABSENCE_IS_A_DEFECT` above is kept
# because it names the three guards THIS file converted; the shared module carries the general form.


def test_no_committed_result_carried_a_v2_class_when_the_emitter_was_built() -> None:
    """Had any result carried a reason token, the class could have been read instead of emitted.

    The same claim as before, against the tree it was made about. Pinned, it is permanent.
    """
    assert carriers_in_tree(TREE_THE_CLAIM_IS_ABOUT) == []
    assert "pinned to that tree by sha" in PINNED_NOT_LIVE


def test_the_carrier_probe_sees_the_carriers_of_a_tree_that_has_them() -> None:
    """Non-vacuity, pinned too: at 230509a exactly the v2 result and its registration carry a reason.

    Without this, the pin above could pass because the probe found nothing anywhere.
    """
    assert carriers_in_tree(FIRST_V2_RESULT_COMMIT) == [
        "data/results/respmap_direction_v2.json",
        "data/results/respmap_direction_v2_registration.json",
    ]


def test_every_committed_carrier_of_a_v2_reason_descends_from_the_v2_result() -> None:
    """The live claim, narrowed to what stays invariant while results land: a carrier has a parent.

    A new result may carry a v2 reason -- four have -- but only by digesting the v2 result among its
    manifest inputs, directly or through a result that does. One carrying a v2 reason from nowhere
    is the defect the first version of this assertion was reaching for, and it is still caught.
    """
    assert must_be_committed("data/results").is_dir()
    carriers = carriers_in_tree("HEAD")
    assert carriers, "HEAD commits no carrier at all, so the v2 result is no longer committed"
    orphans = carriers_from_nowhere(
        carriers,
        results_in_tree("HEAD").__contains__,
        lambda name: payload_in_tree("HEAD", name),
    )
    assert orphans == [], f"a committed result carries a v2 reason from nowhere: {orphans}"


def test_a_carrier_from_nowhere_is_caught_and_a_registered_descendant_is_not(tmp_path: Path) -> None:
    """The narrowing, shown to be able to fail rather than asserted to be.

    One orphan per reason token, so every one of the eight is searched; then the shapes the rule has
    to tell apart: a direct child, a grandchild, a parent named inside a group input, two results
    naming only each other, and one that mentions the v2 result in prose without digesting it.
    """
    bodies: dict[str, dict] = {f"orphan_{token}.json": {"note": token} for token in dv.UNRESOLVED_REASONS}
    token = dv.UNRESOLVED_REASONS[3]
    assert token == "one_track_seen_twice"  # the token the first version of this test grepped for
    bodies["child.json"] = {
        "note": token,
        "result_manifest": {"inputs": [{"path": f"data/results/{V2_RESULT_STEM}.json", "sha256": "x"}]},
    }
    bodies["grandchild.json"] = {
        "note": token,
        "result_manifest": {"inputs": [{"path": "data/results/child.json", "sha256": "x"}]},
    }
    bodies["group_child.json"] = {
        "note": token,
        "result_manifest": {
            "inputs": [
                {
                    "path": "a group",
                    "group": True,
                    "members": [{"path": f"data/results/{V2_RESULT_STEM}_registration.json"}],
                }
            ]
        },
    }
    bodies["cycle_a.json"] = {
        "note": token,
        "result_manifest": {"inputs": [{"path": "data/results/cycle_b.json", "sha256": "x"}]},
    }
    bodies["cycle_b.json"] = {
        "note": token,
        "result_manifest": {"inputs": [{"path": "data/results/cycle_a.json", "sha256": "x"}]},
    }
    bodies["mentions_only.json"] = {
        "note": token,
        "result_manifest": {"sources": [{"version": f"the committed data/results/{V2_RESULT_STEM}.json"}]},
    }
    bodies["not_a_carrier.json"] = {"note": "no reason token here"}
    for name, body in bodies.items():
        (tmp_path / name).write_text(json.dumps(body))

    def read(path: str) -> dict | None:
        p = tmp_path / Path(path).name
        return json.loads(p.read_text()) if p.exists() else None

    carriers = sorted(
        f"data/results/{p.name}"
        for p in tmp_path.glob("*.json")
        if any(t in p.read_text() for t in dv.UNRESOLVED_REASONS)
    )
    assert "data/results/not_a_carrier.json" not in carriers
    assert len(carriers) == len(bodies) - 1
    orphans = carriers_from_nowhere(carriers, lambda path: (tmp_path / Path(path).name).exists(), read)
    assert orphans == sorted(
        [f"data/results/orphan_{t}.json" for t in dv.UNRESOLVED_REASONS]
        + ["data/results/cycle_a.json", "data/results/cycle_b.json", "data/results/mentions_only.json"]
    )


# ---- the class comes from the rule -------------------------------------------------------------


def test_unanimous_two_values_resolve_and_agree() -> None:
    a = assertion(action="activates", value=-0.8)
    r = row(max_drop_tissue="K562", max_drop_log2fc=-0.8, by_cell={"K562": -0.5})
    out = rv.classify(a, row=r)
    assert out["v2_class"] == "resolved_agrees_with_published"
    assert out["v2_action"] == dv.ACTIVATES
    assert out["v2_unresolved_reason"] is None
    assert out["published_value_is_among_retained"] is True


def test_disagreeing_signs_are_unresolved_with_their_reason() -> None:
    a = assertion(value=-0.8)
    r = row(max_drop_tissue="K562", max_drop_log2fc=-0.8, by_cell={"K562": 0.6})
    out = rv.classify(a, row=r)
    assert out["v2_class"] == "unresolved"
    assert out["v2_unresolved_reason"] == "signs_disagree_in_cell"
    assert out["v2_action"] is None


def test_one_value_only_is_unresolved_and_absent_row_names_its_own_reason() -> None:
    a = assertion(value=-0.8)
    one = rv.classify(a, row=row(max_drop_tissue="K562", max_drop_log2fc=-0.8))
    assert one["v2_class"] == "unresolved"
    assert one["v2_unresolved_reason"] == "one_value_only"
    none = rv.classify(a, row=None)
    assert none["v2_class"] == "unresolved"
    assert none["v2_unresolved_reason"] == "no_row_for_target"
    assert none["cached_row_found"] is False


def test_two_fields_one_track_is_unresolved_under_amendment_1() -> None:
    a = assertion(value=-0.8)
    r = row(max_drop_tissue="K562", max_drop_log2fc=-0.8, by_cell={"K562": -0.8})
    assert rv.classify(a, row=r)["v2_unresolved_reason"] == "one_track_seen_twice"


def test_below_the_imported_floor_is_unresolved_and_the_floor_is_not_restated() -> None:
    small = dv.MAGNITUDE_FLOOR / 10
    a = assertion(value=-small)
    r = row(max_drop_tissue="K562", max_drop_log2fc=-small, by_cell={"K562": -small / 2})
    assert rv.classify(a, row=r)["v2_unresolved_reason"] == "below_magnitude_floor"
    assert dv.MAGNITUDE_FLOOR == et.MIN_EFFECT


def test_absence_gloss_correction_is_counted_and_not_glossed() -> None:
    """The escape branch: the one retained value is by_cell, not the published value, opposite sign."""
    a = assertion(value=-0.8)
    r = row(by_cell={"K562": 0.6})
    out = rv.classify(a, row=r)
    assert out["v2_unresolved_reason"] == "one_value_only"
    assert out["published_value_is_among_retained"] is False
    assert out["only_retained_is_by_cell_opposite_to_published"] is True
    # and the same-sign case is not counted as the branch
    same = rv.classify(a, row=row(by_cell={"K562": -0.6}))
    assert same["only_retained_is_by_cell_opposite_to_published"] is False


def test_a_flip_would_be_recorded_as_a_flip_if_the_rule_ever_produced_one() -> None:
    """The flip class is reachable in code, so a measured 0 is a measurement and not an impossibility.

    This plants a row whose retained values are unanimous and opposite to the published action,
    which the registered reading says the committed cache cannot present; the point is only that
    the classifier would not swallow it.
    """
    a = assertion(action="activates", value=-0.8)
    r = row(max_drop_tissue="K562", max_drop_log2fc=0.9, by_cell={"K562": 0.4})
    out = rv.classify(a, row=r)
    assert out["v2_class"] == "resolved_opposite_to_published"
    assert out["v2_action"] == dv.INHIBITS
    assert out["published_value_is_among_retained"] is False


def test_represses_and_inhibits_are_one_axis() -> None:
    a = assertion(action="represses", value=0.8)
    r = row(max_rise_tissue="K562", max_rise_log2fc=0.8, by_cell={"K562": 0.4})
    out = rv.classify(a, row=r)
    assert out["published_axis"] == "represses_target"
    assert out["v2_action"] == dv.INHIBITS
    assert out["v2_class"] == "resolved_agrees_with_published"


# ---- the shape of the result -------------------------------------------------------------------


def test_classes_are_exhaustive_and_counts_add_to_the_denominator() -> None:
    rows = [
        rv.classify(
            assertion(), row=row(max_drop_tissue="K562", max_drop_log2fc=-0.8, by_cell={"K562": -0.4})
        ),
        rv.classify(
            assertion(), row=row(max_drop_tissue="K562", max_drop_log2fc=-0.8, by_cell={"K562": 0.4})
        ),
        rv.classify(assertion(), row=None),
    ]
    c = rv.counts(rows)
    assert sum(c["by_class"].values()) == c["denominator"] == 3
    assert set(c["by_class"]) == set(rv.CLASSES)
    assert set(c["unresolved_by_reason"]) <= set(dv.UNRESOLVED_REASONS)
    absent = set(c["unresolved_reasons_with_no_assertion"])
    assert absent | set(c["unresolved_by_reason"]) == set(dv.UNRESOLVED_REASONS)


def test_a_reshaped_id_fails_loudly() -> None:
    a = assertion()
    a["id"] = "chr21:100-200"
    with pytest.raises(ValueError):
        rv.locus_of(a)
    b = assertion(element="EH38E0000001")
    b["quoted"]["id"] = "EH38E0000002"
    with pytest.raises(ValueError):
        rv.locus_of(b)


def test_only_predicted_assertions_are_in_the_population() -> None:
    payload = {
        "assertions": [
            assertion(),
            {"id": "x", "status": "observed"},
            {"id": "y", "status": "inferred"},
        ]
    }
    got = rv.predicted_assertions(payload)
    assert [a["status"] for a in got] == ["predicted"]


def test_the_published_payload_fixes_the_denominator() -> None:
    payload = json.loads(must_be_committed(SOURCE_RELATIVE).read_text())
    assert payload["counts"]["by_status"]["predicted"] == 166
    assert payload["counts"]["by_status"]["observed"] == 361
    assert payload["counts"]["assertions"] == 527
    assert len(rv.predicted_assertions(payload)) == 166


def test_increment_3_is_a_different_population() -> None:
    """The brief attached increment 2's counts to increment 3's name; the code records which is which."""
    counts = json.loads(must_be_committed(INCREMENT_3_RELATIVE).read_text())["counts"]
    assert counts["assertions"] == 574
    assert counts["by_status"]["predicted"] == 73
    assert counts["chains"] == 93
    assert "574 assertions over 337 entities in 93 chains" in rv.INCREMENT_3_IS_NOT_THIS_POPULATION


# ---- the readings this lane may not weaken -----------------------------------------------------


def test_the_registered_readings_are_carried_and_not_restated() -> None:
    assert dv.FLIPS_ARE_STRUCTURALLY_IMPOSSIBLE in rv.FLIP_IS_REPORTED_AS_MEASURED
    assert rv.DIAGNOSTIC_IS_NOT_A_V2_OUTPUT == dv.SELECTION_EXCLUDED_DIAGNOSTIC
    assert rv.UNRESOLVED_MEANS == dv.UNRESOLVED_MEANS
    assert dv.COUNTS_VALIDATE_NOTHING in rv.VALIDATES_NOTHING
    assert dv.NOT_VALIDATION in rv.VALIDATES_NOTHING


def test_the_too_strong_absence_gloss_is_not_repeated() -> None:
    text = " ".join(v for v in vars(rv).values() if isinstance(v, str)) + " ".join(
        s for t in vars(rv).values() if isinstance(t, tuple) for s in t if isinstance(s, str)
    )
    assert "may NOT be glossed" in rv.ABSENCE_GLOSS_CORRECTED
    assert "one value retained for that cell" not in text
    assert "388,997" not in text


def test_the_two_overlapping_tallies_are_not_added_anywhere() -> None:
    assert "may not be added" in rv.EARLIER_TALLIES_NOT_TOUCHED
    assert "14 `signs_disagree_in_cell`" in rv.EARLIER_TALLIES_NOT_TOUCHED
    # no module constant performs arithmetic on them
    src = (ROOT / "genomeos/attribution/respmap_v2.py").read_text()
    assert "14 +" not in src and "+ 10" not in src


def test_the_lane_makes_no_recommendation_about_the_grammar_token() -> None:
    assert "makes no recommendation" in rv.NO_RECOMMENDATION
    for word in ("should", "recommend that", "we recommend"):
        assert word not in rv.NO_RECOMMENDATION.replace("makes no recommendation", "")


def test_nothing_is_rescored_and_no_threshold_is_defined_here() -> None:
    src = (ROOT / "genomeos/attribution/respmap_v2.py").read_text()
    assert "MIN_EFFECT =" not in src
    assert "MAGNITUDE_FLOOR =" not in src
    assert "no threshold of its own" in rv.NOTHING_RESCORED


def test_grouping_is_not_independence() -> None:
    joined = " ".join(rv.CANNOT_ESTABLISH)
    assert "NOT established biological independence" in joined


def test_the_reason_count_is_read_and_not_typed() -> None:
    """Amended additively: the brief said nine reason tokens and the tuple holds eight."""
    assert len(dv.UNRESOLVED_REASONS) == 8
    src = (ROOT / "scripts/respmap_v2_register.py").read_text()
    assert "len(dv.UNRESOLVED_REASONS)" in src
    assert "amendment_reason_count" in src
