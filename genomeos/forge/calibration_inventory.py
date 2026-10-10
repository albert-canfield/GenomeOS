# SPDX-License-Identifier: AGPL-3.0-or-later
"""The inventory that could become BioForge's calibration record, and the rule that decides if it is one.

Area H, 2026-10-10 (docs/BIOFORGE-CONFIDENCE.md). A design answer carries no probability because no
calibration record exists: no population of design answers with published outcomes. This module states,
before anything is counted, what such a population would have to be: a set of (design question,
BioForge's answer, published outcome) triples whose outcomes the model did not already hold. It fixes the
smallest population at which a calibration is worth computing, from the width of an exact binomial
interval, and it fixes the rules that mark a triple circular. The census that applies them follows the
registration in this file; it never fabricates an outcome, and it computes no match rate on a population
the registration calls circular, because an answer that reproduces its own training phenotype is not a
prediction.
"""

from __future__ import annotations

import dataclasses
import math
import re
from pathlib import Path
from typing import Any

REGISTERED_ON = "2026-10-10"
CONFIDENCE_LEVEL = 0.95
MAX_INTERVAL_WIDTH = 1 / 3


def _binom_cdf(k: int, n: int, p: float) -> float:
    """P(X <= k) for X ~ Binomial(n, p), summed exactly in floating point (n is small here)."""
    if k < 0:
        return 0.0
    if k >= n:
        return 1.0
    return min(1.0, sum(math.comb(n, i) * p**i * (1.0 - p) ** (n - i) for i in range(k + 1)))


def _bisect(f, lo: float = 0.0, hi: float = 1.0, steps: int = 80) -> float:
    """The root of a function increasing on [lo, hi]."""
    for _ in range(steps):
        mid = (lo + hi) / 2
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def clopper_pearson(k: int, n: int, level: float = CONFIDENCE_LEVEL) -> tuple[float, float]:
    """The exact (Clopper-Pearson) two-sided interval for a binomial rate, k successes of n."""
    if n <= 0 or not 0 <= k <= n:
        raise ValueError(f"need 0 <= k <= n and n > 0, got k={k}, n={n}")
    a = (1.0 - level) / 2
    # lower: the p at which P(X >= k) = a; P(X >= k) rises with p
    lo = 0.0 if k == 0 else _bisect(lambda p: (1.0 - _binom_cdf(k - 1, n, p)) - a)
    # upper: the p at which P(X <= k) = a; P(X <= k) falls with p, so bisect its negation
    hi = 1.0 if k == n else _bisect(lambda p: a - _binom_cdf(k, n, p))
    return lo, hi


def worst_width(n: int, level: float = CONFIDENCE_LEVEL) -> float:
    """The widest exact interval any outcome count can give at n: the width a record must plan for."""
    return max(hi - lo for lo, hi in (clopper_pearson(k, n, level) for k in range(n + 1)))


def minimum_triples(max_width: float = MAX_INTERVAL_WIDTH, level: float = CONFIDENCE_LEVEL) -> int:
    """The smallest n at which every possible count gives an exact interval narrower than `max_width`."""
    n = 1
    while worst_width(n, level) >= max_width:
        n += 1
    return n


# Fixed here and pinned by tests/test_forge_calibration_inventory.py against `minimum_triples()`.
MINIMUM_INDEPENDENT_TRIPLES = 38

REGISTRATION: dict[str, Any] = {
    "written": (
        "2026-10-10, by lane-forgecal, after an exploratory read of which programs hold design and"
        " experiment blocks and whose citations their modules repeat, and BEFORE the census below was"
        " coded, before any inverse design question was posed to BioForge and before any match was"
        " counted. The exploratory read already suggested that no triple would be independent; the"
        " minimum below is a function of the interval width alone and was not chosen with any count in"
        " view, so it could not be tuned to make the census pass or fail"
    ),
    "the_triple": (
        "(design question, BioForge's answer, published outcome): a question BioForge can be asked as"
        " a `design` block over an organism program (targets and keeps in assert grammar, a candidate"
        " list of knockouts or additions), the perturbation BioForge reports as its best feasible"
        " candidate, and the perturbation a cited publication reports for the same phenotype"
    ),
    "the_match": (
        "exact set equality between the knockouts and additions BioForge reports and those the"
        " publication reports. A tie between equally good candidates that BioForge breaks by label"
        " order counts as what BioForge reported; the tie is recorded beside the triple"
    ),
    "the_minimum": {
        "n": MINIMUM_INDEPENDENT_TRIPLES,
        "rule": (
            "the smallest n at which the exact 95% Clopper-Pearson interval on the match rate is"
            " narrower than 1/3 for EVERY possible count k in 0..n (worst case k = n/2)"
        ),
        "why_one_third": (
            "a probability on a design answer is worth stating only if the record behind it can tell"
            " a reader at least whether the answer is more often right than wrong. An interval"
            " narrower than 1/3 cannot contain both 1/3 and 2/3, so whatever the count, the record"
            " rules out at least one of 'mostly wrong' (rate <= 1/3) and 'mostly right' (rate >= 2/3)."
            " Below that n some count leaves both readings open and the record says nothing a reader"
            " can act on. A narrower width (0.25 needs 67, 0.2 needs 104) would be better and is not"
            " required: this is the floor, not the target"
        ),
        "what_it_buys": (
            "calibration-in-the-large only: one rate for one population of design answers. A"
            " confidence that varies across answers would need at least this many independent"
            " triples PER reliability bin, and discrimination needs both outcomes present"
        ),
    },
    "circular_if_any": {
        "same_source": (
            "the outcome's publication (first author's surname and year) is cited by an element of"
            " the model the question runs over (a decision, timer, rule, cell type or other evidence"
            " bearing element, imports included), or by a script the program names in its own text as"
            " the source of its gates (scripts/intestine_registration.py for the crypt module). The"
            " search then traverses a rule written from the answer"
        ),
        "acceptance_check": (
            "the outcome is shipped as an `experiment` block run against the same model by `bio test`"
            " in CI, so the model is held to reproduce it: an answer that matches it reproduces a"
            " phenotype the model was built and tested to give"
        ),
        "authored_frame": (
            "the question's candidate list was written by hand by an author who held the answer, and"
            " contains it (every shipped `design` block)"
        ),
        "evidence_is_the_model": (
            "the design block's own evidence names the model's requirements, not a publication, as the"
            " source of its expected answer"
        ),
    },
    "independent": (
        "a triple none of the four circularity rules marks. A triple that no rule marks is not proven"
        " independent: the rules catch a shared citation and a shipped check, not an author who read"
        " the paper and cited another. The count of independent triples is therefore an upper bound"
    ),
    "eligible_shapes": (
        "an `experiment` with a published citation (evidence experimental) can be posed as an inverse"
        " design question only when BioForge's grammar can express it: knockouts at the zygote over a"
        " program that grows an organism (stages declared). Additions forced at a later time, ablations"
        " given as an environment, and programs with no organism are counted and set aside"
    ),
    "decision": (
        f"if the independent triples number at least {MINIMUM_INDEPENDENT_TRIPLES}, register a"
        " calibration-in-the-large on them (outcome, population and method named, as genomeos.certainty"
        " .Calibration requires) and compute it. Otherwise record the shortfall and compute no rate:"
        " no inverse design is run and no match is counted on circular triples, because a rate on them"
        " measures retrieval of the model's own rules"
    ),
    "never": (
        "an outcome is never written that is not in the repository with its citation; a design answer"
        " that matches the mutant it was derived from is never counted as an independent outcome"
    ),
}


# ---------------------------------------------------------------------------------------------------------
# The census, coded after the registration above was committed (5df76f2). It applies the four rules as
# written; where a rule needs an operational test, the test is named beside it.
# ---------------------------------------------------------------------------------------------------------

_BLOCK = re.compile(r"^[ \t]*(design|experiment)[ \t]+([A-Za-z_][A-Za-z0-9_]*)[ \t]*\{", re.M)
_YEAR = re.compile(r"\b(1[89]\d\d|20\d\d)\b")
_SURNAME = re.compile(r"((?:(?:van|der|de|von|da|di|del|le) )*[A-Z][A-Za-z'\-]+)")
_NAMED_SCRIPT = re.compile(r"scripts/[A-Za-z0-9_./-]+\.py")
_NOT_EVIDENCE = ("experiments", "designs")  # the outcomes and the questions, not the model


def citation_keys(source: str) -> list[tuple[str, str]]:
    """(first author's surname, year) for each `;`-separated citation in an evidence string."""
    out: list[tuple[str, str]] = []
    for segment in source.split(";"):
        year = _YEAR.search(segment)
        surname = _SURNAME.search(segment[: year.start()]) if year else None
        if year and surname and (surname.group(1), year.group(1)) not in out:
            out.append((surname.group(1), year.group(1)))
    return out


def cites(text: str, key: tuple[str, str]) -> bool:
    """The surname followed, on the same line and within 80 characters, by the year."""
    surname, year = key
    return re.search(rf"\b{re.escape(surname)}\b[^\n]{{0,80}}?\b{year}\b", text) is not None


def model_evidence(module: Any) -> list[tuple[str, str, str]]:
    """(field, element, evidence source) for every evidence-bearing element of the model the program runs:
    decisions, timers, rules, cell types and the rest, imports included, experiments and designs excluded."""
    out = []
    for f in dataclasses.fields(module):
        if f.name in _NOT_EVIDENCE:
            continue
        value = getattr(module, f.name)
        items = value.values() if isinstance(value, dict) else value if isinstance(value, list) else [value]
        for item in items:
            ev = getattr(item, "evidence", None)
            if ev is not None and getattr(ev, "source", ""):
                name = getattr(item, "name", "") or getattr(item, "id", "") or f.name
                out.append((f.name, str(name), ev.source))
    return out


def named_scripts(root: Path, program: str) -> list[str]:
    """Scripts the program's own text names (its header names the registration that placed its gates)."""
    text = (root / program).read_text(encoding="utf-8", errors="replace")
    return sorted({s for s in _NAMED_SCRIPT.findall(text) if (root / s).is_file()})


def same_source_hits(
    keys: list[tuple[str, str]], evidence: list[tuple[str, str, str]], scripts: dict[str, str]
):
    """Where the model, or a gate registration the program names, cites one of the outcome's publications."""
    hits = []
    for key in keys:
        hits += [
            f"{field}:{name} cites {key[0]} {key[1]}" for field, name, src in evidence if cites(src, key)
        ]
        hits += [f"{path} cites {key[0]} {key[1]}" for path, text in scripts.items() if cites(text, key)]
    return hits


def programs_with_blocks(root: Path) -> dict[str, dict[str, int]]:
    """Every `.bio` under the CI `bio test` roots that holds a design or experiment block, with counts."""
    from genomeos.forge.design_answers import PROGRAM_ROOTS

    out: dict[str, dict[str, int]] = {}
    for top in PROGRAM_ROOTS:
        for path in sorted((root / top).rglob("*.bio")):
            kinds = [m.group(1) for m in _BLOCK.finditer(path.read_text(encoding="utf-8", errors="replace"))]
            if kinds:
                out[path.relative_to(root).as_posix()] = {k: kinds.count(k) for k in ("design", "experiment")}
    return out


def _model_key(module: Any) -> tuple[str, ...]:
    """Programs that import the same model files run over the same model."""
    return tuple(sorted(module.imports))


def _bio_test_runs(module: Any) -> bool:
    """`bio test` runs a program's experiments when it grows an organism or is a located program."""
    return bool(getattr(module, "located", False)) or bool(module.stages and module.organism is not None)


def experiment_record(program: str, module: Any, ex: Any, evidence: list, scripts: dict[str, str]) -> dict:
    keys = citation_keys(ex.evidence.source)
    published = ex.evidence.kind.value == "experimental" and bool(keys)
    set_aside = []
    if not ex.knockouts:
        set_aside.append("no knockout")
    if ex.adds or ex.add_at:
        set_aside.append("adds a factor" + (" at a later time" if ex.add_at else ""))
    if ex.environment:
        set_aside.append(
            "changes the environment (" + ", ".join(f"{k}={v}" for k, v in ex.environment.items()) + ")"
        )
    if module.organism is None or not module.stages:
        set_aside.append("the program grows no organism")
    hits = same_source_hits(keys, evidence, scripts)
    reasons = []
    if hits:
        reasons.append("same_source")
    if _bio_test_runs(module):
        reasons.append("acceptance_check")
    return {
        "program": program,
        "experiment": ex.name,
        "perturbation": {"knockouts": list(ex.knockouts), "adds": list(ex.adds)},
        "citation": ex.evidence.source,
        "citation_keys": [f"{a} {y}" for a, y in keys],
        "published": published,
        "expressible_as_a_design_question": published and not set_aside,
        "set_aside_because": set_aside,
        "circular": reasons,
        "same_source_where": hits[:6],
        "independent": published and not set_aside and not reasons,
    }


def design_record(
    program: str, module: Any, result: Any, evidence: list, scripts: dict, outcomes: list
) -> dict:
    """One shipped design block as a triple: the question, BioForge's answer, and the published outcome the
    repository records for it, found either through the design's own citation or as an experiment block over
    the same model whose knockouts are the answer and whose asserts include every one of the design's targets.
    """
    from genomeos.forge.design_answers import answer_record

    d = result.design
    rec = answer_record(program, result)
    answer = sorted(result.best.knockouts) if result.best is not None and result.best.feasible else None
    keys = citation_keys(d.evidence.source)
    by_citation = [o for o in outcomes if set(o["citation_keys"]) & {f"{a} {y}" for a, y in keys}]
    by_phenotype = [
        o
        for o in outcomes
        if sorted(o["perturbation"]["knockouts"]) == answer and set(d.targets) <= set(o["asserts"])
    ]
    found = by_citation or by_phenotype
    outcome = found[0] if found else None
    # the rules the answer acts through: decisions whose condition requires a factor the answer removes
    through = [
        (dec.id, dec.evidence.source)
        for dec in module.decisions
        if answer and any(dec.when.get(f) == "present" for f in answer) and dec.evidence.source
    ]
    rule_keys = [k for _, src in through for k in citation_keys(src)]
    reasons = []
    if same_source_hits(keys, evidence, scripts) or (
        outcome and outcome["circular"] and "same_source" in outcome["circular"]
    ):
        reasons.append("same_source")
    elif outcome is None and not keys and rule_keys:
        reasons.append("same_source")  # the only published record of the answer is the rule it acts through
    if outcome is not None and "acceptance_check" in outcome["circular"]:
        reasons.append("acceptance_check")
    # operational test of a hand-written frame: the candidate list is not the blind rule "every factor the
    # organism declares", and it contains the recorded outcome (or, with none recorded, the answer)
    declared = set(getattr(module.organism, "factors", None) or [])
    wanted = sorted(outcome["perturbation"]["knockouts"]) if outcome is not None else answer
    if set(d.knockout_any_of) != declared and wanted and set(wanted) <= set(d.knockout_any_of):
        reasons.append("authored_frame")
    if not keys:
        reasons.append("evidence_is_the_model")  # operational test: the design cites no publication
    return {
        **rec,
        "question": {
            "targets": list(d.targets),
            "keeps": list(d.keeps),
            "candidates": list(d.knockout_any_of),
        },
        "design_evidence": d.evidence.source,
        "published_outcome": (
            {
                "experiment": None,
                "knockouts": None,
                "citation": d.evidence.source,
                "found_by": "the design's own citation; no experiment block records its perturbation",
            }
            if keys
            else None
        )
        if outcome is None
        else {
            "experiment": f"{outcome['program']} {outcome['experiment']}",
            "knockouts": outcome["perturbation"]["knockouts"],
            "citation": outcome["citation"],
            "found_by": "the design's own citation"
            if by_citation
            else "an experiment over the same model "
            "whose knockouts are the answer and whose asserts include every target",
        },
        "outcome_recorded_only_as_rule_evidence": None
        if outcome is not None
        else [f"{i}: {src}" for i, src in through if citation_keys(src)],
        "answer_matches_the_recorded_outcome": None
        if outcome is None
        else sorted(outcome["perturbation"]["knockouts"]) == answer,
        "circular": reasons,
        "independent": not reasons,
    }


def build_inventory(root: Path | None = None) -> dict[str, Any]:
    """The census: every design and experiment block under the CI roots, as triples, with the four rules
    applied, the independent count, and the registered decision. Runs each shipped design once (cached)."""
    from genomeos.forge.design_answers import repo_root, results_for
    from genomeos.lang import parse_file

    root = repo_root() if root is None else root
    programs = programs_with_blocks(root)
    parsed = {p: parse_file(str(root / p)) for p in programs}
    scripts = {
        p: {s: (root / s).read_text(encoding="utf-8") for s in named_scripts(root, p)} for p in programs
    }
    evidence = {p: model_evidence(m) for p, m in parsed.items()}
    experiments = []
    for p, m in parsed.items():
        for ex in m.experiments:
            rec = experiment_record(p, m, ex, evidence[p], scripts[p])
            rec["asserts"] = list(ex.asserts)
            rec["_model"] = _model_key(m)
            experiments.append(rec)
    designs = []
    for p, m in parsed.items():
        if not m.designs:
            continue
        same_model = [o for o in experiments if o["_model"] == _model_key(m)]
        designs += [design_record(p, m, r, evidence[p], scripts[p], same_model) for r in results_for(p)]
    for rec in experiments:
        rec.pop("_model")
    pub = [e for e in experiments if e["published"]]
    expressible = [e for e in pub if e["expressible_as_a_design_question"]]
    independent = sum(1 for e in expressible if e["independent"]) + sum(
        1 for d in designs if d["independent"]
    )

    def by_rule(recs: list[dict]) -> dict[str, int]:
        return {
            rule: sum(1 for r in recs if rule in r["circular"]) for rule in REGISTRATION["circular_if_any"]
        }

    shortfall = max(0, MINIMUM_INDEPENDENT_TRIPLES - independent)
    return {
        "what": (
            "the inventory that could become a calibration record for BioForge's design answers: every"
            " (design question, BioForge's answer, published outcome) triple the repository holds or could"
            " pose, with the registered circularity rules applied (genomeos/forge/calibration_inventory.py)"
        ),
        "registered": {
            "commit": "5df76f2",
            "on": REGISTERED_ON,
            "minimum_independent_triples": MINIMUM_INDEPENDENT_TRIPLES,
        },
        "programs_scanned": programs,
        "counts": {
            "design_blocks": len(designs),
            "design_triples_with_a_published_outcome_in_the_repository": sum(
                1 for d in designs if d["published_outcome"] is not None
            ),
            "design_triples_whose_outcome_is_only_a_rule_s_evidence": sum(
                1
                for d in designs
                if d["published_outcome"] is None and d["outcome_recorded_only_as_rule_evidence"]
            ),
            "design_triples_circular_by_rule": by_rule(designs),
            "design_triples_independent": sum(1 for d in designs if d["independent"]),
            "design_answers_tied": sum(1 for d in designs if d["tied_with_the_answer"]),
            "experiment_blocks": len(experiments),
            "experiments_with_a_published_citation": len(pub),
            "expressible_as_inverse_design_questions": len(expressible),
            "set_aside_published_experiments": len(pub) - len(expressible),
            "inverse_triples_circular_by_rule": by_rule(expressible),
            "inverse_triples_independent": sum(1 for e in expressible if e["independent"]),
            "independent_triples": independent,
            "minimum_independent_triples": MINIMUM_INDEPENDENT_TRIPLES,
            "shortfall": shortfall,
        },
        "decision": (
            "enough independent triples: a calibration is to be registered before it is computed, and"
            " this census computes none"
            if shortfall == 0
            else f"shortfall of {shortfall}: no calibration is registered and no match rate is computed;"
            " no inverse design question was posed to BioForge"
        ),
        "what_a_calibration_would_need": (
            f"at least {MINIMUM_INDEPENDENT_TRIPLES} design questions whose published outcome (a) is"
            " cited by no element of the model and by no gate registration the program names, (b) is not"
            " shipped as an experiment block or `# test:` line that bio test runs against the same model,"
            " and (c) is posed over a candidate list fixed by a rule written before the question (every"
            " factor the organism declares), with (d) a publication as the design's evidence. In this"
            " repository every published outcome is also an acceptance check of the model it would be"
            " scored on, so such questions can only come from publications held out before a module is"
            " written; none is held out today"
        ),
        "design_triples": designs,
        "experiment_triples": experiments,
        "not_a_claim": (
            "a census of what could be checked, not a check: no answer here is scored as right or wrong for"
            " calibration, and the one field comparing an answer with a recorded outcome"
            " (answer_matches_the_recorded_outcome) reports a circular match, which reproduces the model's"
            " own acceptance phenotype and is not a prediction"
        ),
    }


#: keys the run and save_result add beside the census (manifest, wall time, spend, date, name); not the census
RUN_ONLY_KEYS = ("result_manifest", "seconds", "money", "date", "result")


def census_of(payload: dict[str, Any]) -> dict[str, Any]:
    """A written inventory without what the run added, as `--check` and the tests compare it."""
    return {k: v for k, v in payload.items() if k not in RUN_ONLY_KEYS}
