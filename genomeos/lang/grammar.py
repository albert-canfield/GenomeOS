# SPDX-License-Identifier: Apache-2.0
"""The BioLang grammar and the BioIR type list, generated from the parser and the IR.

`bio grammar` (or `python -m genomeos.lang.grammar`) prints docs/BIOLANG-GRAMMAR.md so that
the written grammar can never drift from what the parser accepts: the block kinds come from
the parser's table, the properties per block from the table below (which the parser's tests
hold to), and the IR types from the dataclasses themselves.
"""

from __future__ import annotations

import dataclasses
import sys

from genomeos.ir import model as ir
from genomeos.lang import parser as _parser

VERSION = "0.4"

# property -> (form, meaning), per block kind; header form per kind
BLOCKS: dict[str, dict] = {
    "gene": {
        "header": "gene <Id> { ... }",
        "props": {
            "symbol": ("text", "gene symbol (defaults to the id)"),
            "locus": ("chrN:start-end[(+|-)]", "genomic locus"),
            "max": ("number", "maximal transcription rate for the network engine"),
            "basal": ("number", "basal transcription rate"),
            "produces": ("Id, Id", "proteins this gene produces (one `produces` rule each)"),
            "location": ("Id", "v0.4: the compartment where it is read"),
            "cost": ("POOL number per UNIT", "repeatable; v0.4: a draw on a declared pool"),
        },
        "nested": "transcript",
    },
    "transcript": {
        "header": "transcript <Id> { ... }   (inside a gene)",
        "props": {"exons": ("locus, locus", "exon loci"), "cds": ("locus, locus", "coding segments")},
    },
    "protein": {
        "header": "protein <Id> { ... }",
        "props": {
            "location": ("Id, Id", "v0.4: compartments it occupies when it works"),
            "signals": ("name, name", "v0.4: targeting signals a transport recognises"),
            "initial": ("number", "v0.4: amount in each declared location at time 0"),
            "half_life": ("hours", "protein half-life"),
            "sequence": ("residues", "amino-acid sequence"),
            "accession": ("text", "UniProt accession"),
            "isoforms": ("Id, Id", "isoform accessions"),
            "domains": ("name, name", "domain names"),
            "pathways": ("Id, Id", "Reactome pathway ids"),
            "interactions": ("Id, Id", "interaction partners"),
            "structures": ("Id, Id", "PDB ids or AF- AlphaFold ids"),
            "cost": ("POOL number per UNIT", "repeatable; v0.4: a draw on a declared pool"),
        },
    },
    "region": {
        "header": "region <Id> { ... }",
        "props": {
            "locus": ("locus", ""),
            "role": (
                "text | unknown",
                "the budget's tier summary, kept verbatim (Module.unknowns counts `unknown`); the five "
                "axes below are the authoritative reading (R7)",
            ),
        },
    },
    "element": {
        "header": "element <Id> { ... }",
        "props": {
            "class": (
                "promoter | enhancer | insulator | open_chromatin | unknown",
                "a summary derived from the registry role (R7), never from the activity; the five axes "
                "below are the authoritative reading",
            ),
            "locus": ("locus", ""),
            "domain": ("Id", "the node (domain) it lies in"),
            "targets": ("Id, Id", "genes it reaches"),
            "basis": ("text", "how the targets were assigned"),
        },
    },
    "domain": {
        "header": "domain <Id> { ... }",
        "props": {
            "locus": ("locus", ""),
            "genes": ("Id, Id", ""),
            "boundaries": ("Id, Id", "insulator elements"),
        },
    },
    "rule": {
        "header": "rule <Source> activates|inhibits|modifies|produces|binds|degrades <Target> { ... }",
        "props": {
            "strength": ("0..1", ""),
            "threshold": (
                "number [M|mM|uM|nM|pM|fM]",
                "source level giving half-maximal effect: an amount, or v0.4 a concentration with a "
                "molar unit, divided by the absolute_volume of the compartment where the rule acts",
            ),
            "hill": ("number", "cooperativity"),
            "when": ("k = v, k = v", "context in which the rule applies"),
            "id": ("text", "explicit rule id"),
        },
    },
    "param": {"header": "param <name> = <number>|unknown [unit] { ... }", "props": {}},
    "cell_type": {
        "header": "cell_type <Id> { ... }",
        "props": {
            "name": ("text", ""),
            "parent": ("Id", "declared parent type"),
            "expresses": ("Id, Id", "genes expressed; the rest are silenced in this context"),
            "ontology": ("CL:0000000", "Cell Ontology id"),
        },
    },
    "event": {
        "header": "event <Id> { ... }",
        "props": {
            "rate": ("number /unit", "rate with unit"),
            "when": ("k = v, ...", "guards"),
            "effect": ("var op number [unit]", "repeatable; op in += -= *= ="),
            "cost": ("POOL number per UNIT", "repeatable; v0.4: a draw on a declared pool"),
            "partition": (
                "duplicate | contents | binomial",
                "v0.4 stage 4: what a division does to the contents; default duplicate "
                "(contents = binomial below the regime threshold, exact halves above)",
            ),
        },
    },
    "organism": {
        "header": "organism <Id> { ... }   (one per program)",
        "props": {
            "species": ("text", ""),
            "genome": ("text", "assembly"),
            "root": ("Id", "name of the first cell (default Zygote)"),
            "resolution": ("cells | populations", "named cells or counted populations"),
            "seed": ("integer", "timers run with their spread from this seed; omit for means"),
            "cell_type": ("Id", "bootstrap cell type"),
            "factors": ("Id, Id", "maternal factors in the first cell"),
            "environment": ("k = v, ...", ""),
            "tempo": ("number", "multiplies every timer"),
            "space": ("W x H", "grid; cells then hold sites"),
            "origin": ("x,y", "site of the first cell"),
            "sense": ("time", "interval between gradient readings"),
            "observe": ("count, deaths, fates, lineage L", "repeatable"),
            "assert": (
                "count|deaths|fate T|type T|lineage L at N unit (in a..b | = n | >= n | <= n); "
                "with replicates: "
                "exactly one of A, B is T [at N unit] in >= P% | A is T [at N unit] in lo..hi%",
                "repeatable",
            ),
            "reference": ("name", "ground truth to diff against"),
            "contacts": ("path", "v0.4: time-resolved contact table (time, cell, cell[, area]); neighbours"),
            "cell_network": ("time", "v0.4: step every cell's own network together at this interval"),
            "replicates": ("integer", "v0.4: runs over seeds 0..N-1 for replicate asserts"),
            "placement": (
                "nearest | names",
                "v0.4: second daughter at the nearest site or along its name's axis",
            ),
        },
    },
    "stage": {"header": "stage <Id> { ... }", "props": {"from": ("time", ""), "to": ("time", "")}},
    "timer": {
        "header": "timer <Id> { ... }",
        "props": {
            "duration": ("time", "required"),
            "sd": ("number", "spread, same unit"),
            "lengthening": ("number", "multiplier per generation after the first"),
            "when": ("k = v, ...", "which cells it times"),
        },
    },
    "field": {
        "header": "field <Id> { ... }",
        "props": {
            "diffusion": ("grid^2 per hour", ""),
            "decay": ("per hour", ""),
            "source": ("x,y = rate", "repeatable; amount per hour"),
        },
    },
    "signal": {
        "header": "signal <Id> { ... }",
        "props": {
            "mode": ("contact | gradient | systemic", ""),
            "ligand": ("Id", ""),
            "receptor": ("Id", ""),
            "from": ("k = v, ...", "sender condition (contact)"),
            "to": ("k = v, ...", "receiver condition"),
            "sets": ("FACTOR [= value]", "factor the receiver gains"),
            "field": ("Id", "gradient: field read at the cell's site"),
            "threshold": ("number", "gradient: factor set at or above this"),
            "reads": (
                "presence | amount",
                "v0.4 contact: a named sender, or the ligand summed over neighbours",
            ),
        },
    },
    "decision": {
        "header": "decision <Id> { ... }",
        "props": {
            "action": ("divide | differentiate | migrate | quiesce | die | express", "required"),
            "when": ("k = v, ...", "`any`, `absent`, `a|b`, `>=n` `<=n` `>n` `<n`"),
            "daughters": ("Id, Id", "divide: names; omitted = a/p then l/r suffixes"),
            "lineages": ("Id = L, ...", "divide: daughters that found a lineage"),
            "asymmetric": ("F -> Id, ...", "divide: which daughter keeps which factor"),
            "direction": ("+x | -x | +y | -y", "divide: where the second daughter goes; migrate: step"),
            "to": ("Id", "differentiate: cell type"),
            "name": ("text", "differentiate: terminal name"),
            "sets": ("F, F", "express: factors the cell carries"),
            "toward": ("Id", "migrate: climb this field"),
            "steps": ("integer", "migrate: sites per move"),
            "timer": ("Id", "divide: explicit timer; omitted = first matching"),
            "after": ("time", "delay from birth (cells) or from now (populations); recurring for flows"),
            "fraction": ("number", "populations: a share of what is LEFT when it runs (divide may exceed 1)"),
            "share": ("number", "populations: a share of the whole at this decision point, order-free"),
            "priority": (
                "integer",
                "highest wins among matching decisions of one action; ties: first in module order",
            ),
            "competence": ("Id", "differentiate: the window this fate change needs open"),
        },
    },
    "experiment": {
        "header": "experiment <Id> { ... }",
        "props": {
            "knockout": ("F, F", "factors never present; signals by id, ligand or receptor never sent"),
            "add": ("F[, F] [at T]", "factors added to the first cell, or forced at a stated time"),
            "environment": ("k = v, ...", "overrides"),
            "until": ("time", "horizon; omitted = the last stage"),
            "expect": ('"text"', "the published phenotype"),
            "assert": ("assert grammar", "repeatable; checked on the mutant"),
        },
    },
    "order": {
        "header": "order <Id> { ... }",
        "props": {
            "members": ("A, B|C, D", "required; the steps in order, `|` for members a source does not order"),
            "axis": ("position | time", "position is checked when it compiles, time against a run"),
            "direction": ("increasing | decreasing", "position only: along the chromosome"),
            "observe": ("birth", "time only: what counts as a step having happened"),
            "threshold": ("number", "time: the level at which a member counts as on"),
        },
    },
    "competence": {
        "header": "competence <Id> { ... }",
        "props": {
            "allows": ("T, T", "required; the fates the window permits"),
            "closes": ("at T | after generation N | on commitment", "required; a window always closes"),
            "closed_by": ("F", "factor whose presence is needed to close it"),
            "when": ("k = v, ...", "the cells the window governs; omitted = every cell"),
        },
    },
    "commitment": {
        "header": "commitment <Id> { ... }",
        "props": {
            "establish": ("k = v, ...", "required; when the cell is committed"),
            "programme": ("T", "what it commits to; omitted = the cell_type it holds then"),
            "locks": ("cell_type", "what can no longer change"),
            "inherit": ("daughters | no", "the lock passes to the daughters"),
            "release": ("never", "only never is implemented"),
            "maintain": ("-", "rejected: specified in BIOLANG-v0.4-ECONOMY.md §7.2a, not implemented"),
            "excludes": ("-", "rejected: specified in BIOLANG-v0.4-ECONOMY.md §7.2a, not implemented"),
            "hysteresis": ("-", "rejected: specified in BIOLANG-v0.4-ECONOMY.md §7.2a, not implemented"),
        },
    },
    "design": {
        "header": "design <Id> { ... }",
        "props": {
            "knockout_any_of": ("F, F", "BioForge may knock these out"),
            "add_any_of": ("F, F", "BioForge may add these"),
            "at_most": ("integer", "perturbations combined per candidate"),
            "vary": ("timer NAME duration LO..HI | decision ID fraction|after LO..HI", "repeatable knobs"),
            "until": ("time", ""),
            "target": ("assert grammar", "repeatable; loss = distance from holding"),
            "keep": ("assert grammar", "repeatable; must hold"),
        },
    },
    "compartment": {
        "header": "compartment <Id> { ... }",
        "props": {
            "parent": ("Id", "the compartment that contains this one; one root"),
            "membrane": ("yes | no", "a membrane faces its parent and its children"),
            "volume": ("number", "fraction of the cell's volume"),
            "absolute_volume": (
                "number fL|pL|nL|uL|mL|L|um3 | unknown",
                "v0.4: this compartment's own volume (1 um3 = 1 fL); what a concentration divides by",
            ),
            "genome": ("chrM, ... | nuclear", "chromosomes read here"),
            "translation": ("yes | no", "ribosomes are present"),
            "copies": ("integer", "copies per cell"),
        },
    },
    "transport": {
        "header": "transport <Id> { ... }",
        "props": {
            "from": ("Id", "required; adjacent to `to`, or across one membrane"),
            "to": ("Id", "required"),
            "cargo": ("Id, Id | mRNA | signal = S", "required; what it carries"),
            "capacity": ("amount per hour", "required; maximal total flux"),
            "affinity": ("amount", "cargo level giving half-maximal flux"),
            "via": ("Id", "protein whose presence gates the capacity"),
            "via_threshold": ("number", "via level giving half capacity"),
        },
    },
    "pool": {
        "header": "pool <Id> { ... }",
        "props": {
            "location": ("Id", "the compartment the pool is in"),
            "size": ("number | unknown", "capacity in molecules; unknown refuses a burden claim"),
            "regenerates": (
                "number /h | from Id, Id",
                "a rate, or the processes that refill it when a rate is the wrong shape (ATP)",
            ),
            "returns_as": ("Id", "what a drawn unit becomes when the work ends (ATP returns as ADP)"),
        },
    },
    "allocation": {
        "header": "allocation <Id> { ... }",
        "props": {
            "pool": ("Id", "the pool this policy divides; required"),
            "policy": (
                "proportional | priority | competitive | optimise",
                "the scientific claim, named in the program and never hidden in the engine",
            ),
            "order": ("Id, Id", "priority: the order demands are served in"),
            "objective": ("text", "optimise: the declared objective; always a modelling device"),
        },
    },
    "regime": {
        "header": "regime <Id> { ... }   (one per program)",
        "props": {
            "treatment": ("continuous | stochastic | auto", "per species; auto uses the threshold"),
            "threshold": ("copies", "auto: stochastic below this"),
            "units": ("au | copies", "stochastic treatment needs copies"),
            "update": ("continuous | synchronous | asynchronous | event", "recorded in every result"),
            "fates": (
                "first | last",
                "Body: one fate per decision point by precedence (default), or the last match (legacy)",
            ),
            "allocation": ("competitive | proportional | priority | optimise", "shared capacities"),
            "seed": ("integer", ""),
            "recheck": ("crossings | none", "a cell decides again when a read it names crosses a threshold"),
        },
    },
}
#: Review R7 (2026-09-28): five axes an `element` or `region` may state, each on its own key, so that
#: where a sequence came from, what it is biochemically, what it was seen or predicted to do, how it
#: relates to a gene and what the claim rests on are never folded into one label. A value may carry a
#: qualifier after `/` (`repeat_derived/LINE`). In a property, `,` joins values that all hold (roles
#: overlap: a CTCF-bound enhancer-like element is both) and `|` joins alternatives of which one holds
#: and none is chosen (`silencer|insulator_like|competing_promoter|unknown`). `unknown` is a value on
#: every axis. No axis has a value meaning "no function": constraint is evidence of selection, and its
#: absence proves nothing, so it lives on `evidence_status`, never on `molecular_role`.
AXES: dict[str, dict[str, str]] = {
    "origin": {
        "unique": "no interspersed repeat over the interval",
        "repeat_derived": "at least half the interval is interspersed repeat (RepeatMasker); /CLASS",
        "partly_repeat_derived": "some but under half of the interval is interspersed repeat; /CLASS",
        "satellite": "satellite or centromeric array",
        "tandem_repeat": "simple or tandem repeat array",
        "segmental_duplication": "at least half the interval is a curated segmental duplication",
        "assembly_gap": "no sequence to attribute",
        "unknown": "origin not read",
    },
    "molecular_role": {
        "promoter_like": "promoter-like biochemical signature (ENCODE PLS, or a CpG island)",
        "enhancer_like": "enhancer-like biochemical signature (ENCODE pELS or dELS)",
        "insulator_like": "CTCF signature (CTCF-only, or a CTCF-bound cCRE)",
        "open_chromatin": "open chromatin with a promoter mark (ENCODE DNase-H3K4me3)",
        "silencer": "a repressive element; never inferred from a direction of effect alone",
        "competing_promoter": "represses a gene by competing for its regulators",
        "structural": "a mechanical role (centromere, satellite array)",
        "coding_candidate": "a long open reading frame that may code",
        "unknown": "no role read",
    },
    "activity": {
        "activates_target": "removing it lowers its target's expression (predicted or measured)",
        "represses_target": "removing it raises its target's expression (predicted or measured)",
        "no_effect_measured": "a well-powered perturbation measured no effect on the genes tested",
        "active_in_reporter": "a reporter assay read it active in at least one cell or tissue",
        "inactive_in_reporter": "a reporter assay read it inactive everywhere it was tested",
        "unknown": "no activity read",
    },
    "target_relation": {
        "predicted_deletion_target": "a model deletion names the gene",
        "nearest_tss_in_domain": "the gene is the nearest TSS in its domain",
        "measured_perturbation_target": "silencing it in place changed the gene (training split)",
        "tested_no_effect": "genes were tested in place and none changed",
        "unassigned": "no gene assigned",
    },
    "evidence_status": {
        "curated_annotation": "an annotation (RepeatMasker, the assembly) states the origin",
        "registry_biochemical": "the ENCODE cCRE registry states the biochemical signature",
        "predicted_model": "a model prediction states the activity or target",
        "measured": "an assay measured it",
        "measured_negative": "an assay measured the absence of an effect",
        "conflicting": "a measurement disagrees with the prediction",
        "under_selection": "constrained across mammals (>= 5% of bases at phyloP >= 2.27)",
        "selection_weak": "some constraint (3% to 5% of bases)",
        "selection_not_detected": (
            "under 3% of bases constrained: no sign of selection, not a sign of no function"
        ),
        "selection_not_measured": "constraint not read over this interval",
        "unknown": "nothing stated",
    },
}
#: the five axes as `element` and `region` properties, after each block's own keys
AXIS_PROPS: dict[str, tuple[str, str]] = {
    axis: (" | ".join(vals), "R7 axis; `,` values all hold, `|` unresolved alternatives, `/` a qualifier")
    for axis, vals in AXES.items()
}
for _kind in ("region", "element"):
    BLOCKS[_kind]["props"].update(AXIS_PROPS)

COMMON = {
    "evidence": ('kind "source" [note]', "kind in experimental, curated, predicted, inferred, none"),
    "confidence": ("0..1", ""),
}
DIRECTIVES = {
    "module": "module <dotted.name>",
    "import": "import <dotted.name> | <relative/path.bio> | bio.std.<name>",
}
IR_TYPES = [
    "Evidence",
    "Entity",
    "Compartment",
    "Transport",
    "Regime",
    "Region",
    "RegulatoryElement",
    "Domain",
    "Signal",
    "Gene",
    "Transcript",
    "Protein",
    "CellType",
    "Effect",
    "Event",
    "Parameter",
    "Rule",
    "Field",
    "Timer",
    "Stage",
    "Decision",
    "Experiment",
    "Design",
    "Organism",
    "Module",
]


def check() -> list[str]:
    """Differences between this table and the parser: block kinds missing on either side."""
    problems = []
    for kind in _parser._KINDS:
        if kind not in BLOCKS:
            problems.append(f"parser accepts {kind!r} but the grammar table does not describe it")
    for kind in BLOCKS:
        if kind not in _parser._KINDS:
            problems.append(f"grammar table describes {kind!r} but the parser does not accept it")
    for name in IR_TYPES:
        if not hasattr(ir, name):
            problems.append(f"IR type {name!r} is listed but not exported")
    return problems


def render() -> str:
    out = [
        f"# BioLang v{VERSION} grammar and BioIR {VERSION} types",
        "",
        "Generated by `python -m genomeos.lang.grammar` from the parser's block table and the IR",
        "dataclasses; `tests/test_grammar.py` fails when this file is out of date. Prose and examples",
        "live in BIOLANG-v0.1.md, v0.2.md and v0.3.md.",
        "",
        "## Lexical rules",
        "",
        "- A file is directives followed by blocks. `#` starts a comment. Blocks are",
        "  `kind header { props }`; props are `key: value`, one per line or `;`-separated; a block may sit",
        "  on one line.",
        "- Only `transcript` nests (inside `gene`). Repeatable keys: "
        + ", ".join(sorted(_parser._REPEATABLE))
        + ".",
        "- Times take a unit: min, h, d, wk, yr. Loci are `chrN:start-end` with an optional strand.",
        "- `when` clauses: `k = v, k = v`; `v` may be `any`, `absent`, alternatives `a|b`, or a comparison",
        "  `>=n` `<=n` `>n` `<n`. `unknown` states that the context was not recorded: it matches no",
        "  context, so a rule gated on it runs in no cell rather than in every cell.",
        "- Current behaviour: every block that has a `when` (decisions, timers, rules, events,",
        "  competence windows, commitments and signals) reads it with the one matcher,",
        "  `genomeos.ir.model.matches`, so `any`, `absent`, `a|b`, the comparisons and `unknown` mean",
        "  the same on each.",
        "- `!=` is refused: `k != v` is a parse error naming the clause. Write the values that do match",
        "  instead (`a|b`, `absent`, or a comparison); what `!=` would mean on a missing key is",
        "  undecided, and no program in the repo uses it.",
        "- Every block accepts " + " and ".join(f"`{k}: {v[0]}`" for k, v in COMMON.items()) + ".",
        "",
        "## Directives",
        "",
    ]
    out += [f"- `{v}`" for v in DIRECTIVES.values()]
    out += ["", "## Blocks", ""]
    for kind in _parser._KINDS:
        b = BLOCKS[kind]
        out.append(f"### `{kind}`")
        out.append("")
        out.append(f"`{b['header']}`")
        out.append("")
        if b.get("props"):
            out.append("| property | form | meaning |")
            out.append("|---|---|---|")
            for k, (form, meaning) in b["props"].items():
                out.append(f"| `{k}` | `{form}` | {meaning} |")
            out.append("")
        if b.get("nested"):
            out.append(f"May contain `{b['nested']}` blocks.")
            out.append("")
    out += [
        "## BioIR types",
        "",
        "Dataclasses in `genomeos.ir`; `Module.to_dict()` / `from_dict()` round-trip them as JSON",
        f"(`bioir_version` {VERSION}).",
        "",
        "A confidence a block leaves out is `UNSTATED`: 0.0 in any calculation, told apart from a stated",
        "0.0 by `confidence_stated()`. BioIR JSON writes it as `null` (a stated 0.0 stays `0.0`) and marks",
        f"the file `{ir.RECORDS_UNSTATED}: true`; a file without that mark predates the difference.",
        "",
    ]
    for name in IR_TYPES:
        cls = getattr(ir, name)
        fields = [f"`{f.name}`" for f in dataclasses.fields(cls)]
        out.append(f"- **{name}**: " + ", ".join(fields))
    out.append("")
    return "\n".join(out)


def main(argv: list[str] | None = None) -> int:
    problems = check()
    for p in problems:
        print("grammar:", p, file=sys.stderr)
    sys.stdout.write(render())
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
