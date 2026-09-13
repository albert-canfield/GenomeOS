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

VERSION = "0.3"

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
            "half_life": ("hours", "protein half-life"),
            "sequence": ("residues", "amino-acid sequence"),
            "accession": ("text", "UniProt accession"),
            "isoforms": ("Id, Id", "isoform accessions"),
            "domains": ("name, name", "domain names"),
            "pathways": ("Id, Id", "Reactome pathway ids"),
            "interactions": ("Id, Id", "interaction partners"),
            "structures": ("Id, Id", "PDB ids or AF- AlphaFold ids"),
        },
    },
    "region": {
        "header": "region <Id> { ... }",
        "props": {"locus": ("locus", ""), "role": ("text | unknown", "")},
    },
    "element": {
        "header": "element <Id> { ... }",
        "props": {
            "class": ("promoter | enhancer | insulator | open_chromatin", "regulatory element class"),
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
            "threshold": ("number", "source level giving half-maximal effect"),
            "hill": ("number", "cooperativity"),
            "when": ("k = v, k = v", "context in which the rule applies"),
            "id": ("text", "explicit rule id"),
        },
    },
    "param": {"header": "param <name> = <number> [unit] { ... }", "props": {}},
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
                "count|deaths|fate T|type T|lineage L at N unit (in a..b | = n | >= n | <= n)",
                "repeatable",
            ),
            "reference": ("name", "ground truth to diff against"),
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
            "fraction": ("number", "populations: share (divide may exceed 1)"),
        },
    },
    "experiment": {
        "header": "experiment <Id> { ... }",
        "props": {
            "knockout": ("F, F", "factors never present; signals by id, ligand or receptor never sent"),
            "add": ("F, F", "factors added to the first cell"),
            "environment": ("k = v, ...", "overrides"),
            "until": ("time", "horizon; omitted = the last stage"),
            "expect": ('"text"', "the published phenotype"),
            "assert": ("assert grammar", "repeatable; checked on the mutant"),
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
}
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
        "  `>=n` `<=n` `>n` `<n`.",
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
