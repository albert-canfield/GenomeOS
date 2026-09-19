# SPDX-License-Identifier: AGPL-3.0-or-later
# ruff: noqa: E501  (BioLang templates keep one fact per line)
"""The stage 1 gate of BioLang v0.4 (docs/BIOLANG-v0.4-ECONOMY.md §4.4): does a located program put
the mitochondrial proteome where it belongs, and does breaking a route change the prediction?

Reads MitoCarta3.0 (Rath et al. 2021; the human inventory table, downloaded once into
data/knowledge/mitocarta, git-ignored), the chrM gene models already committed with GENCODE v50,
and UniProt's transit-peptide keyword from the local protein cache when it is there. Generates a
located BioLang program in memory for all 1,136 MitoCarta genes and runs it on the engine's located
runtime under four conditions: wild type, import closed, mitochondrial genome removed (rho0), and
one subunit's presequence removed. Writes data/results/located_mitochondrion.json and the committed
OXPHOS program data/organisms/human/oxphos.bio, which tests itself under `bio test`.

    uv run python scripts/located_mitochondrion.py [--hours 24] [--skip-full]
"""

from __future__ import annotations

import argparse
import gzip
import html
import json
import re
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MITOCARTA_DIR = ROOT / "data" / "knowledge" / "mitocarta"
MITOCARTA_URL = "https://personal.broadinstitute.org/scalvo/MitoCarta3.0/human.mitocarta3.0.html"
PROTEIN_CACHE = ROOT / "data" / "knowledge" / "proteins"
GENCODE_CHRM = ROOT / "data" / "results" / "gencode_v50_chr21_chrM.gff3.gz"
OXPHOS_PROGRAM = ROOT / "data" / "organisms" / "human" / "oxphos.bio"
COMPLEXES = {
    "CI": "Complex I",
    "CII": "Complex II",
    "CIII": "Complex III",
    "CIV": "Complex IV",
    "CV": "Complex V",
}


def mitocarta_rows() -> list[dict]:
    """Symbol, evidence, sub-compartment and pathways for every MitoCarta3.0 human gene."""
    path = MITOCARTA_DIR / "human.mitocarta3.0.html"
    if not path.exists():
        MITOCARTA_DIR.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(MITOCARTA_URL, path)  # noqa: S310  (a fixed https URL)
    rows = []
    for tr in re.findall(r"<tr>(.*?)</tr>", path.read_text(errors="replace"), flags=re.S):
        cells = [
            html.unescape(re.sub(r"<[^>]+>", " ", c)).replace("\xa0", " ")
            for c in re.findall(r"<td[^>]*>(.*?)</td>", tr, flags=re.S)
        ]
        if len(cells) < 7:
            continue
        rows.append(
            {
                "symbol": cells[0].strip(),
                "evidence": " ".join(cells[4].split()),
                "subcompartment": cells[5].strip(),
                "pathways": [" ".join(p.split()) for p in cells[6].split(";") if p.strip()],
            }
        )
    return rows


def uniprot_keywords(symbol: str) -> tuple[str, list[str]] | None:
    """(accession, keywords) from the local UniProt cache, or None when the record is not cached."""
    path = PROTEIN_CACHE / f"{symbol}.json"
    if not path.exists():
        return None
    items = ((json.loads(path.read_text()).get("sections") or {}).get("identity") or {}).get("items") or {}
    return (items.get("accession", ""), list(items.get("keywords") or [])) if items else None


def chrm_loci() -> dict[str, str]:
    """GENCODE v50 protein-coding gene loci on chrM, 0-based half-open as BioLang writes them."""
    loci = {}
    with gzip.open(GENCODE_CHRM, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[0] != "chrM" or f[2] != "gene" or "gene_type=protein_coding" not in f[8]:
                continue
            name = re.search(r"gene_name=([^;]+)", f[8]).group(1)
            loci[name] = f"chrM:{int(f[3]) - 1}-{f[4]}({f[6]})"
    return loci


def signal_of(row: dict) -> tuple[str, str, str, float]:
    """(signal, evidence kind, source, confidence) for a nuclear-encoded mitochondrial protein."""
    uni = uniprot_keywords(row["symbol"])
    if uni and "Transit peptide" in uni[1]:
        return "presequence", "curated", f"UniProt {uni[0]} transit peptide", 0.9
    if "targetP" in row["evidence"]:
        return "presequence", "predicted", "TargetP, as reported by MitoCarta3.0", 0.6
    return "internal", "inferred", "no presequence annotated; internal signal (Wiedemann & Pfanner 2017)", 0.4


CELL = """
compartment Extracellular { evidence: curated "Alberts et al., Molecular Biology of the Cell 6e"; confidence: 0.9 }
compartment PlasmaMembrane { parent: Extracellular; membrane: yes; evidence: curated "Alberts MBoC 6e"; confidence: 0.9 }
compartment Cytosol { parent: PlasmaMembrane; volume: 0.54; translation: yes
  evidence: experimental "Alberts MBoC 6e Table 12-1, hepatocyte volume fractions"; confidence: 0.6 }
compartment Nucleus { parent: Cytosol; volume: 0.06; genome: nuclear
  evidence: experimental "Alberts MBoC 6e Table 12-1, hepatocyte volume fractions"; confidence: 0.6 }
compartment Mitochondrion { parent: Cytosol; volume: 0.22; genome: chrM; translation: yes
  evidence: experimental "Alberts MBoC 6e Table 12-1, hepatocyte volume fractions"; confidence: 0.6 }
compartment ER { parent: Cytosol; volume: 0.09
  evidence: experimental "Alberts MBoC 6e Table 12-1, rough ER cisternae (hepatocyte)"; confidence: 0.6 }

# capacities are set far above demand: stage 1 tests routes, not rates, and says so
transport NuclearPore { from: Nucleus; to: Cytosol; cargo: mRNA; capacity: 1e6; affinity: 1e5
  evidence: inferred "route: mRNA export through the nuclear pore (Alberts MBoC 6e); capacity not measured"; confidence: 0.3 }
transport TOM_TIM23 { from: Cytosol; to: Mitochondrion; cargo: signal = presequence; capacity: 1e6; affinity: 1e5
  evidence: inferred "route: presequence pathway (Wiedemann & Pfanner 2017, Annu Rev Biochem 86:685); capacity not measured"; confidence: 0.3 }
transport TOM_internal { from: Cytosol; to: Mitochondrion; cargo: signal = internal; capacity: 1e6; affinity: 1e5
  evidence: inferred "routes: carrier, beta-barrel and IMS pathways lumped (Wiedemann & Pfanner 2017); capacity not measured"; confidence: 0.3 }

param mrna_half_life = 0.693 h { evidence: inferred "time scale only; stage 1 claims places, not amounts"; confidence: 0.1 }
param protein_half_life = 0.139 h { evidence: inferred "time scale only; stage 1 claims places, not amounts"; confidence: 0.1 }
"""


def complex_members(rows: list[dict]) -> dict[str, list[str]]:
    return {
        cid: [r["symbol"] for r in rows if any(p.endswith(f"{cid} subunits") for p in r["pathways"])]
        for cid in COMPLEXES
    }


def program(rows: list[dict], name: str) -> tuple[str, dict]:
    """A located BioLang program for these MitoCarta genes; returns (text, counts by evidence)."""
    loci = chrm_loci()
    out = [f"module {name}", CELL]
    counts: dict[str, int] = {}
    for r in rows:
        s = r["symbol"]
        if s.startswith("MT-"):
            out.append(
                f"gene {s} {{ locus: {loci[s]}; basal: 1; produces: {s}p\n"
                f'  evidence: curated "GENCODE v50; translated and matched to UniProt, 13 of 13 (translation_vs_uniprot_chrM)"; confidence: 0.95 }}'
            )
            out.append(
                f"protein {s}p {{ location: Mitochondrion\n"
                f'  evidence: curated "MitoCarta3.0 (Rath et al. 2021, NAR 49:D1541)"; confidence: 0.9 }}'
            )
            counts["mtDNA"] = counts.get("mtDNA", 0) + 1
            continue
        sig, kind, src, conf = signal_of(r)
        counts[f"{sig} ({kind})"] = counts.get(f"{sig} ({kind})", 0) + 1
        out.append(
            f"gene {s} {{ location: Nucleus; basal: 1; produces: {s}p\n"
            f'  evidence: curated "MitoCarta3.0: nuclear-encoded"; confidence: 0.9 }}'
        )
        out.append(
            f"protein {s}p {{ location: Mitochondrion; signals: {sig}\n"
            f'  evidence: {kind} "MitoCarta3.0 location ({r["subcompartment"]}); signal: {src}"; confidence: {conf} }}'
        )
    return "\n".join(out) + "\n", counts


def complexes(rows: list[dict]) -> str:
    """Each OXPHOS complex as a readout assembled from its MitoCarta3.0 subunits, all required."""
    out = ["", "# every listed subunit is required, tissue paralogues included: a stated simplification"]
    for cid, members in complex_members(rows).items():
        src = f"MitoCarta3.0 OXPHOS > {COMPLEXES[cid]} > {cid} subunits"
        out.append(f'protein {cid} {{ location: Mitochondrion; evidence: curated "{src}"; confidence: 0.8 }}')
        out += [f'rule {s}p binds {cid} {{ evidence: curated "{src}"; confidence: 0.8 }}' for s in members]
    return "\n".join(out) + "\n"


def full_gate(rows: list[dict], hours: float) -> dict:
    """Gates 1 and 2 on the whole MitoCarta3.0 inventory."""
    from genomeos.lang import parse
    from genomeos.runtime.located import LocatedRuntime

    text, counts = program(rows, "human.mitocarta")
    module = parse(text)
    mt = [r["symbol"] for r in rows if r["symbol"].startswith("MT-")]
    nuclear = [r["symbol"] for r in rows if not r["symbol"].startswith("MT-")]
    runs, seconds = {}, {}
    for label, ko in (
        ("wild_type", set()),
        ("import_closed", {"TOM_TIM23", "TOM_internal"}),
        ("presequence_route_closed", {"TOM_TIM23"}),
        ("internal_route_closed", {"TOM_internal"}),
    ):
        t0 = time.time()
        runs[label] = LocatedRuntime(module, knockouts=ko).run(hours=hours, dt=0.1, record_every=10**9)
        seconds[label] = round(time.time() - t0, 1)
    wt = runs["wild_type"]
    lay = LocatedRuntime(module).layout
    in_mito = lambda res, s: res.final(f"{s}p@Mitochondrion") > 1e-6  # noqa: E731
    mt_ok = [
        s
        for s in mt
        if wt.where(f"{s}p") == ["Mitochondrion"]
        and lay.synthesis.get(f"{s}p") == ["Mitochondrion"]
        and in_mito(runs["import_closed"], s)
    ]
    return {
        "module": module,
        "runs": runs,
        "seconds": seconds,
        "counts": counts,
        "mt": mt,
        "nuclear": nuclear,
        "mt_ok": mt_ok,
        "in_mito": in_mito,
    }


def summarise_full(g: dict) -> dict:
    runs, nuclear, in_mito = g["runs"], g["nuclear"], g["in_mito"]
    signal = {p.id[:-1]: (p.signals or [""])[0] for p in g["module"].proteins() if p.id.endswith("p")}
    reach = {k: sum(1 for s in nuclear if in_mito(r, s)) for k, r in runs.items()}
    by_route = {
        k: {
            sig: sum(1 for s in nuclear if signal[s] == sig and in_mito(runs[k], s))
            for sig in ("presequence", "internal")
        }
        for k in ("presequence_route_closed", "internal_route_closed")
    }
    return {
        "gate_1_mtDNA_proteins": {
            "proteins": len(g["mt"]),
            "born_and_kept_in_mitochondrion_without_transport": len(g["mt_ok"]),
            "pass": len(g["mt_ok"]) == len(g["mt"]) == 13,
        },
        "gate_2_nuclear_encoded_import": {
            "proteins": len(nuclear),
            "reach_mitochondrion": reach,
            "by_signal_when_one_route_closed": by_route,
            "stranded_when_import_closed": len(runs["import_closed"].stranded),
            "pass": reach["wild_type"] == len(nuclear) and reach["import_closed"] == 0,
        },
        "signals": g["counts"],
        "seconds": g["seconds"],
        "regime": runs["wild_type"].regime,
    }


OXPHOS_HEADER = """# The five respiratory complexes, located: BioLang v0.4 stage 1 (docs/BIOLANG-v0.4-ECONOMY.md).
# Generated by scripts/located_mitochondrion.py from MitoCarta3.0 (Rath et al. 2021, NAR 49:D1541),
# GENCODE v50 chrM models and UniProt transit-peptide annotation. Two genomes in one cell: the 13
# chrM proteins are made inside the mitochondrion, every nuclear-encoded subunit needs a transport.
# The published test is the rho0 cell (King & Attardi 1989, Science 246:500): without mitochondrial
# DNA, complexes I, III, IV and V have no activity and complex II, entirely nuclear-encoded, keeps it.
"""
OXPHOS_TESTS = """
# test: stranded == 0
# test: MT-ND1p@Mitochondrion final > 0.05
# test: SDHBp@Mitochondrion final > 0.05
# test: CI@Mitochondrion final > 0.05
# test: CII@Mitochondrion final > 0.05
# test: CIII@Mitochondrion final > 0.05
# test: CIV@Mitochondrion final > 0.05
# test: CV@Mitochondrion final > 0.05
"""

OXPHOS_EXPERIMENTS = """
experiment rho0 {
  knockout: chrM
  expect: "no complex I, III, IV or V activity; complex II intact (King & Attardi 1989, Science 246:500)"
  assert: CI@Mitochondrion == 0
  assert: CII@Mitochondrion > 0.05
  assert: CIII@Mitochondrion == 0
  assert: CIV@Mitochondrion == 0
  assert: CV@Mitochondrion == 0
  evidence: experimental "King & Attardi 1989, Science 246:500"; confidence: 0.9
}
experiment import_closed {
  knockout: TOM_TIM23, TOM_internal
  expect: "no nuclear-encoded subunit reaches the mitochondrion, so no complex assembles"
  assert: CII@Mitochondrion == 0
  assert: CV@Mitochondrion == 0
  assert: stranded >= 80
  evidence: inferred "consequence of the routes (Wiedemann & Pfanner 2017); mitochondrial translation's own dependence on import is not modelled in stage 1"; confidence: 0.4
}
experiment sdhb_without_presequence {
  knockout: SDHBp:presequence
  expect: "SDHB stranded in the cytosol; complex II alone is lost"
  assert: CII@Mitochondrion == 0
  assert: CI@Mitochondrion > 0.05
  assert: stranded == 1
  evidence: inferred "a mislocalisation must change the prediction"; confidence: 0.3
}
"""


def write_oxphos(rows: list[dict]) -> tuple[Path, str]:
    members = {s for v in complex_members(rows).values() for s in v}
    subset = [r for r in rows if r["symbol"] in members]
    text, _ = program(subset, "human.oxphos")
    text = OXPHOS_HEADER + text + complexes(rows) + OXPHOS_EXPERIMENTS + OXPHOS_TESTS
    OXPHOS_PROGRAM.parent.mkdir(parents=True, exist_ok=True)
    OXPHOS_PROGRAM.write_text(text)
    return OXPHOS_PROGRAM, text


def oxphos_gate(rows: list[dict], hours: float) -> dict:
    """Gate 3: the rho0 phenotype and a mislocalisation, on the committed OXPHOS program."""
    from genomeos.lang import parse_file
    from genomeos.runtime.located import LocatedRuntime

    module = parse_file(OXPHOS_PROGRAM)
    conditions = {
        "wild_type": set(),
        "rho0": {"chrM"},
        "import_closed": {"TOM_TIM23", "TOM_internal"},
        "sdhb_without_presequence": {"SDHBp:presequence"},
    }
    out = {}
    for label, ko in conditions.items():
        res = LocatedRuntime(module, knockouts=ko).run(hours=hours, dt=0.05, record_every=10**9)
        out[label] = {
            "complexes": {cid: round(res.final(f"{cid}@Mitochondrion"), 6) for cid in COMPLEXES},
            "stranded": [
                f"{r['protein']} (declared {r['declared']}, found in {r['found_in']})" for r in res.stranded
            ][:5],
            "stranded_count": len(res.stranded),
        }
    wt, rho0 = out["wild_type"]["complexes"], out["rho0"]["complexes"]
    out["pass"] = bool(
        all(v > 0 for v in wt.values())
        and rho0["CII"] > 0
        and all(rho0[c] == 0 for c in ("CI", "CIII", "CIV", "CV"))
        and all(v == 0 for v in out["import_closed"]["complexes"].values())
        and out["sdhb_without_presequence"]["complexes"]["CII"] == 0
        and out["sdhb_without_presequence"]["complexes"]["CI"] > 0
        and out["sdhb_without_presequence"]["stranded_count"] == 1
    )
    out["subunits"] = {cid: len(v) for cid, v in complex_members(rows).items()}
    return out


def refusals() -> dict:
    """What the compiler must refuse rather than run: each program here is one modelling error."""
    from genomeos.lang import BioLangError, parse

    cell = "compartment Cytosol { translation: yes }\ncompartment Nucleus { parent: Cytosol; genome: nuclear }\ncompartment Mitochondrion { parent: Cytosol; genome: chrM; translation: yes }\n"
    cases = {
        "gene without a location": "gene X { produces: Xp }\nprotein Xp { location: Cytosol }",
        "protein without a location": "gene X { location: Nucleus; produces: Xp }\nprotein Xp { }",
        "chrM gene declared nuclear": "gene MT-ND1 { locus: chrM:3306-4262(+); location: Nucleus }",
        "nuclear gene with no mRNA export": "gene X { location: Nucleus; produces: Xp }\nprotein Xp { location: Cytosol }",
        "rule across compartments": "gene X { location: Nucleus }\nprotein TF { location: Mitochondrion }\nrule TF activates X { }",
        "transport between non-adjacent places": "transport T { from: Nucleus; to: Mitochondrion; cargo: mRNA; capacity: 1 }",
    }
    out = {}
    for label, body in cases.items():
        try:
            parse(cell + body)
            out[label] = "ACCEPTED (should have been refused)"
        except BioLangError as e:
            out[label] = f"refused: {e}"
    return out


def evidence_share(module) -> dict:
    """How much of a located program is measured or curated, and how much is predicted or guessed."""
    facts = list(module.entities.values()) + list(module.rules) + list(module.parameters.values())
    kinds: dict[str, int] = {}
    for f in facts:
        kinds[f.evidence.kind.value] = kinds.get(f.evidence.kind.value, 0) + 1
    layer = [e for e in module.entities.values() if e.kind in ("compartment", "transport")]
    layer_kinds: dict[str, int] = {}
    for e in layer:
        layer_kinds[e.evidence.kind.value] = layer_kinds.get(e.evidence.kind.value, 0) + 1
    grounded = kinds.get("experimental", 0) + kinds.get("curated", 0)
    return {
        "facts": len(facts),
        "by_evidence": kinds,
        "measured_or_curated_share": round(grounded / len(facts), 3) if facts else None,
        "structure_constructs": layer_kinds,
    }


def main(argv: list[str] | None = None) -> None:
    from genomeos.bio import evaluate, parse_tests
    from genomeos.lang import parse_file
    from genomeos.results import save_result

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--hours", type=float, default=24.0)
    ap.add_argument("--skip-full", action="store_true", help="only the OXPHOS and red-cell programs")
    args = ap.parse_args(argv)
    t0 = time.time()
    rows = mitocarta_rows()
    write_oxphos(rows)
    payload: dict = {"stage": "BioLang v0.4 stage 1: compartments, locations, transports"}
    if not args.skip_full:
        full = full_gate(rows, args.hours)
        payload["full_mitocarta"] = summarise_full(full)
        payload["evidence_full_program"] = evidence_share(full["module"])
    payload["gate_3_rho0_and_mislocalisation"] = oxphos_gate(rows, args.hours)
    red = ROOT / "data" / "organisms" / "human" / "erythrocyte.bio"
    red_module = parse_file(red)
    checks = evaluate(red_module, parse_tests(red.read_text()), hours=args.hours)
    payload["gate_4_no_genome"] = {
        "genes": len(red_module.genes()),
        "checks": [{"test": c["test"], "got": c["got"], "ok": c["ok"]} for c in checks],
        "pass": red_module.genes() == [] and all(c["ok"] for c in checks),
    }
    payload["refusals"] = refusals()
    payload["evidence_oxphos_program"] = evidence_share(parse_file(OXPHOS_PROGRAM))
    payload["evidence_erythrocyte_program"] = evidence_share(red_module)
    payload["seconds"] = round(time.time() - t0, 1)
    payload["sources"] = [
        "MitoCarta3.0, Rath et al. 2021, Nucleic Acids Res 49:D1541 (human inventory, sub-compartments, OXPHOS subunits, TargetP)",
        "UniProtKB/Swiss-Prot keywords (transit peptide), local cache",
        "GENCODE v50 chrM gene models",
        "King & Attardi 1989, Science 246:500 (rho0 cells)",
        "Wiedemann & Pfanner 2017, Annu Rev Biochem 86:685 (import routes)",
        "Alberts et al., Molecular Biology of the Cell 6e, Table 12-1 (compartment volumes)",
    ]
    payload["note"] = (
        "stage 1 claims places, not amounts: levels are arbitrary units, capacities are set far above demand and"
        " every OXPHOS subunit is required (tissue paralogues included); mitochondrial translation does not yet"
        " depend on imported ribosomal proteins, which is a resource dependence for stage 2"
    )
    path = save_result("located_mitochondrion", payload)
    print(
        json.dumps(
            {k: v for k, v in payload.items() if k.startswith(("gate", "full"))}, indent=1, default=str
        )[:4000]
    )
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
