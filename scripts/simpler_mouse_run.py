# SPDX-License-Identifier: AGPL-3.0-or-later
"""Simpler genomes, area I: the domain rule in a MAMMAL -- mouse, the between-gene stratum.

    uv run python scripts/simpler_mouse_run.py fetch       # inputs into the git-ignored store
    uv run python scripts/simpler_mouse_run.py inventory   # data/results/simpler_mouse_inventory.json
    uv run python scripts/simpler_mouse_run.py register    # census + registration, no RNA value read
    uv run python scripts/simpler_mouse_run.py score       # refuses unless the registration is committed

`inventory` answers the lane's first question -- what mouse accessibility, RNA and boundary data
exist, under what licence, at what assembly and stage, and how many BETWEEN-GENE elements change
gene under the boundary rule -- entirely from positions, the ATAC peak calls and the boundaries. If
that count is under `simpler_mouse.STOP_RULE` the inventory IS the result and the verdict is
NO_SOURCE_QUALIFIES. `register` fixes the test before any RNA value is read. `score` runs it once.

docs/SIMPLER-GENOMES.md section 12 holds the reasoning. No model request is sent, no key is read,
no money is spent, and no human reader biosample is touched or widened.
"""

from __future__ import annotations

import datetime as dt
import hashlib
import json
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import simpler  # noqa: E402
from genomeos.attribution import simpler_flysource as fs  # noqa: E402
from genomeos.attribution import simpler_mouse as sm  # noqa: E402
from genomeos.results import save_result  # noqa: E402

INVENTORY = "simpler_mouse_inventory"
REGISTRATION = "simpler_mouse_registration"
RESULT = "simpler_mouse_result"

ENTRY = "scripts/simpler_mouse_run.py"
OWN_CODE = (
    "genomeos/attribution/simpler_mouse.py",
    "scripts/simpler_mouse_run.py",
    "genomeos/attribution/simpler.py",
    "genomeos/attribution/simpler_fly.py",
    "genomeos/attribution/simpler_flynode.py",
    "genomeos/attribution/simpler_flysource.py",
)

BUDGET_BYTES = 2 * 1024**3

PARAMETERS = {
    "assembly": sm.ASSEMBLY,
    "liftover": sm.LIFTOVER,
    "labels": list(sm.LABELS),
    "detection_tpm": sm.DETECTION_TPM,
    "boundary_bin_bp": sm.BIN_BP,
    "primary_boundary_set": sm.PRIMARY_SET,
    "cluster_bp": sm.CLUSTER_BP,
    "bootstrap_reps": sm.BOOTSTRAP,
    "control_draws": sm.CONTROL_DRAWS,
    "rotation_draws": sm.ROTATION_DRAWS,
    "seed": sm.SEED,
    "alpha": sm.ALPHA,
    "stop_rule_min_between_gene_elements_that_change_gene": sm.STOP_RULE,
    "min_label_for_a_per_label_row": sm.MIN_LABEL,
}

EXCLUSIONS = [
    "cCREs whose SCREEN class is PLS, pELS or DNase-H3K4me3 (kept only as the positive control)",
    "cCREs that overlap an ATAC replicated peak in more than one of the twelve labels",
    "cCREs that overlap an ATAC replicated peak in none of the twelve labels",
    "cCREs whose midpoint lies in a protein-coding gene body (reported beside, not the population)",
    "cCREs on unplaced and alternate mm10 contigs, where a gene order and a boundary are undefined",
    "Bonev Table S2's Compartment column: which compartment a domain sits in is not what the rule uses",
    "ENCODE's e14.5 `brain` experiments, because brain contains forebrain, midbrain and hindbrain",
    "the narrowPeak score, p-value, q-value and summit columns: the peak list IS the source's call",
    "every RNA value, until the registration is committed",
]


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ------------------------------------------------------------------------------------------- inputs


def _encode(acc: str, ext: str) -> str:
    return sm.ENCODE_FILE_URL.format(acc=acc, ext=ext)


def _inputs() -> dict[str, dict[str, Any]]:
    """Every file this lane reads, with the URL it came from and what is read out of it."""
    spec: dict[str, dict[str, Any]] = {
        "gencode.vM21.annotation.gff3.gz": {
            "url": sm.GENE_MODELS_URL,
            "reads": "protein_coding gene and transcript records: chrom, start, end, strand, id",
            "licence": sm.GENE_MODELS_LICENCE,
        },
        "mm10-cCREs.bed": {
            "url": sm.ELEMENT_URL,
            "reads": "columns 1-3, 5 and 6: the element span, its id and SCREEN's own class",
            "licence": sm.ELEMENT_LICENCE,
        },
        "bonev2017_TableS2.xlsx": {
            "url": sm.BOUNDARY_URL,
            "reads": (
                "the chrom, start and end columns of all five sheets; the Compartment column is never read"
            ),
            "licence": sm.BOUNDARY_LICENCE,
        },
    }
    for label, (exp, acc) in sorted(sm.ATAC_FILES.items()):
        spec[f"atac_{acc}.bed.gz"] = {
            "url": _encode(acc, "bed.gz"),
            "reads": f"columns 1-3 only: the replicated ATAC peaks of {label} at e14.5",
            "licence": sm.ENCODE_LICENCE,
            "encode_experiment": exp,
            "label": label,
        }
    for label, (exp, acc) in sorted(sm.RNA_FILES.items()):
        spec[f"rna_{acc}.tsv"] = {
            "url": _encode(acc, "tsv"),
            "reads": f"the gene_id and TPM columns only: {label} at e14.5",
            "licence": sm.ENCODE_LICENCE,
            "encode_experiment": exp,
            "label": label,
        }
    for name, (exp, acc, who) in sorted(sm.CTCF_FILES.items()):
        spec[f"ctcf_{acc}.bed.gz"] = {
            "url": _encode(acc, "bed.gz"),
            "reads": f"columns 1-3 only: the conservative-IDR CTCF peaks of {who} at e14.5",
            "licence": sm.ENCODE_LICENCE,
            "encode_experiment": exp,
            "boundary_set": name,
        }
    return spec


def atac_paths() -> dict[str, Path]:
    return {lab: sm.STORE / f"atac_{acc}.bed.gz" for lab, (_e, acc) in sm.ATAC_FILES.items()}


def rna_paths() -> dict[str, Path]:
    return {lab: sm.STORE / f"rna_{acc}.tsv" for lab, (_e, acc) in sm.RNA_FILES.items()}


def ctcf_paths() -> dict[str, Path]:
    return {name: sm.STORE / f"ctcf_{acc}.bed.gz" for name, (_e, acc, _w) in sm.CTCF_FILES.items()}


def declared_inputs() -> list[tuple[str, str]]:
    """(repo-relative path, what is read) for the manifest. RELATIVE, never absolute (d72f959)."""
    out = [(str(sm.FETCH_RECORD), "each input's registered sha256 and byte count")]
    for name, spec in sorted(_inputs().items()):
        out.append((str(sm.STORE / name), spec["reads"]))
    # THE COMMITTED REGISTRATION IS A READER OF THE SCORING RUN AND IS DECLARED AS ONE. `score`
    # opens it to check that the census recomputes to the registered one, so the manifest contract
    # counts it among the files read under data/ -- and quarantined the first write because it was
    # not declared. It is NOT declared for `register`, which WRITES it: a result may not hash itself
    # as its own input. The step is read from argv exactly as `main` dispatches on it.
    if sys.argv[-1] == "score":
        out.append(
            (
                f"data/results/{REGISTRATION}.json",
                "the registration this run is scored against: the design, and the census the run "
                "must recompute before any value is read",
            )
        )
    return out


def fetch() -> None:
    sm.STORE.mkdir(parents=True, exist_ok=True)
    free = shutil.disk_usage(".").free
    files: dict[str, Any] = {}
    moved = 0
    for name, spec in sorted(_inputs().items()):
        path = sm.STORE / name
        if not path.exists():
            req = urllib.request.Request(spec["url"], headers={"User-Agent": "GenomeOS/1.0"})
            with urllib.request.urlopen(req, timeout=1800) as r:
                body = r.read()
            moved += len(body)
            if moved > BUDGET_BYTES:
                sys.exit("refused: over the 2 GB budget")
            path.write_bytes(body)
            downloaded = len(body)
        else:
            downloaded = 0
        files[name] = {
            **spec,
            "bytes": path.stat().st_size,
            "sha256": _sha(path),
            "downloaded_bytes": downloaded,
            "downloaded_now": downloaded > 0,
        }
    record = {
        "fetched_utc": dt.datetime.now(dt.UTC).isoformat(timespec="seconds"),
        "disk_free_bytes_before": free,
        "budget_bytes": BUDGET_BYTES,
        "download_bytes": moved,
        "kept_bytes": sum(v["bytes"] for v in files.values()),
        "files": files,
    }
    sm.FETCH_RECORD.write_text(json.dumps(record, indent=1) + "\n")
    print(json.dumps({k: v for k, v in record.items() if k != "files"}, indent=1))


# ---------------------------------------------------------------------------------- building the rows


def _elements() -> dict[str, Any]:
    """Positions, SCREEN's class and the ATAC peak calls only. NOTHING here opens an RNA file."""
    chroms, tss, names = sm.read_annotation(sm.GFF3)
    ccres = sm.read_ccres(sm.CCRE_FILE)
    peaks, peak_counts = sm.read_atac_peaks(atac_paths())
    label_of, call_funnel = sm.single_label_elements(ccres, peaks, sm.LABELS)
    els, el_funnel = sm.elements(ccres, label_of, chroms)
    rows = simpler.attribute(els, chroms)
    return {
        "chroms": chroms,
        "tss": tss,
        "names": names,
        "rows": rows,
        "funnel": {
            "protein_coding_genes": sum(len(c.genes) for c in chroms.values()),
            "transcript_5p_ends": sum(len(v) for v in tss.values()),
            "atac_peaks_per_label": peak_counts,
            **call_funnel,
            **el_funnel,
            "attributed_rows": len(rows),
        },
        "lengths": sm.chrom_lengths(sm.CCRE_FILE),
    }


def _boundary_sets() -> tuple[dict[str, dict[str, list[int]]], dict[str, Any]]:
    """Every boundary set, author-called: Bonev's five Hi-C domain sheets and ENCODE's three e14.5
    CTCF peak sets. The Hi-C sheets become boundaries through their domain EDGES; a CTCF peak
    becomes a boundary through its span's midpoint, because a peak is already point-like."""
    bounds: dict[str, dict[str, list[int]]] = {}
    summary: dict[str, Any] = {}
    for name, (sheet, expect, who) in sorted(sm.BONEV_SHEETS.items()):
        rows = sm.read_domain_sheet(sm.DOMAIN_FILE, sheet, expect)
        b, s = sm.boundaries_from_domains(rows)
        bounds[name] = b
        summary[name] = {
            **s,
            "kind": "author-called Hi-C domains, read as their edges",
            "sheet": sheet,
            "rows_the_registration_recorded": expect,
            "biosample": who,
        }
    for name, path in sorted(ctcf_paths().items()):
        b, s = sm.boundaries_from_spans(sm.read_bed_spans(path))
        bounds[name] = b
        exp, acc, who = sm.CTCF_FILES[name]
        summary[name] = {
            **s,
            "kind": "author-called CTCF peaks, read as their span midpoints",
            "encode_experiment": exp,
            "encode_file": acc,
            "biosample": who,
        }
    return bounds, summary


def _manifest(inputs: list[tuple[str, str]], partitions: Any) -> dict:
    return {
        "sources": [
            {"accession": sm.ELEMENT_SOURCE, "version": "SCREEN V3, as served on 2026-10-12"},
            {"accession": sm.ACCESSIBILITY_SOURCE, "version": "as served on 2026-10-12"},
            {"accession": sm.RNA_SOURCE, "version": "as served on 2026-10-12"},
            {"accession": sm.BOUNDARY_SOURCE, "version": "as served on 2026-10-12"},
            {"accession": sm.GENE_MODELS, "version": "GENCODE release M21"},
        ],
        "inputs": [mf.input_entry(path, partition=part) for path, part in inputs],
        "assembly": sm.ASSEMBLY,
        "coordinates": {
            "base": 0,
            "interval": "half-open",
            "note": (
                "the cCRE BED and the narrowPeak files are 0-based half-open as served, the GFF3 is "
                "1-based closed, and a boundary is its peak span's midpoint; positions move by at "
                "most 1 bp"
            ),
        },
        "parameters": PARAMETERS,
        "exclusions": EXCLUSIONS,
        "partitions": partitions,
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


def _census() -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    """The position-only census under every boundary set. No RNA value is read on this path."""
    built = _elements()
    bounds, bsummary = _boundary_sets()
    main = sm.attribute(built["rows"], built["chroms"], bounds[sm.PRIMARY_SET], sm.HALF_BIN)
    alternatives = {
        name: sm.attribute(built["rows"], built["chroms"], b, sm.HALF_BIN)
        for name, b in bounds.items()
        if name != sm.PRIMARY_SET
    }
    cen = sm.census(main, alternatives, built["funnel"], sm.LABELS)
    moves = fs.between_gene_moves(built["rows"], built["chroms"], bounds[sm.PRIMARY_SET], sm.HALF_BIN)
    cen["position_only_move_count_cross_check"] = moves
    cen["boundary_sets"] = bsummary
    return cen, built, {"bounds": bounds, "main": main, "alternatives": alternatives}


# ---------------------------------------------------------------------------------------- the entries


def inventory() -> None:
    cen, _built, _parts = _census()
    out = {
        "result": "the mouse inventory, counted from positions, the ATAC peak calls and the "
        "boundaries alone -- before any RNA value is read",
        "date": "2026-10-12",
        "read_before_any_rna_value_was_read": True,
        "question": sm.REGISTRATION["question"],
        "stop_list_of_2026_09_27": sm.REGISTRATION["stop_list_of_2026_09_27"],
        "what_the_two_invertebrates_found": {"worm": sm.WORM, "fly": sm.FLY},
        "source_chosen": {
            "assembly": sm.ASSEMBLY,
            "liftover": sm.LIFTOVER,
            "gene_models": sm.GENE_MODELS,
            "gene_models_licence": sm.GENE_MODELS_LICENCE,
            "element_registry": sm.ELEMENT_SOURCE,
            "element_registry_licence": sm.ELEMENT_LICENCE,
            "accessibility": sm.ACCESSIBILITY_SOURCE,
            "rna": sm.RNA_SOURCE,
            "boundaries": sm.BOUNDARY_SOURCE,
            "boundary_caveat": sm.BOUNDARY_CAVEAT,
            "encode_licence": sm.ENCODE_LICENCE,
            "one_source_for_both_modalities": sm.REGISTRATION["one_source_for_both_modalities"],
            "labels": list(sm.LABELS),
            "atac_files": {k: list(v) for k, v in sorted(sm.ATAC_FILES.items())},
            "rna_files": {k: list(v) for k, v in sorted(sm.RNA_FILES.items())},
            "ctcf_files": {k: list(v) for k, v in sorted(sm.CTCF_FILES.items())},
        },
        "inputs_on_disk": json.loads(sm.FETCH_RECORD.read_text()),
        "census": cen,
        "stop_rule": {
            "rule": sm.REGISTRATION["stop_rule"],
            "between_gene_elements": cen["between_genes"]["elements"],
            "between_gene_elements_that_change_gene": cen["between_gene_elements_that_change_gene"],
            "threshold": sm.STOP_RULE,
            "met": cen["stop_rule_met"],
            "verdict": cen["verdict"],
        },
        "alphagenome_requests": 0,
        "money": "none",
    }
    out["result_manifest"] = _manifest(
        declared_inputs(),
        {"census": "every count from positions, the peak calls and the boundaries"},
    )
    path = save_result(INVENTORY, out)
    print(
        f"{path}: {cen['verdict']}, {cen['between_gene_elements_that_change_gene']} change gene "
        f"of {cen['between_genes']['elements']} between-gene elements (stop rule {sm.STOP_RULE})"
    )


def _committed(rel: str) -> bool:
    try:
        tracked = subprocess.run(
            ["git", "ls-files", "--error-unmatch", rel], cwd=ROOT, capture_output=True, text=True
        )
        if tracked.returncode != 0:
            return False
        diff = subprocess.run(["git", "diff", "--quiet", "HEAD", "--", rel], cwd=ROOT)
        return diff.returncode == 0
    except OSError:
        return False


def register() -> None:
    cen, _built, _parts = _census()
    if not cen["stop_rule_met"]:
        sys.exit(
            f"refused: {cen['between_gene_elements_that_change_gene']} between-gene elements change "
            f"gene, under the stop rule's {sm.STOP_RULE}. The inventory is the result."
        )
    inv = f"data/results/{INVENTORY}.json"
    if not _committed(inv):
        sys.exit(f"refused: {inv} is not committed and unchanged; commit the inventory first")
    reg = dict(sm.REGISTRATION)
    bg = cen["between_genes"]
    reg["power_in_advance"] = (
        f"{bg['differ_from_nearest_tss']:,} of {bg['elements']:,} population elements change gene "
        f"under the primary boundary set, across {bg['clusters_1mb']:,} 1-Mb clusters: about "
        f"{bg['mde_points_unclustered_vs_nearest_tss']:.2f} points at 80% power, two-sided 5%, "
        f"before clustering. The fly's section 11 had 216 of 3,265 and about 1.3 points."
    )
    reg["stop_rule_counts"] = {
        name: alt["differ_from_nearest_tss"]
        for name, alt in sorted(cen["alternative_boundary_sets_between_genes"].items())
    } | {sm.PRIMARY_SET: bg["differ_from_nearest_tss"]}
    out = {
        "result": "the registration of the mammalian between-gene test, written before any RNA "
        "value was read",
        "date": "2026-10-12",
        "written_before_any_rna_value_was_read": True,
        "registration": reg,
        "census": cen,
        "inventory_commit_checked": True,
        "reading_template_is_committed_in": "genomeos/attribution/simpler_mouse.reading",
        "alphagenome_requests": 0,
        "money": "none",
    }
    out["result_manifest"] = _manifest(
        declared_inputs(), {"registration": "the design, fixed before any RNA value was read"}
    )
    path = save_result(REGISTRATION, out)
    print(f"{path}: registered, {bg['differ_from_nearest_tss']} change gene of {bg['elements']}")


def score() -> None:
    for rel in (*OWN_CODE, f"data/results/{REGISTRATION}.json", f"data/results/{INVENTORY}.json"):
        if not _committed(rel):
            sys.exit(f"refused: {rel} is not committed and unchanged")
    record = json.loads(sm.FETCH_RECORD.read_text())
    for name in sorted(_inputs()):
        path = sm.STORE / name
        got = _sha(path)
        want = record["files"][name]["sha256"]
        if got != want:
            sys.exit(f"refused: {path} is {got}, the registration recorded {want}")
    cen, built, parts = _census()
    if not cen["stop_rule_met"]:
        sys.exit("refused: the stop rule is not met; the inventory is the result")
    reg = json.loads((ROOT / "data/results" / f"{REGISTRATION}.json").read_text())
    if reg["census"]["between_genes"] != cen["between_genes"]:
        sys.exit("refused: the census does not recompute to the registered one")

    expr, expr_summary = sm.read_expression(rna_paths())
    chroms, lengths = built["chroms"], built["lengths"]
    main = parts["main"]
    bounds = parts["bounds"]
    bg = sm.between_genes(main)
    ad = sm.distal(main)
    prom = [r for r in main if r["kind"] == sm.PROMOTER]

    def block(rows: list[dict[str, Any]], with_controls: bool) -> dict[str, Any]:
        out: dict[str, Any] = {
            "boundary_vs_nearest": sm.compare(rows, expr, "boundary", "nearest", sm.LABELS),
            "boundary_vs_gene_structure": sm.compare(rows, expr, "boundary", "domain", sm.LABELS),
            "gene_structure_vs_nearest": sm.compare(rows, expr, "domain", "nearest", sm.LABELS),
        }
        if with_controls:
            out["falsifier"] = sm.matched_distance_control(rows, chroms, expr, sm.LABELS, "boundary")
            out["rotation"] = sm.rotation_null(
                rows, chroms, bounds[sm.PRIMARY_SET], lengths, expr, sm.LABELS, sm.HALF_BIN
            )
        return out

    differ = [r for r in bg if r["boundary_gene"] != r["nearest_gene"]]
    slots = [r[k] for r in bg for k in ("boundary_gene", "nearest_gene", "domain_gene")]
    res: dict[str, Any] = {
        "result": "the domain rule's between-gene stratum in a mammal, run once on the "
        "committed registration",
        "date": "2026-10-12",
        "registration": f"data/results/{REGISTRATION}.json",
        "configuration": {
            "labels": len(sm.LABELS),
            "label_kind": sm.LABEL_KIND,
            "label_noun": sm.LABEL_NOUN,
            "element_source_short": "ENCODE ATAC-seq, embryonic day 14.5, 12 tissues",
            "rna_source_short": "ENCODE polyA plus RNA-seq, the same 12 tissues at the same stage",
            "boundary_source_short": ("CTCF peaks of embryonic-day-14.5 liver (ENCODE, conservative IDR)"),
            "boundaries": cen["boundary_sets"][sm.PRIMARY_SET]["boundaries"],
            "what_is_worse_here": sm.REGISTRATION["what_is_worse_here"],
        },
        "between_genes": block(bg, True),
        "all_distal": block(ad, False),
        "where_the_rules_differ": sm.compare(differ, expr, "boundary", "nearest", sm.LABELS),
        "positive_control": sm.positive_control(prom, chroms, expr, sm.LABELS),
        "alternative_boundary_sets": {
            name: sm.compare(sm.between_genes(rows), expr, "boundary", "nearest", sm.LABELS)
            for name, rows in sorted(parts["alternatives"].items())
        },
        "by_label": {
            lab: sm.compare([r for r in bg if r["tissue"] == lab], expr, "boundary", "nearest", sm.LABELS)
            for lab in sm.LABELS
            if sum(1 for r in bg if r["tissue"] == lab) >= sm.MIN_LABEL
        },
        "rna": expr_summary,
        "attributed_gene_slots": len(slots),
        "attributed_genes_absent_from_the_rna_rows": sum(1 for g in slots if g not in expr),
        "census": cen,
        "alphagenome_requests": 0,
        "money": "none",
    }
    res["outcome"] = res["between_genes"]["boundary_vs_nearest"]["label"]
    res["second_comparison"] = res["between_genes"]["boundary_vs_gene_structure"]["label"]
    res["falsifier"] = res["between_genes"]["falsifier"]["falsifier"]
    # THE READING'S BOUNDARY SLOT IS CORRECTED HERE, BEFORE THE READING IS BUILT. The literal above
    # at `boundary_source_short` is a leftover from an earlier arrangement of this lane in which an
    # ENCODE CTCF set was primary: it names "CTCF peaks of embryonic-day-14.5 liver" while the count
    # beside it is read from the census as the PRIMARY set's 4,796, so the reading would have
    # described the edges of Bonev's `ncx NPC` domains under the liver CTCF set's name -- that set
    # having 49,312 boundaries, not 4,796. `sm.PRIMARY_SOURCE_SHORT` is DERIVED from
    # BONEV_SHEETS[PRIMARY_SET], so the name, the domain count and the sheet cannot drift from the
    # set actually used. The stale literal is left in place and made INERT by this line, rather than
    # deleted, because deleting a line cf05e75 added needs `commit_own.sh --force` and this lane's
    # brief forbids --force; the assertion below is what makes "inert" a fact and not a hope.
    res["configuration"]["boundary_source_short"] = sm.PRIMARY_SOURCE_SHORT
    assert res["configuration"]["boundary_source_short"] == sm.PRIMARY_SOURCE_SHORT
    res["reading"] = sm.reading(res)
    res["result_manifest"] = _manifest(
        declared_inputs(), {"run": "one run on the committed registration, every reader declared"}
    )
    path = save_result(RESULT, res)
    print(f"{path}")
    print(f"OUTCOME {res['outcome']}; SECOND {res['second_comparison']}; FALSIFIER {res['falsifier']}")
    print()
    print(res["reading"])


def main() -> None:
    steps = {"fetch": fetch, "inventory": inventory, "register": register, "score": score}
    if len(sys.argv) != 2 or sys.argv[1] not in steps:
        sys.exit(f"usage: {ENTRY} {{{'|'.join(steps)}}}")
    steps[sys.argv[1]]()


if __name__ == "__main__":
    main()
