# SPDX-License-Identifier: AGPL-3.0-or-later
"""Build the measured gene-rate table the R3 bridge needs, under the registration in ff2062a.

    uv run python scripts/build_gene_rates.py [--no-save]

Reads the cached sources named in `attribution.bridge.MEASURED_RATE_SOURCES`, carries each mouse
measurement to a human symbol through an MGI homology class holding exactly one gene of each species,
and writes one row per human symbol: the total transcription rate T in molecules/(cell*h), the mRNA
and protein half-lives in hours, the translation rate constant, and the measured mRNA copy number
(held for the acceptance check, never supplied as a parameter). Every rate is mouse fibroblast; the
result says so on every row's behalf in `species`.

It also runs the registered species falsifier once: Schofield et al. 2018 measured transcript
half-lives in mouse fibroblasts and in human K562 by one method, and the Spearman correlation of the
two over the same 1:1 homology classes decides whether the mouse constants are declared transferable.

No network, no AlphaGenome request: the three inputs are already in data/cache/rates/.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import pandas as pd

from genomeos import manifest as mf
from genomeos.attribution import bridge
from genomeos.results import save_result

#: rate source key -> the species and cell every row from it carries. One entry today; the survey in
#: `bridge.HUMAN_RATE_SOURCE_SURVEY` says why there is no human second entry.
SPECIES_CELL = {"schwanhausser2011": "Mus musculus, NIH 3T3 fibroblast"}

CACHE = Path("data/cache/rates")
SCHWAN = CACHE / "schwanhausser2011_TableS3.xls"
SCHOFIELD = CACHE / "schofield2018_TableS2_halflives.xlsx"
MGI = CACHE / "MGI_HOM_MouseHumanSequence.rpt"

COLS = {
    "transcription_rate": "transcription rate (vsr) average [molecules/(cell*h)]",
    "mrna_half_life_h": "mRNA half-life average [h]",
    "protein_half_life_h": "Protein half-life average [h]",
    "translation_rate": "translation rate constant (ksp) average [molecules/(mRNA*h)]",
    "mrna_copies": "mRNA copy number average [molecules/cell]",
}
#: The source's rows are protein groups, so one row can name several mouse symbols (shared peptides,
#: histone clusters, paralogues). A row is assigned to a human gene only when EVERY mouse symbol on it
#: that has a 1:1 homology class points at the same human gene; a row whose symbols point at two or
#: more human genes names no single gene and is dropped, never assigned to its first symbol. A human
#: symbol receiving measurements from two different rows is dropped as well rather than averaged. The
#: stricter variant (keep only rows naming exactly one mouse symbol) is reported beside it as
#: `single_symbol_rows_only`, so the cost of the choice is visible.
AMBIGUITY = (
    "a row whose mouse symbols map to more than one human gene, or a human symbol receiving more than"
    " one row, is dropped; never assigned to the first symbol and never averaged"
)


def homology() -> dict[str, str]:
    """mouse symbol -> human symbol, only within a class holding exactly one gene of each species."""
    d = pd.read_csv(MGI, sep="\t", dtype=str)
    key, org, sym = "DB Class Key", "Common Organism Name", "Symbol"
    d = d[[key, org, sym]].dropna()
    out: dict[str, str] = {}
    for _, grp in d.groupby(key):
        mouse = grp[grp[org].str.startswith("mouse")][sym].tolist()
        human = grp[grp[org].str.startswith("human")][sym].tolist()
        if len(mouse) == 1 and len(human) == 1:
            out[mouse[0]] = human[0]
    return out


def rate_table(h2m: dict[str, str]) -> tuple[dict, dict]:
    """One row per human symbol, and the counts of what was dropped on the way."""
    d = pd.read_excel(SCHWAN)
    counts = {"rows": len(d)}
    d = d[d["Gene Names"].notna()].copy()
    counts["with_a_symbol"] = len(d)
    counts["single_symbol_rows_only"] = int((~d["Gene Names"].astype(str).str.contains(";")).sum())

    def target(names: str) -> str | None:
        hits = {h2m[m] for m in names.split(";") if m in h2m}
        return hits.pop() if len(hits) == 1 else None

    d["human"] = d["Gene Names"].astype(str).map(target)
    counts["rows_naming_several_human_genes"] = int(
        sum(len({h2m[m] for m in str(n).split(";") if m in h2m}) > 1 for n in d["Gene Names"])
    )
    mapped = d[d["human"].notna()]
    counts["mapped_1to1_to_human"] = len(mapped)
    rated = mapped[mapped[COLS["transcription_rate"]].notna()]
    counts["with_a_transcription_rate"] = len(rated)
    rows: dict[str, dict] = {}
    duplicated = {h for h, n in rated["human"].value_counts().items() if n > 1}
    counts["human_symbols_with_several_mouse_rows"] = len(duplicated)
    for _, r in rated.iterrows():
        h = r["human"]
        if h in duplicated:
            continue
        row = {k: (None if pd.isna(r[c]) else round(float(r[c]), 4)) for k, c in COLS.items()}
        row["mouse_symbol"] = str(r["Gene Names"])
        # provenance per ROW, registered in bridge.RATE_ROW_PROVENANCE: a second source's rows must be
        # distinguishable from these ones, and a human rate must never be pooled with a borrowed one.
        # Today every row says the same thing, and that is the point: it says it per row, not in a note
        row["source"] = "schwanhausser2011"
        row["species_cell"] = SPECIES_CELL["schwanhausser2011"]
        rows[h] = row
    counts["human_genes_in_the_table"] = len(rows)
    return rows, counts


def _spearman(xs: list[float], ys: list[float]) -> float:
    n = len(xs)

    def ranks(v: list[float]) -> list[float]:
        order = sorted(range(n), key=lambda i: v[i])
        out = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j + 1 < n and v[order[j + 1]] == v[order[i]]:
                j += 1
            r = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                out[order[k]] = r
            i = j + 1
        return out

    rx, ry = ranks(xs), ranks(ys)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry, strict=True))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    return num / den if den else float("nan")


def species_falsifier(h2m: dict[str, str]) -> dict:
    """The registered single test: mouse fibroblast against human K562 half-lives, one method."""
    k = pd.read_excel(SCHOFIELD, sheet_name="Table S2_K562")
    m = pd.read_excel(SCHOFIELD, sheet_name="Table S2_MEF")
    human = dict(zip(k["transcript"].astype(str), k["mean_half_life"].astype(float), strict=True))
    pairs = []
    for mouse_sym, hl in zip(m["transcript"].astype(str), m["mean_half_life"].astype(float), strict=True):
        h = h2m.get(mouse_sym)
        if h is not None and h in human:
            pairs.append((hl, human[h]))
    rho = _spearman([a for a, _ in pairs], [b for _, b in pairs]) if len(pairs) > 2 else float("nan")
    threshold = 0.5
    return {
        "test": bridge.RATE_FALSIFIER["species_transfer"],
        "pairs": len(pairs),
        "spearman": round(rho, 4),
        "threshold": threshold,
        "mouse_median_half_life_h": round(float(m["mean_half_life"].median()), 3),
        "human_median_half_life_h": round(float(k["mean_half_life"].median()), 3),
        "verdict": "transferable"
        if rho >= threshold
        else "NOT transferable: the mouse constants are borrowed",
    }


def manifest() -> dict:
    return {
        "sources": [
            {"accession": f"{k}: {v}", "version": "as cached 2026-09-28"}
            for k, v in bridge.MEASURED_RATE_SOURCES.items()
            if k != "sender_milo_2021"
        ],
        "inputs": [mf.input_entry(p) for p in (SCHWAN, SCHOFIELD, MGI) if p.exists()],
        "assembly": mf.NOT_APPLICABLE + ": gene symbols, no coordinate is used or reported",
        "coordinates": mf.NOT_APPLICABLE + ": no interval is reported",
        "parameters": {
            "ambiguity_rule": AMBIGUITY,
            "units": bridge.MEASURED_RATE_UNITS,
            "species": bridge.RATE_SPECIES,
            "declared_strength": bridge.DECLARED_STRENGTH,
            "internal_consistency_disclosed": bridge.RATE_SOURCE_INTERNAL_CONSISTENCY,
        },
        "exclusions": [AMBIGUITY, "a row with no measured transcription rate"],
        "partitions": mf.NOT_APPLICABLE + ": a parameter table, not an evaluation",
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--no-save", action="store_true")
    args = ap.parse_args(argv)
    h2m = homology()
    rows, counts = rate_table(h2m)
    counts["mgi_one_to_one_classes"] = len(h2m)
    ts = sorted(r["transcription_rate"] for r in rows.values())
    median = ts[len(ts) // 2]
    out = {
        "result": "gene_rates",
        "review_item": "R3",
        "registered_in": "ff2062a",
        "species": bridge.RATE_SPECIES,
        "units": bridge.MEASURED_RATE_UNITS,
        "counts": counts,
        "borrowed_median_transcription_rate": round(median, 4),
        "transcription_rate_range": [ts[0], ts[-1]],
        "species_falsifier": species_falsifier(h2m),
        "row_provenance": bridge.RATE_ROW_PROVENANCE,
        "rates_by_source": {k: sum(1 for r in rows.values() if r["source"] == k) for k in SPECIES_CELL},
        "human_measured_rates": sum(1 for r in rows.values() if r["species_cell"].startswith("Homo")),
        "human_rate_survey": {
            "verdict": bridge.HUMAN_RATE_SURVEY_VERDICT,
            "falsifier": bridge.HUMAN_RATE_SURVEY_FALSIFIER,
            "construction_not_taken": bridge.HUMAN_RATE_CONSTRUCTION_NOT_TAKEN,
            "candidates": bridge.HUMAN_RATE_SOURCE_SURVEY,
        },
        "rates": rows,
    }
    if not args.no_save:
        print("saved", save_result("gene_rates", out, manifest=manifest()))
    print({k: v for k, v in out.items() if k != "rates" and k != "units"})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
