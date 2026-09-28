# SPDX-License-Identifier: AGPL-3.0-or-later
"""Read allele fractions of HG002's trio de novo candidates against inherited heterozygous calls
(docs/DATA.md, "Allele fractions of the trio candidates (2026-09-28)", pre-registered before this ran).

Inputs are local and git-ignored. Read depths: GIAB's v4.2.1 HG002 GRCh38 benchmark VCF, FORMAT ADALL
("net allele depths across all datasets"), fetched once from
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/
into data/cache/v421 (156 MB). Assembly: the v5.0q dipcall files under data/cache/q100, as used by
scripts/trio_q100_phase.py. Data: NIST, not subject to US copyright (17 USC 105), provided AS IS;
acknowledgement: National Institute of Standards and Technology, Genome in a Bottle.

Writes data/results/trio_vaf.json: counts, medians and histograms per group, never a position, an
allele or a genotype.
"""

from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.genome import individuals as ind  # noqa: E402
from genomeos.results import save_result  # noqa: E402

DEPTHS = Path("data/cache/v421/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz")
Q100 = Path("data/cache/q100")
VCF = Q100 / "GRCh38_HG2-T2TQ100-V1.1_dipcall-z2k.dip.vcf.gz"
DIP = Q100 / "GRCh38_HG2-T2TQ100-V1.1_dipcall-z2k.dip.bed"
MOSAIC = Q100 / "mergedQ100denovomosaicexclusions_240506_GRCh38_withoverlaprepeats.bed.gz"

# Registered 2026-09-28, before the run (docs/DATA.md).
REGISTERED_CANDIDATES = 1430
FIELD = "ADALL"
MIN_DEPTH = 30
ALPHA = 0.001
LOW_VAF = 0.40
HIGH_VAF = 0.60
C1_MEDIAN = (0.45, 0.55)  # control SNV median VAF
C2_MAX_LOW = 0.05  # control SNV share classed low
C3_MIN_ASSESSED = 0.80  # control SNVs with depth >= MIN_DEPTH
P1_MAX_LOW = 0.25  # one-parent candidates classed low
P2_ALPHA = 0.01  # candidates' low share above the control's, one-sided binomial
SUBCLONAL_MAJORITY = 0.50
CLONAL_MARGIN = 0.05


def verdicts(r: dict) -> dict:
    g = r["groups"]
    cs, one = g.get("control_snv", {}), g.get("candidates_one_parent", {})
    out: dict = {}
    c_assessed = cs.get("assessed", 0) / cs["total"] if cs.get("total") else 0.0
    out["C1_control_snv_median_vaf"] = {
        "value": cs.get("median_vaf"),
        "holds": cs.get("median_vaf") is not None and C1_MEDIAN[0] <= cs["median_vaf"] <= C1_MEDIAN[1],
    }
    out["C2_control_snv_low_share"] = {
        "value": cs.get("low_share"),
        "holds": cs.get("low_share") is not None and cs["low_share"] <= C2_MAX_LOW,
    }
    out["C3_control_snv_assessed"] = {"value": round(c_assessed, 4), "holds": c_assessed >= C3_MIN_ASSESSED}
    out["void"] = not all(out[k]["holds"] for k in list(out))
    low, n = one.get("low", 0), one.get("assessed", 0)
    share = low / n if n else None
    out["P1_one_parent_low_share_at_most"] = {
        "limit": P1_MAX_LOW,
        "value": share,
        "holds": share is not None and share <= P1_MAX_LOW,
    }
    # P2: more candidates classed low than the control's rate predicts
    # the control's rate, weighted by the candidates' own SNV:indel mix
    p0 = 0.0
    for t in ("snv", "indel"):
        nt = g.get(f"candidates_one_parent_{t}", {}).get("assessed", 0)
        p0 += (nt / n if n else 0.0) * (g.get(f"control_{t}", {}).get("low_share") or 0.0)
    p0 = round(p0, 5)
    p2 = ind.binomial_tail(low, n, p0, lower=False) if n and 0 < p0 < 1 else None
    out["P2_low_excess_over_control"] = {
        "control_low_share": p0,
        "p_one_sided": p2,
        "holds": p2 is not None and p2 < P2_ALPHA,
    }
    if share is None:
        out["outcome"] = None
    elif share >= SUBCLONAL_MAJORITY:
        out["outcome"] = "subclonal majority"
    elif share <= p0 + CLONAL_MARGIN:
        out["outcome"] = "clonal"
    else:
        out["outcome"] = "clonal majority with a subclonal minority"
    if n:  # Wilson 95% interval on the low share
        z, ph = 1.96, low / n
        c = (ph + z * z / (2 * n)) / (1 + z * z / n)
        h = z * math.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / (1 + z * z / n)
        out["one_parent_low_share_95"] = [round(c - h, 4), round(c + h, 4)]
    return out


def main() -> int:
    t0 = time.time()
    r = ind.trio_allele_fractions(
        "HG002",
        "HG003",
        "HG004",
        DEPTHS,
        VCF,
        DIP,
        hap1_parent="father",
        other_beds={"nist_de_novo_mosaic_exclusions": MOSAIC},
        field=FIELD,
        min_depth=MIN_DEPTH,
        alpha=ALPHA,
        low_vaf=LOW_VAF,
        high_vaf=HIGH_VAF,
    )
    if r["candidates_total"] != REGISTERED_CANDIDATES:
        print(f"refused: {r['candidates_total']} candidates, registered {REGISTERED_CANDIDATES}")
        return 1
    r["verdicts"] = verdicts(r)
    r["registration"] = "docs/DATA.md, section 'Allele fractions of the trio candidates (2026-09-28)'"
    r["evidence"] = (
        "measured: pooled read allele depths (GIAB v4.2.1 ADALL, all datasets) per candidate against "
        "inherited heterozygous calls of the same child; counts only"
    )
    r["seconds"] = round(time.time() - t0, 1)
    manifest = {
        "sources": [
            {
                "accession": "GIAB HG002 NISTv4.2.1 GRCh38 benchmark VCF",
                "version": "v4.2.1",
                "url": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/"
                "HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
                "licence": "NIST, not subject to US copyright (17 USC 105), AS IS",
            },
            {
                "accession": "GIAB HG002 v5.0q dipcall of T2T HG002 Q100",
                "version": "v5.0q / assembly v1.1",
                "licence": "NIST, not subject to US copyright (17 USC 105), AS IS",
            },
            {"accession": "GIAB HG003 and HG004 benchmark calls and regions", "version": "v4.2.1"},
        ],
        "inputs": [mf.input_entry(p) for p in (DEPTHS, VCF, DIP, MOSAIC)],
        "assembly": "GRCh38",
        "coordinates": "n/a: counts only, no coordinates in the result",
        "parameters": {
            "field": FIELD,
            "min_depth": MIN_DEPTH,
            "alpha": ALPHA,
            "low_vaf": LOW_VAF,
            "high_vaf": HIGH_VAF,
            "window": 10,
            "registered_candidates": REGISTERED_CANDIDATES,
        },
        "exclusions": [
            "chrX, chrY and chrM (the trio is autosomal)",
            "child calls outside any trio member's trusted regions",
            "v4.2.1 already excluded heterozygous calls with net allele fraction below 0.2 or above 0.8",
        ],
        "partitions": "n/a: a descriptive measurement, no model is trained or evaluated",
    }
    save_result("trio_vaf", r, manifest=manifest)
    keep = ("control", "control_snv", "control_indel", "candidates", "candidates_one_parent")
    print(json.dumps({"groups": {k: r["groups"].get(k) for k in keep}, "verdicts": r["verdicts"]}, indent=1))
    print(f"{r['seconds']} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
