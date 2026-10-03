# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""N1b: is the derived per-gene statistic T = X * sqrt(n) CALIBRATED on the non-targeting rows?

WHAT THIS STUDY IS, and the whole of it. One question: does T behave like a standard normal variate
where nothing was perturbed? Nothing about factors. Nothing about whether anything responds. No
factor row is read, at all, by any code in this file.

WHAT THIS STUDY IS NOT. It is NOT a resumption of N1. N1 was closed on 2026-10-01 under its own
registered stop rule, on the external reviewer's advice relayed by Albert, and nothing here reopens
that closure or revisits its finding. N1's committed registration, amendments and result are inputs
to nothing in this file and are not altered by it. This file imports N1's SOURCE constant only so
that the file identity (URL, size, published md5) has ONE copy in the project and cannot drift; it
reads none of N1's results and applies none of N1's thresholds. N1's own gate used |T| >= 3 against
a 1.5x ratio ceiling over ten deciles. This study's rule is different and is registered below on its
own terms: a variance band and a 1.96 tail band over five quintiles. Neither rule is derived from the
other and neither inherits the other's verdict.

OUTCOME EXPOSURE. The non-targeting rows are the NULL. They carry no factor response, by construction
of the experiment: the sgRNA targets nothing. Reading them is therefore "null-calibration exposure,
no outcome row" -- it spends no confirmatory outcome and leaves the project's independent outcome data
untouched, because no row that could show a response is opened.

WHAT A FAIL MEANS, and this is written before the run. If either band is missed, anywhere, the
derivation T = X * sqrt(n) is UNUSABLE and this study says so and stops. It does not then propose a
repair, a wider band or a second statistic. The fallbacks -- the authors' own published per-gene
results, or the single-cell file -- are ALBERT'S decision and not this lane's. Nothing about factors
is read either way.

THE INSTRUMENT'S ASSUMPTIONS. Three, and every figure this module produces carries them:
  (1) X averages exactly `num_cells_filtered` cells. UNCONFIRMED. The producers document X as a
      per-row mean of per-cell gemgroup z-scores but do not document that the divisor is the same
      `num_cells_filtered` the obs column reports. If it is not, sqrt(n) is the wrong scale factor
      and every figure here moves. Nothing in the file settles it.
  (2) The cells pooled into a row are independent.
  (3) The normal approximation holds in the tail the rule reads, at |T| = 1.96.
A figure that fails only because (1) is false would look exactly like a derivation that is wrong, and
this study cannot tell those apart. It says so rather than choosing.
"""

from __future__ import annotations

import hashlib
import io
import math
import urllib.request
from collections.abc import Iterable, Sequence
from typing import Any

from genomeos.attribution.n1_perturb_response import SOURCE

# --- the frozen rules: the pass rule, the strata, the bands ------------------------------------------------

STUDY = "n1b_calibration"

#: The variance band. RELATIVE to the derivation's own nominal value, which is exactly 1.
VAR_LO, VAR_HI = 0.80, 1.25
NOMINAL_VAR = 1.0

#: The tail band. The share of values with |T| strictly greater than TAIL_Z.
TAIL_Z = 1.96
TAIL_LO, TAIL_HI = 0.04, 0.06
NOMINAL_TAIL = math.erfc(TAIL_Z / math.sqrt(2.0))  # two-sided N(0, 1) tail at 1.96: 0.0499958

#: Quintiles of control expression. Five, because the band is read in every one of them.
STRATA = 5

#: A non-targeting row is usable when its cell count is finite and at least one cell.
MIN_CELLS = 1.0

#: The cluster is the non-targeting ROW. Cells are pooled within a row and genes are correlated
#: within a row, so the row is the unit that resamples.
BOOTSTRAP_B = 2000
BOOTSTRAP_SEED = 20261003
BOOTSTRAP_BLOCK = 200  # resamples per matrix product, so no block exceeds a few MB

#: Both bands are read as ratios to a nominal value that is a THEORETICAL CONSTANT, not an empirical
#: control group. Under the three assumptions above, T ~ N(0, 1) exactly, so the reference is 1 for the
#: variance and 0.0499958 for the tail, each with no sampling error of its own. Relative and absolute
#: therefore coincide here, and the bands are stated relative because that is what a reader checks.
#: There is NO control group in this study and none is possible: the null IS the reference, and a
#: second null would be another draw from the same rows.
CONTROL_KIND = "theoretical nominal value, not an empirical control group"

CONSTANTS: dict[str, Any] = {
    "study": STUDY,
    "statistic": "T = X * sqrt(num_cells_filtered), per (non-targeting row, gene)",
    "var_band": [VAR_LO, VAR_HI],
    "nominal_var": NOMINAL_VAR,
    "tail_z": TAIL_Z,
    "tail_band": [TAIL_LO, TAIL_HI],
    "nominal_tail": NOMINAL_TAIL,
    "strata": STRATA,
    "min_cells": MIN_CELLS,
    "bootstrap_b": BOOTSTRAP_B,
    "bootstrap_seed": BOOTSTRAP_SEED,
    "cluster": "non-targeting row",
    "control_kind": CONTROL_KIND,
}

#: The summary that DECIDES, named before the run. "The per-gene empirical variance of T" is a
#: quantity per gene; a band is read against one number, so the summary is named here and only this
#: one decides. The two others are reported beside it and decide nothing.
PRIMARY_VAR_SUMMARY = "median over genes of the per-gene variance of T about the gene's own mean, ddof=1"
SECONDARY_VAR_SUMMARIES = (
    "mean over genes of the same per-gene variance",
    "pooled variance of T over all (row, gene) values in the stratum",
    "mean over genes of the per-gene SECOND MOMENT about zero, E[T^2] = variance + mean^2, which is "
    "the quantity the tail share actually depends on",
)
PRIMARY_TAIL_SUMMARY = "share of all (row, gene) values in the stratum with |T| > 1.96"
SECONDARY_TAIL_SUMMARIES = ("median over genes of the per-gene share with |T| > 1.96",)

PASS_RULE = (
    "Over the usable non-targeting rows, in EVERY ONE of the five control-expression quintiles and "
    f"overall: the primary variance summary lies in [{VAR_LO}, {VAR_HI}] times the nominal 1, AND the "
    f"primary tail summary lies in [{TAIL_LO}, {TAIL_HI}]. Both must hold; either failing, in any one "
    "stratum, fails the whole rule. A stratum whose figure cannot be formed fails it too, and says why."
)

INCONCLUSIVE_WORDING = "the data cannot tell"

#: WHAT A PASS LICENSES, and it is narrower than a pass will read to someone skimming. Registered
#: before any byte of X on the coordinator's reading of the producers' normalisation, which is right
#: and which this lane had not seen: Replogle's per-cell z-scores are computed per gemgroup AGAINST
#: THE NON-TARGETING CELLS -- the mean and SD of the control population. So every non-targeting row
#: is a per-guide SUBSET of the very reference that defined the scale. A pass therefore tests the
#: sqrt(n) scaling, cell independence and normality WITHIN the reference population. It does not show
#: that T is calibrated for PERTURBED rows, whose cell counts, variances and knockdown-induced shifts
#: all differ. `decide` cannot emit a verdict without this sentence attached, so the caveat travels
#: with the figures rather than sitting once in a file nobody re-reads.
PASS_LICENCE = (
    "a PASS licenses T only as calibrated on non-targeting subsets of the normalisation reference; "
    "its transfer to perturbed rows is an ASSUMPTION, and step 2, if it ever runs, states it again "
    "beside every perturbed-row figure"
)

#: What this study cannot establish, written before the run. None of these is a band and none of them
#: is a reason to widen one; they bound what a pass MEANS, not what counts as one.
CANNOT_ESTABLISH = (
    PASS_LICENCE,
    "whether X's divisor is the num_cells_filtered the obs column reports (instrument assumption 1, "
    "UNCONFIRMED). A figure that failed only because of that would look exactly like a derivation "
    "that is wrong, and this study cannot tell those apart",
    "whether the normal approximation holds for PERTURBED rows at |T| > 1.96. If perturbed rows carry "
    "far fewer cells than non-targeting rows, it is weaker there than this check can possibly show. "
    "The cell-count DISTRIBUTION over all perturbed rows is reported as a property of the file, as "
    "quantiles with NO perturbation identity attached, because a perturbed row's num_cells_filtered "
    "is outcome-adjacent: knocking down an essential gene lowers its cell count, which is a fitness "
    "phenotype and not metadata. PER-FACTOR COUNTS ARE NOT READ IN N1b, AT ALL -- not sampled, not "
    "spot-checked, not printed for a sanity check",
    "anything whatever about which genes respond to which factor. No factor row is read on a pass or a fail",
)


# --- the usable rows and the uncalibratable genes ----------------------------------------------------------


def _finite(v: Any) -> bool:
    return isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v)


def usable_rows(cell_counts: Sequence[float]) -> dict[str, Any]:
    """Which non-targeting rows carry a cell count T can be formed from, and why the others do not.

    No row is dropped by anything measured on it. The only test is the obs field the scale factor
    needs: finite and at least MIN_CELLS.
    """
    used, dropped = [], []
    for i, n in enumerate(cell_counts):
        if not _finite(n):
            dropped.append({"row": i, "reason": "num_cells_filtered is not finite"})
        elif n < MIN_CELLS:
            dropped.append({"row": i, "reason": f"num_cells_filtered {n} is below {MIN_CELLS}"})
        else:
            used.append(i)
    return {"rows_given": len(cell_counts), "rows_used": used, "rows_dropped": dropped}


# --- the strata: quintiles of control expression -----------------------------------------------------------


def expression_quintiles(expression: Sequence[float], strata: int = STRATA) -> list[int]:
    """Genes ranked by control expression (ties broken by position), cut into `strata` groups whose
    sizes differ by at most one. Stratum 0 holds the LEAST expressed genes, where the derivation is
    expected to fail: X is a mean of z-scores over few non-zero counts there, so the normal
    approximation (assumption 3) is weakest exactly where an unstratified pass would hide it."""
    n = len(expression)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda j: (expression[j], j))
    out = [0] * n
    for rank, j in enumerate(order):
        out[j] = rank * strata // n
    return out


# --- the two figures, and the intervals around them --------------------------------------------------------


def wilson_interval(k: int, n: int, z: float = 1.96) -> dict[str, Any]:
    """A Wilson score interval on a share, used where a cluster bootstrap is DEGENERATE.

    At 0 or n successes every resample of the clusters returns the identical value, so the bootstrap
    interval is a point and decides nothing. The route out is an interval that does not resample: a
    Wilson score interval on k of n CLUSTERS. n is the number of row-clusters and never the number of
    (row, gene) values, because values within one row are not independent, and using the larger
    denominator would shrink the interval by a factor the data does not earn.
    """
    if n <= 0:
        return {"kind": "wilson", "low": None, "high": None, "note": "no cluster to bound"}
    p = k / n
    d = 1.0 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return {
        "kind": "wilson",
        "clusters": n,
        "successes": k,
        "low": max(0.0, centre - half),
        "high": min(1.0, centre + half),
        "z": z,
        "note": "a bound on the CLUSTER share, not on the (row, gene) share; the bootstrap was degenerate",
    }


def _percentile(sorted_values: Sequence[float], q: float) -> float:
    n = len(sorted_values)
    if n == 1:
        return float(sorted_values[0])
    pos = q * (n - 1)
    lo = int(math.floor(pos))
    hi = min(lo + 1, n - 1)
    frac = pos - lo
    return float(sorted_values[lo] * (1 - frac) + sorted_values[hi] * frac)


def bootstrap_interval(
    draws: Sequence[float], point: float, *, degenerate_fallback: dict[str, Any] | None = None
) -> dict[str, Any]:
    """A 95% percentile interval over cluster resamples, with the identical-resample share beside it.

    `identical_resample_share` is the share of resamples whose statistic equals the point estimate
    exactly. It is reported for EVERY interval, degenerate or not, because it is the one number that
    says whether the interval came from variation in the data or from the absence of it. At 1.0 the
    interval is a point and the fallback bound is what a reader should use.
    """
    if not draws:
        return {"kind": "cluster_bootstrap", "low": None, "high": None, "identical_resample_share": None}
    same = sum(1 for d in draws if d == point)
    share = same / len(draws)
    s = sorted(float(d) for d in draws)
    out: dict[str, Any] = {
        "kind": "cluster_bootstrap",
        "resamples": len(draws),
        "low": _percentile(s, 0.025),
        "high": _percentile(s, 0.975),
        "identical_resample_share": share,
        "degenerate": share >= 1.0,
    }
    if share >= 1.0:
        out["decides"] = "nothing: every resample returned the identical value"
        out["fallback"] = degenerate_fallback
    return out


def ratio_report(observed: float | None, nominal: float, name: str) -> dict[str, Any]:
    """A ratio NEVER alone. Its absolute level and its control sit beside it, in one record.

    On 2026-10-02 a result in this project reported a 2.9-4.1x relative excess whose absolute level
    was 6.0%, so 94% of pairs disagreed and the ratio read as a finding the level did not support.
    Every ratio this module emits carries `absolute` and `control` in the same object so the two
    cannot be separated by a reader or by a later summary.
    """
    return {
        "name": name,
        "absolute": observed,
        "control": nominal,
        "control_kind": CONTROL_KIND,
        "ratio_to_control": (observed / nominal) if (observed is not None and nominal) else None,
    }


# --- the pass rule ------------------------------------------------------------------------------------------


def band_verdict(value: float | None, lo: float, hi: float) -> dict[str, Any]:
    if value is None:
        return {"in_band": False, "value": None, "band": [lo, hi], "reason": "the figure could not be formed"}
    return {
        "in_band": lo <= value <= hi,
        "value": value,
        "band": [lo, hi],
        "reason": None if lo <= value <= hi else "outside the registered band",
    }


def decide(overall: dict[str, Any], strata: Sequence[dict[str, Any]]) -> dict[str, Any]:
    """The registered rule applied, and nothing else. Returns the verdict and every leg that decided it."""
    legs: list[dict[str, Any]] = []
    for scope, rec in [("overall", overall), *[(f"quintile {s['stratum']}", s) for s in strata]]:
        legs.append(
            {
                "scope": scope,
                "variance": band_verdict(rec.get("var_primary"), VAR_LO, VAR_HI),
                "tail": band_verdict(rec.get("tail_primary"), TAIL_LO, TAIL_HI),
            }
        )
    failed = [
        {"scope": leg["scope"], "legs": [k for k in ("variance", "tail") if not leg[k]["in_band"]]}
        for leg in legs
        if not (leg["variance"]["in_band"] and leg["tail"]["in_band"])
    ]
    passed = not failed
    return {
        "pass_rule": PASS_RULE,
        "passed": passed,
        "pass_licence": PASS_LICENCE,
        "cannot_establish": list(CANNOT_ESTABLISH),
        "legs": legs,
        "failed_scopes": failed,
        "reading": (
            "T = X * sqrt(n) is calibrated on the non-targeting rows under the registered rule: both "
            "bands hold overall and in every control-expression quintile. The reading is about the NULL "
            "and about nothing else, and it licenses no more than this: " + PASS_LICENCE + "."
            if passed
            else "T = X * sqrt(n) is UNUSABLE. The registered rule is missed in "
            f"{len(failed)} of {len(legs)} scopes. The derivation is not repaired here and no second "
            "statistic is proposed; the fallbacks are Albert's decision and not this lane's."
        ),
        "strata_count": len(strata),
    }


# --- the remote file, read by exact byte ranges -------------------------------------------------------------

BLOCK = 1 << 16


class RangeFile(io.RawIOBase):
    """A remote file read in 64 KiB blocks by HTTP range requests, logging every range and its sha256.

    Used for the file's HEADERS only -- the superblock, the object headers and the small obs columns --
    so that h5py can report where `X` lives without any of `X` being decoded. The `X` rows themselves
    are fetched by `fetch_rows`, which asks for the exact bytes of a row and nothing else.
    """

    def __init__(self, url: str, size: int, agent: str = "GenomeOS/n1b-calibration"):
        self.url, self.size, self.pos, self.agent = url, size, 0, agent
        self.fetched = 0
        self.cache: dict[int, bytes] = {}
        self.log: list[dict[str, Any]] = []

    def readable(self) -> bool:
        return True

    def seekable(self) -> bool:
        return True

    def seek(self, off: int, whence: int = 0) -> int:
        self.pos = off if whence == 0 else self.pos + off if whence == 1 else self.size + off
        return self.pos

    def tell(self) -> int:
        return self.pos

    def _get(self, start: int, end: int, purpose: str) -> bytes:
        req = urllib.request.Request(
            self.url, headers={"Range": f"bytes={start}-{end}", "User-Agent": self.agent}
        )
        with urllib.request.urlopen(req, timeout=180) as r:  # noqa: S310
            if r.status != 206:
                raise RuntimeError(f"the server did not honour the range request (HTTP {r.status})")
            data = r.read()
        if len(data) != end - start + 1:
            raise RuntimeError(f"range {start}-{end} returned {len(data)} bytes, not {end - start + 1}")
        self.fetched += len(data)
        self.log.append(
            {
                "url": self.url,
                "range": f"bytes={start}-{end}",
                "bytes": len(data),
                "sha256": hashlib.sha256(data).hexdigest(),
                "purpose": purpose,
            }
        )
        return data

    def _block(self, i: int) -> bytes:
        if i not in self.cache:
            start, end = i * BLOCK, min(self.size, (i + 1) * BLOCK) - 1
            self.cache[i] = self._get(start, end, "header block")
        return self.cache[i]

    def readinto(self, b) -> int:
        n = min(len(b), self.size - self.pos)
        if n <= 0:
            return 0
        out, p = bytearray(), self.pos
        while len(out) < n:
            blk = self._block(p // BLOCK)
            off = p % BLOCK
            take = blk[off : off + n - len(out)]
            out += take
            p += len(take)
        b[:n] = out
        self.pos += n
        return n

    def fetch_rows(self, base: int, stride: int, rows: Iterable[int], purpose: str) -> dict[int, bytes]:
        """The exact bytes of each named row of a CONTIGUOUS dataset, one range request per row.

        One row at a time, and the bytes of one row are handed back and then dropped by the caller: no
        whole file and no whole matrix is opened. A row of this file is 32,992 bytes.
        """
        out = {}
        for r in rows:
            start = base + r * stride
            out[r] = self._get(start, start + stride - 1, f"{purpose} row {r}")
        return out


def x_layout(h: Any, name: str = "X") -> dict[str, Any]:
    """Where `X` lives, from the object header alone. Raises unless it is contiguous and float32."""
    ds = h[name]
    offset = ds.id.get_offset()
    if offset is None:
        raise RuntimeError(f"{name} is not contiguous: a chunked dataset has no single row range")
    if ds.dtype.itemsize != 4 or ds.dtype.kind != "f":
        raise RuntimeError(f"{name} is {ds.dtype}, not float32: the row stride would be wrong")
    rows, cols = ds.shape
    return {
        "dataset": name,
        "offset": int(offset),
        "rows": int(rows),
        "cols": int(cols),
        "dtype": str(ds.dtype),
        "row_stride_bytes": int(cols) * 4,
        "layout": "contiguous",
    }


def source_file(kind: str) -> dict[str, Any]:
    """The file identity, from the ONE copy in the project (N1's SOURCE). Nothing of N1 is altered."""
    return dict(SOURCE["files"][kind])


MD5_NOT_VERIFIABLE = (
    "the whole-file md5 Figshare publishes cannot be checked without fetching all 374,587,922 bytes, "
    "which is the cost this study exists to avoid. It is therefore NOT verified for any file read by "
    "range, and what stands in its place is the sha256 of every range read, logged with its byte offsets. "
    "A file whose bytes were swapped between two ranges would not be caught by this."
)


# --- the figures, from the T matrix -----------------------------------------------------------------------
#
# T is (usable rows x genes). At 514 rows and 8,248 genes that is 16.9 MB as float32, which is the whole
# of what this study holds in memory at once, beside two matrices of the same shape for the bootstrap.
# No file is opened whole: rows arrive one range request at a time and are written into this matrix.
#
# If a memory ceiling is wanted, INJECT the figure (`rss_ceiling_bytes`) and read the live reading only in
# a real run. A live reading taken inside a test process is valid only in a process that has done nothing
# substantial first, and a guard that reads a process-wide RSS figure measures the test process rather than
# the work -- which is why one such guard is being removed from this project.


def gene_figures(t: Any) -> dict[str, Any]:
    """Per gene, over the usable rows: the variance of T about the gene's own mean, the second moment
    about zero, and the count of |T| > 1.96. Nothing is summarised over genes here."""
    import numpy as np

    t = np.asarray(t, dtype=np.float64)
    m = t.shape[0]
    if m < 2:
        raise ValueError("a variance about the gene's own mean needs at least two rows")
    s1 = t.sum(axis=0)
    s2 = (t * t).sum(axis=0)
    return {
        "rows": m,
        "variance": (s2 - s1 * s1 / m) / (m - 1),
        "second_moment": s2 / m,
        "mean": s1 / m,
        "tail_count": (np.abs(t) > TAIL_Z).sum(axis=0),
    }


def summarise(figures: dict[str, Any], cols: Sequence[int], t: Any) -> dict[str, Any]:
    """The primary figure for one scope, and every secondary beside it. Only `var_primary` and
    `tail_primary` are read by `decide`."""
    import numpy as np

    cols = np.asarray(list(cols), dtype=np.int64)
    m = figures["rows"]
    if cols.size == 0:
        return {
            "genes": 0,
            "values": 0,
            "var_primary": None,
            "tail_primary": None,
            "reason": "no gene in this scope",
        }
    var = np.asarray(figures["variance"])[cols]
    second = np.asarray(figures["second_moment"])[cols]
    tail = np.asarray(figures["tail_count"])[cols]
    block = np.asarray(t, dtype=np.float64)[:, cols]
    values = int(m) * int(cols.size)
    grand = block.mean()
    pooled = float(((block - grand) ** 2).sum() / (values - 1))
    tail_share = float(tail.sum()) / values
    return {
        "genes": int(cols.size),
        "values": values,
        "rows": int(m),
        "var_primary": float(np.median(var)),
        "var_primary_summary": PRIMARY_VAR_SUMMARY,
        "var_mean_over_genes": float(var.mean()),
        "var_pooled_over_values": pooled,
        "second_moment_mean_over_genes": float(second.mean()),
        "t_mean_over_values": float(grand),
        "var_secondary_summaries": list(SECONDARY_VAR_SUMMARIES),
        "tail_primary": tail_share,
        "tail_primary_summary": PRIMARY_TAIL_SUMMARY,
        "tail_count": int(tail.sum()),
        "tail_median_over_genes": float(np.median(tail / m)),
        "tail_secondary_summaries": list(SECONDARY_TAIL_SUMMARIES),
        "rows_with_any_exceedance": int((np.abs(block) > TAIL_Z).any(axis=1).sum()),
        "ratios": [
            ratio_report(float(np.median(var)), NOMINAL_VAR, "variance of T, median over genes"),
            ratio_report(tail_share, NOMINAL_TAIL, f"share of |T| > {TAIL_Z}"),
        ],
    }


def cluster_bootstrap(
    t: Any, scopes: dict[str, Sequence[int]], points: dict[str, dict[str, Any]]
) -> dict[str, dict[str, Any]]:
    """A 95% percentile interval for every scope's two primary figures, resampling ROW-CLUSTERS.

    The row is the cluster: cells are pooled within a row and genes are correlated within a row, so a
    resample that drew (row, gene) values would treat 4.2 million correlated values as independent and
    return an interval the data has not earned. Rows resample with replacement, m of m.
    """
    import numpy as np

    t = np.asarray(t, dtype=np.float64)
    m, _ = t.shape
    t2 = t * t
    exceed = (np.abs(t) > TAIL_Z).astype(np.float64)
    names = list(scopes)
    colsets = {nm: np.asarray(list(scopes[nm]), dtype=np.int64) for nm in names}
    # per row, the exceedance count inside each scope: the tail share needs nothing per gene
    counts = np.stack([exceed[:, colsets[nm]].sum(axis=1) for nm in names], axis=1)  # (m, scopes)
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    draws: dict[str, dict[str, list[float]]] = {nm: {"var": [], "tail": []} for nm in names}
    done = 0
    while done < BOOTSTRAP_B:
        b = min(BOOTSTRAP_BLOCK, BOOTSTRAP_B - done)
        w = rng.multinomial(m, np.full(m, 1.0 / m), size=b).astype(np.float64)  # (b, m)
        s1 = w @ t
        s2 = w @ t2
        var_b = (s2 - s1 * s1 / m) / (m - 1)  # (b, genes)
        tail_b = w @ counts  # (b, scopes)
        for p, nm in enumerate(names):
            cs = colsets[nm]
            if cs.size == 0:
                continue
            draws[nm]["var"].extend(np.median(var_b[:, cs], axis=1).tolist())
            draws[nm]["tail"].extend((tail_b[:, p] / (m * cs.size)).tolist())
        done += b
    out: dict[str, dict[str, Any]] = {}
    for nm in names:
        pt = points[nm]
        rows_any = pt.get("rows_with_any_exceedance", 0)
        out[nm] = {
            "cluster": "non-targeting row",
            "clusters": int(m),
            "variance": bootstrap_interval(
                draws[nm]["var"],
                pt.get("var_primary"),
                degenerate_fallback={
                    "kind": "none available",
                    "note": "a variance has no exact binomial analogue; a degenerate interval here "
                    f"means {INCONCLUSIVE_WORDING} and the figure stands with no interval",
                },
            ),
            "tail": bootstrap_interval(
                draws[nm]["tail"],
                pt.get("tail_primary"),
                degenerate_fallback=wilson_interval(int(rows_any), int(m)),
            ),
            "assumptions": list(ASSUMPTIONS),
        }
    return out


ASSUMPTIONS = (
    "X averages exactly num_cells_filtered cells -- UNCONFIRMED by the producers' documentation; if the "
    "divisor differs, sqrt(n) is the wrong scale factor and every figure beside this line moves",
    "the cells pooled into a non-targeting row are independent",
    f"the normal approximation holds at |T| = {TAIL_Z}, the point the tail band reads",
    "the row is the resampling cluster; genes within a row are correlated and do not resample separately",
)


def t_from_rows(x: Any, cell_counts: Sequence[float]) -> Any:
    """T = X * sqrt(n), row by row. `x` is (rows x genes) of normalized X for the USABLE rows only and
    `cell_counts` their num_cells_filtered, in the same order. This is the whole derivation under test."""
    import numpy as np

    x = np.asarray(x, dtype=np.float64)
    n = np.asarray(list(cell_counts), dtype=np.float64)
    if x.shape[0] != n.shape[0]:
        raise ValueError(f"{x.shape[0]} rows of X against {n.shape[0]} cell counts")
    return x * np.sqrt(n)[:, None]


def uncalibratable_genes(t: Any, expression: Sequence[float]) -> list[int]:
    """Genes a band cannot be read on: any non-finite T, or a non-finite control expression. They leave
    before any figure is formed and are counted in the result."""
    import numpy as np

    t = np.asarray(t, dtype=np.float64)
    e = np.asarray(list(expression), dtype=np.float64)
    bad = ~np.isfinite(t).all(axis=0) | ~np.isfinite(e)
    return [int(j) for j in np.flatnonzero(bad)]
