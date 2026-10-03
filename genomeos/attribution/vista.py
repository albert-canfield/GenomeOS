# SPDX-License-Identifier: AGPL-3.0-or-later
"""The attribution held against measured enhancers: VISTA's in-vivo transgenic assays.

The VISTA Enhancer Browser (Lawrence Berkeley National Laboratory) has tested about 2,400
human sequences in e11.5 mouse embryos: each is positive in a named set of tissues (forebrain,
heart, limb, ...) or negative. It is the largest set of enhancers measured in a living animal,
and it is the ground truth the parser never had: the elements were chosen for conservation, so
constraint alone cannot separate the positives from the negatives, and whatever separates them
has learnt something about enhancers.

The loci table (120 kB, hg38 coordinates) is streamed once into data/knowledge/vista, local
and never committed. For every VISTA element on a chromosome this module records what GenomeOS
already says about that stretch of sequence, with no knowledge of the assay: does an ENCODE
enhancer-like element sit there; which CTCF node and which nearest coding TSS; how constrained
it is (Zoonomia phyloP); and, when a scorer is given, what AlphaGenome predicts on deleting the
whole element, gene and tissue. The summary compares positives with negatives on each of those,
and, for the positives, asks whether the predicted tissue falls in the measured one (neural,
heart, limb, ...). Everything measured is `experimental` evidence; the predictions stay
`predicted`.
"""

from __future__ import annotations

import bisect
import gzip
import json
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any

VISTA_LOCUS_URL = "https://gitlab.com/egsb-mfgl/vista-data/-/raw/main/locus.tsv.gz"
KNOWLEDGE = Path("data/knowledge/vista")
EVIDENCE = "experimental: VISTA Enhancer Browser, transgenic mouse e11.5 (LBNL); hg38 coordinates"
MAX_DELETION = 5_000  # longer elements are not deleted in the model
CONSTRAINED = 0.2  # constrained fraction that counts as a constrained element (as in the budget)

TISSUES = {
    "ba": "branchial arch",
    "hb": "hindbrain",
    "eye": "eye",
    "fb": "forebrain",
    "nt": "neural tube",
    "mb": "midbrain",
    "nose": "nose",
    "lb": "limb",
    "cn": "cranial nerve",
    "drg": "dorsal root ganglion",
    "ht": "heart",
    "tri": "trigeminal ganglion",
    "som": "somite",
    "gen": "genital tubercle",
    "mel": "melanocytes",
    "tail": "tail",
    "fm": "facial mesenchyme",
    "other": "other",
    "ear": "ear",
    "lv": "liver",
    "bv": "blood vessels",
    "pan": "pancreas",
}
GROUP_OF = {
    **dict.fromkeys(("fb", "mb", "hb", "nt", "cn", "drg", "tri", "eye"), "neural"),
    **dict.fromkeys(("ht", "bv"), "heart"),
    **dict.fromkeys(("lb", "som", "tail", "fm", "ba"), "limb and mesenchyme"),
    "lv": "liver",
    "pan": "pancreas",
}
# words in an AlphaGenome RNA-seq track name that place it in a VISTA tissue group
TRACK_WORDS = {
    "neural": (
        "brain",
        "cerebell",
        "cortex",
        "lobe",
        "hippocamp",
        "amygdala",
        "caudate",
        "putamen",
        "substantia",
        "hypothalam",
        "spinal cord",
        "accumbens",
        "neur",
        "purkinje",
        "astrocyte",
        "glia",
        "retina",
    ),
    "heart": ("heart", "ventricle", "atrium", "cardiac", "aorta", "coronary"),
    "limb and mesenchyme": ("skeletal muscle", "muscle", "chondrocyte", "osteoblast", "fibroblast"),
    "liver": ("liver", "hepat"),
    "pancreas": ("pancrea",),
}


@dataclass(frozen=True, slots=True)
class VistaElement:
    id: str
    chrom: str
    start: int  # 0-based half-open
    end: int
    status: str  # positive | negative
    tissues: tuple[str, ...]
    experiments: int

    @property
    def groups(self) -> set[str]:
        return tissue_groups(self.tissues)


def tissue_groups(tissues) -> set[str]:
    return {GROUP_OF[t] for t in tissues if t in GROUP_OF}


def group_of_track(tissue: str | None) -> str | None:
    """The VISTA tissue group an AlphaGenome track name belongs to, or None."""
    if not tissue:
        return None
    t = tissue.lower()
    for group, words in TRACK_WORDS.items():
        if any(w in t for w in words):
            return group
    return None


def locus_path(knowledge: Path = KNOWLEDGE) -> Path:
    return knowledge / "locus.tsv.gz"


def fetch_loci(knowledge: Path = KNOWLEDGE, url: str = VISTA_LOCUS_URL) -> Path:
    """Stream the loci table once (about 120 kB)."""
    dest = locus_path(knowledge)
    if dest.exists():
        return dest
    knowledge.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": "GenomeOS/0.1 (stream)"})
    with urllib.request.urlopen(req, timeout=300) as resp:  # noqa: S310
        dest.write_bytes(resp.read())
    return dest


def parse_loci(fh, chrom: str | None = None) -> list[VistaElement]:
    """Human elements (hs*) with hg38 coordinates from the loci table."""
    header = None
    out = []
    for line in fh:
        f = [x.strip('"') for x in line.rstrip("\n").split("\t")]
        if header is None:
            header = {k: i for i, k in enumerate(f)}
            continue
        vid = f[header["vista_id"]]
        coord = f[header["coordinate_hg38"]]
        if not vid.startswith("hs") or ":" not in coord:
            continue
        c, _, rng = coord.partition(":")
        if chrom and c != chrom:
            continue
        a, _, b = rng.partition("-")
        tissues = tuple(t for t in f[header["tissue"]].split(";") if t)
        out.append(
            VistaElement(
                vid,
                c,
                int(a) - 1,
                int(b),
                f[header["curation_status"]],
                tissues,
                int(f[header["experiments"]] or 0),
            )
        )
    out.sort(key=lambda e: (e.chrom, e.start))
    return out


def load_loci(chrom: str | None = None, knowledge: Path = KNOWLEDGE) -> list[VistaElement]:
    p = fetch_loci(knowledge)
    with gzip.open(p, "rt") as fh:
        return parse_loci(fh, chrom)


def nearest_coding_tss(ctx, start: int, end: int) -> dict[str, Any] | None:
    """The nearest coding TSS inside the element's CTCF node (the node inference for a bare region)."""
    from genomeos.genome.regulation import _tss

    mid = (start + end) // 2
    starts = [d.start for d in ctx.domains]
    i = bisect.bisect_right(starts, mid) - 1
    dom = (
        ctx.domains[i]
        if 0 <= i < len(ctx.domains) and ctx.domains[i].start <= mid < ctx.domains[i].end
        else None
    )
    if dom is None:
        return None
    best = None
    for g in ctx.annotation.genes.values():
        if g.type != "protein_coding" or g.locus.chrom != ctx.chrom:
            continue
        t = _tss(g)
        if dom.start <= t < dom.end:
            d = abs(t - mid)
            if best is None or d < best[1]:
                best = (g.symbol, d)
    if best is None:
        return {"domain": dom.id, "gene": None, "distance": None}
    return {"domain": dom.id, "gene": best[0], "distance": best[1], "basis": "nearest TSS in domain"}


def annotate(ctx, elements: list[VistaElement], stats=None) -> list[dict[str, Any]]:
    """What the registry, the node model and constraint say about each VISTA element, blind to the assay."""
    ccres = sorted(ctx.ccres, key=lambda c: c.start)
    ccre_starts = [c.start for c in ccres]
    rows = []
    for k, e in enumerate(elements):
        i = bisect.bisect_left(ccre_starts, e.start - 5_000)
        over = []
        while i < len(ccres) and ccres[i].start < e.end:
            c = ccres[i]
            if c.end > e.start:
                over.append(c)
            i += 1
        classes = sorted({c.cls for c in over})
        node = nearest_coding_tss(ctx, e.start, e.end)
        reg = next((ctx.element_by_id[c.id] for c in over if c.cls in ("dELS", "pELS")), None)
        s = stats[k] if stats else None
        rows.append(
            {
                "id": e.id,
                "start": e.start,
                "end": e.end,
                "length": e.end - e.start,
                "status": e.status,
                "tissues": list(e.tissues),
                "groups": sorted(e.groups),
                "ccres": [c.id for c in over],
                "ccre_classes": classes,
                "enhancer_like": any(c in ("dELS", "pELS") for c in classes),
                "any_ccre": bool(over),
                "domain": node["domain"] if node else "",
                "inferred": (
                    reg.targets[0]
                    if reg is not None and reg.targets
                    else (node if node and node.get("gene") else None)
                ),
                "constrained_fraction": round(s.fraction_above or 0.0, 3)
                if s is not None and s.bases
                else None,
                "phylop_mean": round(s.mean, 3) if s is not None and s.bases else None,
            }
        )
    return rows


def score(ctx, scorer, row: dict[str, Any], cache: Path | None = None) -> dict[str, Any]:
    """Delete the whole VISTA element in AlphaGenome; predicted target and tissue join the row."""
    from genomeos.predict.enhancer_target import CACHE, score_element

    if row["length"] > MAX_DELETION:
        row["predicted"] = None
        row["skipped"] = f"longer than {MAX_DELETION} bp"
        return row
    dom = ctx.domain_by_id.get(row["domain"])
    inferred = [row["inferred"]] if row.get("inferred") and row["inferred"].get("gene") else None
    r = score_element(
        scorer,
        ctx.genome.fetch,
        ctx.chrom,
        row["id"],
        row["start"],
        row["end"],
        inferred,
        list(dom.genes) if dom is not None else [],
        coding=ctx.coding,
        cache=cache or CACHE,
    )
    row["predicted"] = r.get("predicted")
    row["predicted_coding"] = r.get("predicted_coding")
    row["verdict_coding"] = r.get("verdict_coding")
    p = row["predicted"]
    row["predicted_group"] = group_of_track(p["tissue"]) if p else None
    row["tissue_agrees"] = (
        (row["predicted_group"] in row["groups"]) if p and row["predicted_group"] and row["groups"] else None
    )
    return row


#: THE REGISTERED FALSIFIER, word for word from 24adf33: "a VISTA row whose `tissue_agrees` is False
#: on the head and True for some other gene at the bar. One such row falsifies the stored agreement
#: rate as an agreement rate of the ELEMENT; zero of them leaves the stored rate standing unchanged."
WINDOW_FALSIFIER = (
    "a row whose tissue_agrees is False on the head and True for some other gene at the bar; one "
    "such row falsifies the stored rate as a rate of the ELEMENT, zero leave it standing unchanged"
)
#: What the window CANNOT answer here, stated at the move. `tissue_agrees` is None for a row with no
#: observed groups and for a row with no predicted group, and those two silences are different; the
#: window can turn a False into a True and can give a group to a row that had none, but it cannot
#: supply observed groups, so a row with no `groups` stays unjudged under the window as well.
WINDOW_CANNOT_SUPPLY_OBSERVED_GROUPS = (
    "a VISTA row with no observed expression groups is unjudged with a window as without one: the "
    "window supplies a predicted group per gene and never an observed one"
)
#: An over-length element is SKIPPED by `score` and its head is written as None (the registered read
#: at the skip). The window cannot rescue it either: nothing was ever scored, so there is no cached
#: record, and such rows are counted as `skipped` and never as "no gene at the bar".
SKIPPED_IS_NOT_UNPREDICTED = (
    "an element longer than MAX_DELETION was never scored, so it has no window; it is counted as "
    "skipped and never as an element whose window named nothing"
)


def window_tissue_agreement(
    rows: list[dict[str, Any]], chrom: str, responses=None, min_effect: float = 0.1
) -> dict[str, Any]:
    """vista.py's own falsifier: does any row agree on a gene at the bar that its head disagreed on.

    `{}` without a reader, so `summarise` contributes no new key and every committed figure of this
    module is produced by the unchanged code path. `tissue_agrees` is NOT recomputed and NOT
    overwritten; the window's answer is reported beside it.

    The group of a gene at the bar is `group_of_track` applied to that gene's own strongest track,
    which `attribution/targets.tracks_at_bar` returns -- the amendment recorded in
    `onetarget2.AMENDMENTS`. A gene whose track has no VISTA group contributes nothing, exactly as a
    head with no group does today.
    """
    if responses is None:
        return {}
    from genomeos.attribution import onetarget2 as ot
    from genomeos.attribution.targets import genes_at_bar, tracks_at_bar

    flipped, not_cached, skipped, unjudged, disagreements, judged = [], 0, 0, 0, 0, 0
    for r in rows:
        if r.get("skipped"):
            skipped += 1
            continue
        groups = r.get("groups")
        if groups is None and r.get("tissues") is not None:
            #: A COMMITTED vista_chr*.json row carries `tissues` and not `groups`, because `groups`
            #: is a property of VistaElement and was never serialised. It is re-derived here with
            #: the module's OWN rule, `tissue_groups`, rather than treating a stored row as
            #: unjudgeable -- which is what the first run of this arm did, reporting 2,179 of 2,223
            #: rows as having no observed groups while 91 of them carried a tissue_agrees of False.
            groups = sorted(tissue_groups(r["tissues"]))
        if not groups:
            unjudged += 1
            continue
        record = responses.element(chrom, str(r["id"]))
        bar = tracks_at_bar(record, min_effect)
        if bar is None:
            not_cached += 1
            continue
        plain = genes_at_bar(record, min_effect) or []
        if [(g, v) for g, v, _ in bar] != list(plain):
            raise ValueError(
                f"{r['id']}: tracks_at_bar and genes_at_bar disagree, so the amendment has become a "
                "second size rule. Nothing is reported."
            )
        head = ((r.get("predicted") or {}).get("gene")) or None
        if bar and (bar[0][0] != head):
            disagreements += 1
        judged += 1
        if r.get("tissue_agrees") is not False:
            continue
        for gene, signed, tissue in bar:
            if gene == head:
                continue
            group = group_of_track(tissue)
            if group and group in groups:
                flipped.append(
                    {
                        "id": r["id"],
                        "status": r.get("status"),
                        "head": head,
                        "head_group": r.get("predicted_group"),
                        "observed_groups": sorted(groups),
                        "gene_at_the_bar": gene,
                        "its_group": group,
                        "its_log2_fold_change": signed,
                    }
                )
                break
    ot.check_head_invariant({"head_disagreements": disagreements}, f"attribution/vista.py {chrom}")
    return {
        "window": {
            "falsifier": WINDOW_FALSIFIER,
            "cannot_supply_observed_groups": WINDOW_CANNOT_SUPPLY_OBSERVED_GROUPS,
            "skipped_is_not_unpredicted": SKIPPED_IS_NOT_UNPREDICTED,
            "rows": len(rows),
            "rows_judged_under_the_window": judged,
            "rows_skipped_as_over_length": skipped,
            "rows_with_no_observed_groups": unjudged,
            "rows_not_cached": not_cached,
            "rows_agreeing_on_a_gene_the_head_dropped": len(flipped),
            "flipped": flipped[:25],
            "verdict_against_its_own_falsifier": (
                f"FIRES: {len(flipped)} row(s) disagree on the head and agree on another gene at the "
                "bar, so the stored rate is a rate of the HEAD and not of the element"
                if flipped
                else "does NOT fire: no row agrees on a gene other than its head, so the stored rate "
                "stands unchanged as a rate of the element"
            ),
        }
    }


def _mean(xs: list[float]) -> float | None:
    xs = [x for x in xs if x is not None]
    return round(sum(xs) / len(xs), 3) if xs else None


def _expected_agreement(judged: list[dict[str, Any]], pool: list[dict[str, Any]]) -> float | None:
    """Agreement if the measured tissues were shuffled among the elements: the chance level."""
    labelled = [set(q["groups"]) for q in pool if q["groups"]]
    if not judged or len(labelled) < 2:
        return None
    total = 0.0
    for r in judged:
        others = [g for g in labelled if g != set(r["groups"])] or labelled
        total += sum(r["predicted_group"] in g for g in others) / len(others)
    return round(total / len(judged), 3)


def summarise(
    rows: list[dict[str, Any]], chrom: str | None = None, responses=None, min_effect: float = 0.1
) -> dict[str, Any]:
    """Positives against negatives on what GenomeOS says blind: registry, constraint, prediction, tissue.

    With `responses=None`, the default and what every committed vista result was produced under,
    nothing is opened and the result is exactly what it was. With a reader (and the chromosome, which
    the rows do not carry) one key is added, `window`, and no stored field moves: see
    `window_tissue_agreement`.
    """
    out: dict[str, Any] = {}
    for status in ("positive", "negative"):
        rs = [r for r in rows if r["status"] == status]
        n = len(rs)
        scored = [r for r in rs if "predicted" in r and not r.get("skipped")]
        named = [r for r in scored if r.get("predicted")]
        judged = [r for r in named if r.get("tissue_agrees") is not None]
        out[status] = {
            "elements": n,
            "enhancer_like": round(sum(r["enhancer_like"] for r in rs) / n, 3) if n else None,
            "any_ccre": round(sum(r["any_ccre"] for r in rs) / n, 3) if n else None,
            "in_a_node": round(sum(bool(r["domain"]) for r in rs) / n, 3) if n else None,
            "mean_constrained_fraction": _mean([r["constrained_fraction"] for r in rs]),
            "constrained": (
                round(sum((r["constrained_fraction"] or 0) >= CONSTRAINED for r in rs) / n, 3) if n else None
            ),
            "scored": len(scored),
            "with_predicted_target": len(named),
            "fraction_with_target": round(len(named) / len(scored), 3) if scored else None,
            "strong": sum(1 for r in named if r["predicted"]["strength"] == "strong"),
            "coding_target_agrees_with_nearest": (
                round(
                    sum(r.get("verdict_coding") == "agrees with nearest TSS in domain" for r in named)
                    / len(named),
                    3,
                )
                if named
                else None
            ),
            "tissue_judged": len(judged),
            "tissue_agrees": round(sum(r["tissue_agrees"] for r in judged) / len(judged), 3)
            if judged
            else None,
            "tissue_agrees_expected": _expected_agreement(judged, rs),
        }
    pos, neg = out["positive"], out["negative"]
    out["separation"] = {
        k: (round(pos[k] - neg[k], 3) if pos.get(k) is not None and neg.get(k) is not None else None)
        for k in ("enhancer_like", "constrained", "fraction_with_target", "mean_constrained_fraction")
    }
    groups: dict[str, int] = {}
    for r in rows:
        if r["status"] == "positive":
            for g in r["groups"]:
                groups[g] = groups.get(g, 0) + 1
    out["positive_groups"] = dict(sorted(groups.items(), key=lambda kv: -kv[1]))
    out["evidence"] = EVIDENCE
    if responses is not None:
        if not chrom:
            raise ValueError(
                "summarise needs the chromosome with a reader: the rows do not carry it and the "
                "response cache is keyed by chromosome, so a window could be read from the wrong one"
            )
        out.update(window_tissue_agreement(rows, chrom, responses, min_effect))
    return out


def save_knowledge(chrom: str, rows: list[dict[str, Any]], knowledge: Path = KNOWLEDGE) -> Path:
    p = knowledge / f"rows_{chrom}.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(rows, indent=1))
    return p


__all__ = [
    "VistaElement",
    "annotate",
    "fetch_loci",
    "group_of_track",
    "load_loci",
    "parse_loci",
    "score",
    "summarise",
    "tissue_groups",
]
