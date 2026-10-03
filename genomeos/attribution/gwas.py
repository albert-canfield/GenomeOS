# SPDX-License-Identifier: AGPL-3.0-or-later
"""Human trait associations (the GWAS Catalog) inside the attributed elements: does a consequence land here?

An eQTL says an element moves a gene; a reporter says it drives transcription; VISTA says it is an
enhancer in an embryo. None says the element matters to a person. The GWAS Catalog's lead variants
do, weakly and statistically: a lead variant is the marker most associated with a trait in a locus,
and it sits inside the causal element only sometimes (linkage carries the signal a kilobase or more).
So the question here is one of enrichment: do the elements the attribution singles out (a named
target, a strong effect, constraint, a VISTA positive) hold a lead variant more often than the rest,
and which traits? The catalog's association table (600 MB unpacked, 1.2 million rows) is streamed
from its zip once; only the rows inside or within MARGIN of an element are kept under
data/knowledge/gwas, local, with the rows on the same elements shifted SHIFT bases along the
chromosome, which is the chance level the enrichment is read against. The catalog's mapped gene is its
own nearest-gene call, recorded next to the predicted and inferred targets for the reader's eye, not as a
measurement. Evidence: `curated`.
"""

from __future__ import annotations

import io
import json
import time
import urllib.request
import zipfile
from pathlib import Path
from typing import Any

from genomeos.attribution.eqtl import Intervals

GWAS_ZIP = "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations-full.zip"
KNOWLEDGE = Path("data/knowledge/gwas")
MARGIN = 1_000  # a lead variant this close counts as landing on the element
SHIFT = 100_000  # the control: the same elements moved this far along the chromosome
SHIFTED = "|shifted"
EVIDENCE = "curated: NHGRI-EBI GWAS Catalog lead variants (GRCh38 positions), all associations"
HIT_COLUMNS = ("chrom", "pos", "rs", "trait", "mapped_gene", "context", "pvalue_mlog", "pubmed", "elements")


def parse_rows(fh):
    """(chrom, pos, rs, trait, mapped_gene, context, -log10 p, pubmed) for rows with one GRCh38 position."""
    header = None
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if header is None:
            header = {k: i for i, k in enumerate(f)}
            continue
        try:
            chrom, pos = f[header["CHR_ID"]], f[header["CHR_POS"]]
        except IndexError:
            continue
        if not chrom or not pos or ";" in pos or "x" in pos or not pos.isdigit():
            continue
        yield (
            f"chr{chrom}" if not chrom.startswith("chr") else chrom,
            int(pos),
            f[header["SNPS"]],
            f[header["DISEASE/TRAIT"]],
            f[header["MAPPED_GENE"]],
            f[header["CONTEXT"]],
            f[header["PVALUE_MLOG"]],
            f[header["PUBMEDID"]],
        )


def distil(
    intervals: Intervals,
    knowledge: Path = KNOWLEDGE,
    source: str | Path = GWAS_ZIP,
    progress=None,
) -> dict[str, Any]:
    """Stream the catalog once; keep the associations that land on an element (± MARGIN)."""
    knowledge.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    if isinstance(source, str) and source.startswith(("http://", "https://")):
        req = urllib.request.Request(source, headers={"User-Agent": "GenomeOS/0.1 (stream)"})
        with urllib.request.urlopen(req, timeout=900) as resp:  # noqa: S310
            blob = io.BytesIO(resp.read())
    else:
        blob = io.BytesIO(Path(source).read_bytes())
    scanned = kept = 0
    dest = knowledge / "hits.tsv"
    with zipfile.ZipFile(blob) as zf, dest.open("w") as out:
        out.write("\t".join(HIT_COLUMNS) + "\n")
        name = next(n for n in zf.namelist() if n.endswith(".tsv"))
        with zf.open(name) as raw, io.TextIOWrapper(raw, encoding="utf-8", errors="replace") as fh:
            for chrom, pos, rs, trait, gene, ctx, mlog, pmid in parse_rows(fh):
                scanned += 1
                if progress and scanned % 200_000 == 0:
                    progress(f"GWAS Catalog: {scanned:,} associations read, {kept:,} on elements")
                ids = intervals.at(chrom, pos)
                if ids:
                    kept += 1
                    out.write(
                        f"{chrom}\t{pos}\t{rs}\t{trait}\t{gene}\t{ctx}\t{mlog}\t{pmid}\t{','.join(ids)}\n"
                    )
    summary = {
        "source": str(source),
        "associations_read": scanned,
        "hits": kept,
        "seconds": round(time.time() - t0),
        "elements_indexed": intervals.n,
        "margin_bp": MARGIN,
        "evidence": EVIDENCE,
    }
    (knowledge / "distil_summary.json").write_text(json.dumps(summary, indent=1))
    return summary


def load_hits(knowledge: Path = KNOWLEDGE) -> dict[str, list[dict[str, Any]]]:
    out: dict[str, list[dict[str, Any]]] = {}
    p = knowledge / "hits.tsv"
    if not p.exists():
        return out
    with p.open() as fh:
        next(fh, None)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            row = {
                "chrom": f[0],
                "pos": int(f[1]),
                "rs": f[2],
                "trait": f[3],
                "mapped_gene": f[4],
                "context": f[5],
                "pvalue_mlog": float(f[6]) if f[6] else None,
                "pubmed": f[7],
            }
            for eid in f[8].split(","):
                out.setdefault(eid, []).append(row)
    return out


def index(elements: list[dict[str, Any]], margin: int = MARGIN, shift: int = SHIFT) -> Intervals:
    """Every element with its margin, and a copy shifted along the chromosome as the chance control."""
    iv = Intervals()
    for e in elements:
        iv.add(e["chrom"], max(0, e["start"] - margin), e["end"] + margin, e["key"])
        iv.add(e["chrom"], max(0, e["start"] - margin + shift), e["end"] + margin + shift, e["key"] + SHIFTED)
    return iv.freeze()


def _frac(rows, pred) -> float | None:
    return round(sum(1 for r in rows if pred(r)) / len(rows), 3) if rows else None


#: THE ONE MODULE IN THE CENSUS WHOSE HEADLINE FIGURE CAN ONLY MOVE ONE WAY, which is why its
#: registered falsifier is "the agreement share not rising": an element agrees if ANY gene at the bar
#: is the catalogue's mapped gene, and agreeing on the HEAD is a special case of that, so the window
#: set is a superset of the head set. A share that falls, or an element that loses its agreement, is
#: a DEFECT IN THE MOVE and not a finding, and `window_agreement` refuses rather than reporting it.
WINDOW_CAN_ONLY_RISE = (
    "agreement on the head is a special case of agreement on any gene at the bar, so the window set "
    "contains the head set; an element that loses its agreement means the window was not formed from "
    "the same element and is a defect, not a result"
)
#: This module reads the ANY-GENE head (`predicted`), not the coding head, so the 2026-09-27
#: coding-head exception does not reach these figures; the any-gene invariant does, and
#: `onetarget2.check_head_invariant` is called on this population.
WINDOW_READS_THE_ANY_GENE_HEAD = (
    "gwas.py reads the any-gene head field and not the coding one, so the registered coding-head "
    "limit is not load-bearing here and the any-gene invariant is checked on this population"
)


class WindowLostAnAgreementError(RuntimeError):
    """An element agreed on its head and not on its window. That cannot happen; the move is wrong."""


def _mapped_genes(h: dict[str, Any]) -> list[str]:
    """The catalogue's mapped gene(s) for one hit, by the same rule `summarise` has always used.

    A new helper rather than a rewrite of the two expressions above it: those lines produce committed
    figures and are left exactly as they are.
    """
    return (h.get("mapped_gene") or "").replace(" - ", ",").split(",")


def window_agreement(
    with_hit: list[dict[str, Any]],
    responses=None,
    coding: set[str] | None = None,
    min_effect: float = 0.1,
) -> dict[str, Any]:
    """Agreement with the catalogue read off the WHOLE window instead of the one predicted gene.

    `{}` without a reader, so `summarise` contributes no new key and every committed figure of this
    module is produced by the unchanged code path.

    The refusal is the point. `WINDOW_CAN_ONLY_RISE`: an element whose head is the mapped gene must
    also have that gene at its bar, because the head is always at its own bar. An element that
    agreed on the head and does not agree on the window is therefore a defect in how the window was
    formed -- the wrong element, the wrong chromosome, a stale cache -- and this function raises and
    names the element rather than publishing a share that fell.
    """
    if responses is None:
        return {}
    from genomeos.attribution import onetarget2 as ot
    from genomeos.attribution.targets import window_reading

    head_agree, bar_agree, gained, not_cached, disagreements = [], [], [], 0, 0
    for e in with_hit:
        mapped = {g for h in e["gwas"] for g in _mapped_genes(h) if g}
        on_head = ((e.get("predicted") or {}).get("gene") or "") in mapped
        w = window_reading(responses, e["chrom"], e, coding=coding, min_effect=min_effect)
        if w.not_cached:
            #: AN ELEMENT THE CACHE DOES NOT HOLD KEEPS ITS COMPACT ANSWER, which is
            #: `WindowReading`'s own contract: not_cached is neither an empty window nor a zero, and
            #: a consumer must keep its compact reading and count it by name. Doing anything else
            #: here would let the share FALL on elements that were never scored rather than on
            #: anything the window said -- the first form of this function raised on exactly such an
            #: element, and the refusal was right about the shape and wrong about the cause.
            not_cached += 1
            on_bar = on_head
        else:
            if w.head_agrees is False:
                disagreements += 1
            on_bar = bool(mapped & set(w.genes))
        if on_head:
            head_agree.append(e)
            if not on_bar:
                raise WindowLostAnAgreementError(
                    f"{e.get('key') or e.get('id')}: the catalogue's mapped gene is this element's "
                    f"predicted head and is NOT at its bar. {WINDOW_CAN_ONLY_RISE}. The share is not "
                    "published; the move is wrong."
                )
        if on_bar:
            bar_agree.append(e)
            if not on_head:
                gained.append(
                    {
                        "key": e.get("key") or e.get("id"),
                        "head": (e.get("predicted") or {}).get("gene"),
                        "mapped": sorted(mapped),
                        "rank_in_the_window": min(
                            (w.genes.index(g) + 1 for g in mapped if g in w.genes), default=None
                        ),
                    }
                )
    ot.check_head_invariant({"head_disagreements": disagreements}, "attribution/gwas.py window arm")
    n = len(with_hit)
    return {
        "catalog_gene_is_at_the_bar": round(len(bar_agree) / n, 3) if n else None,
        "window": {
            "falsifier": "the agreement share not rising",
            "can_only_rise": WINDOW_CAN_ONLY_RISE,
            "limit": WINDOW_READS_THE_ANY_GENE_HEAD,
            "elements_with_a_lead_variant": n,
            "agree_on_the_head": len(head_agree),
            "agree_anywhere_in_the_window": len(bar_agree),
            "elements_gaining_an_agreement": len(gained),
            "not_cached": not_cached,
            "gained": gained[:25],
            "verdict_against_its_own_falsifier": (
                f"FIRES: {len(gained)} element(s) agree with the catalogue on a gene the one-target "
                "projection dropped, so the stored share understates agreement of the ELEMENT"
                if gained
                else "does NOT fire: no element agrees on a gene other than its head, so the stored "
                "share stands unchanged as a share of the element and not only of the head"
            ),
        },
    }


def summarise(
    elements: list[dict[str, Any]],
    hits: dict[str, list[dict[str, Any]]],
    responses=None,
    coding: set[str] | None = None,
) -> dict[str, Any]:
    """Enrichment of lead variants over the attribution's own partitions of the elements.

    With `responses=None`, the default and what every committed figure was produced under, nothing
    is opened and the result is exactly what it was: no `catalog_gene_is_at_the_bar` key and no
    `window` key. With a reader, those two are added and nothing else moves -- see
    `window_agreement`, which refuses if any element loses an agreement.
    """
    for e in elements:
        e["gwas"] = hits.get(e["key"], [])
        e["gwas_shifted"] = len(hits.get(e["key"] + SHIFTED, []))
    has = lambda e: bool(e["gwas"])  # noqa: E731
    shifted = _frac(elements, lambda e: e["gwas_shifted"] > 0)
    named = [e for e in elements if (e.get("predicted") or {}).get("gene")]
    unnamed = [e for e in elements if not (e.get("predicted") or {}).get("gene") and "predicted" in e]
    strong = [e for e in named if e["predicted"].get("strength") == "strong"]
    cons = [e for e in elements if (e.get("constrained_fraction") or 0) >= 0.2]
    uncons = [
        e for e in elements if e.get("constrained_fraction") is not None and e["constrained_fraction"] < 0.2
    ]
    with_hit = [e for e in elements if e["gwas"]]
    agree_pred = [
        e
        for e in with_hit
        if (e.get("predicted") or {}).get("gene")
        and any(
            (e["predicted"]["gene"] in (h["mapped_gene"] or "").replace(" - ", ",").split(","))
            for h in e["gwas"]
        )
    ]
    agree_inf = [
        e
        for e in with_hit
        if (e.get("inferred") or {}).get("gene")
        and any(
            (e["inferred"]["gene"] in (h["mapped_gene"] or "").replace(" - ", ",").split(","))
            for h in e["gwas"]
        )
    ]
    traits: dict[str, int] = {}
    for e in with_hit:
        for t in {h["trait"] for h in e["gwas"]}:
            traits[t] = traits.get(t, 0) + 1
    out = {
        "elements": len(elements),
        "with_lead_variant": len(with_hit),
        "fraction_with_lead_variant": _frac(elements, has),
        "shifted_control_fraction": shifted,
        "enrichment_over_shifted": (
            round(_frac(elements, has) / shifted, 2) if shifted and _frac(elements, has) is not None else None
        ),
        "associations": sum(len(e["gwas"]) for e in elements),
        "by_partition": {
            "named_target": {"elements": len(named), "with_lead_variant": _frac(named, has)},
            "no_target": {"elements": len(unnamed), "with_lead_variant": _frac(unnamed, has)},
            "strong_effect": {"elements": len(strong), "with_lead_variant": _frac(strong, has)},
            "constrained": {"elements": len(cons), "with_lead_variant": _frac(cons, has)},
            "not_constrained": {"elements": len(uncons), "with_lead_variant": _frac(uncons, has)},
        },
        "catalog_gene_is_the_predicted_target": round(len(agree_pred) / len(with_hit), 3)
        if with_hit
        else None,
        "catalog_gene_is_the_inferred_target": round(len(agree_inf) / len(with_hit), 3) if with_hit else None,
        **window_agreement(with_hit, responses, coding),
        "top_traits": dict(sorted(traits.items(), key=lambda kv: -kv[1])[:15]),
        "margin_bp": MARGIN,
        "evidence": EVIDENCE,
    }
    if any("status" in e for e in elements):
        pos = [e for e in elements if e.get("status") == "positive"]
        neg = [e for e in elements if e.get("status") == "negative"]
        out["by_partition"]["vista_positive"] = {"elements": len(pos), "with_lead_variant": _frac(pos, has)}
        out["by_partition"]["vista_negative"] = {"elements": len(neg), "with_lead_variant": _frac(neg, has)}
    return out
