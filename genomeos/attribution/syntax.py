# SPDX-License-Identifier: AGPL-3.0-or-later
"""Syntax against values: which bases of a gene never change, and which ones people differ at.

Albert's framing (2026-09-12): a language has syntax, operators and values. In a
genome the syntax is what selection has held still across mammals, and a value is
a position where people differ. Read together at one gene, over the reference and
every genome imported locally, each base falls in one of three classes:

- syntax: constrained across 241 mammals (Zoonomia phyloP at or above 2.27) and the
  same in everyone we hold;
- value: unconstrained, and people differ there, the free slots;
- value in syntax: constrained *and* variable between people, the rare positions
  where a change lands on something selection has protected. Eye colour's
  rs12913832, inside the HERC2 enhancer that controls OCA2, is the textbook one.

Nothing here decides what a value does; it says where the values are and which of
them sit on syntax, so the next question (what does this slot control) is asked of
the right positions. Variants are read from the persons' own files, never copied.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any

from genomeos.attribution.bigwig import BigWig
from genomeos.attribution.constraint import PHYLOP_241_URL, PHYLOP_THRESHOLD
from genomeos.coords import Strand
from genomeos.genome.annotation import Annotation, default_gencode
from genomeos.genome.individuals import rows_in, sources
from genomeos.genome.regulatory import load_ccres

# a few values the literature has named, so the report can point at them
KNOWN_VALUES = {
    "OCA2": [
        ("rs12913832", "chr15", 28120472, "HERC2 intron 86 enhancer of OCA2; A brown, G blue (Sturm 2008)")
    ],
    "HERC2": [
        ("rs12913832", "chr15", 28120472, "HERC2 intron 86 enhancer of OCA2; A brown, G blue (Sturm 2008)")
    ],
}


def gene_span(ann: Annotation, gene: str, flank: int) -> tuple[Any, int, int]:
    g = ann.gene(gene)
    return g, max(0, g.locus.start - flank), g.locus.end + flank


def features(g, ccres: list, start: int, end: int) -> list[tuple[int, int, str]]:
    """(start, end, label) for exons and UTRs of the canonical transcript and cCREs in the span."""
    out: list[tuple[int, int, str]] = []
    ts = list(g.transcripts.values())
    canon = [t for t in ts if "Ensembl_canonical" in t.tags] or sorted(
        ts, key=lambda t: -(t.locus.end - t.locus.start)
    )
    if canon:
        cds = [(c.start, c.end) for c, _ in canon[0].cds]
        for e in canon[0].exons:
            coding = any(c0 < e.end and e.start < c1 for c0, c1 in cds)
            out.append((e.start, e.end, "exon (coding)" if coding else "exon (UTR or non-coding)"))
    for c in ccres:
        if c.end > start and c.start < end:
            out.append((c.start, c.end, f"cCRE {c.cls}"))
    return out


def label_for(pos0: int, feats: list[tuple[int, int, str]]) -> str:
    hits = [lab for s, e, lab in feats if s <= pos0 < e]
    exon = [h for h in hits if h.startswith("exon")]
    ccre = [h for h in hits if h.startswith("cCRE")]
    if exon:
        return exon[0] + (f", {ccre[0]}" if ccre else "")
    if ccre:
        return "intron, " + ccre[0]
    return "intron"


def read_variants(
    chrom: str, start: int, end: int, names: list[str] | None, root: Path | None = None
) -> dict[int, dict[str, Any]]:
    """SNVs by 1-based position across the persons with rows on this chromosome: ref, alt, and each
    person's genotype (het, hom); indels are counted apart."""
    by_pos: dict[int, dict[str, Any]] = {}
    people: list[str] = []
    for name, path, _ev in sources(chrom, root):
        if names and name not in names:
            continue
        people.append(name)
        for f in rows_in(path, start + 1, end):
            ref, alt = f[3], f[4].split(",")[0]
            gt = f[9].split(":")[0] if len(f) > 9 else "./."
            alleles = [a for a in gt.replace("|", "/").split("/") if a not in (".", "")]
            if not any(a != "0" for a in alleles):
                continue
            zyg = "hom" if alleles and all(a == alleles[0] for a in alleles) and alleles[0] != "0" else "het"
            rec = by_pos.setdefault(
                int(f[1]), {"ref": ref, "alt": alt, "snv": len(ref) == 1 and len(alt) == 1, "carriers": {}}
            )
            rec["carriers"][name] = zyg
    for rec in by_pos.values():
        rec["people"] = people
    return by_pos


def syntax_values(
    gene: str,
    chrom: str,
    names: list[str] | None = None,
    flank: int = 0,
    threshold: float = PHYLOP_THRESHOLD,
    url: str = PHYLOP_241_URL,
    progress=None,
) -> dict[str, Any]:
    t0 = time.time()
    gff = default_gencode({chrom})
    if gff is None:
        raise FileNotFoundError(f"no GENCODE models for {chrom}; run genomeos data fetch --chrom {chrom}")
    ann = Annotation.from_gff3(gff, {chrom})
    g, start, end = gene_span(ann, gene, flank)
    feats = features(g, load_ccres(chrom), start, end)
    variants = read_variants(chrom, start, end, names)
    people = next(
        (v["people"] for v in variants.values()),
        [n for n, _p, _e in sources(chrom) if not names or n in names],
    )
    snv_pos = sorted(p for p, v in variants.items() if v["snv"])
    bw = BigWig(url)
    try:
        span = bw.summarise(chrom, [(start, end)], threshold, progress=progress)[0]
        # one-base intervals at every SNV position (1-based -> 0-based)
        per_base = bw.summarise(chrom, [(p - 1, p) for p in snv_pos], threshold) if snv_pos else []
        mb = bw.src.bytes_fetched / 1e6
    finally:
        bw.close()
    rows = []
    in_syntax = 0
    for p, st in zip(snv_pos, per_base, strict=True):
        v = variants[p]
        phy = st.maximum if st.bases else None
        constrained = phy is not None and phy >= threshold
        in_syntax += int(constrained)
        rows.append(
            {
                "pos": p,
                "ref": v["ref"],
                "alt": v["alt"],
                "phylop": round(phy, 2) if phy is not None else None,
                "class": "value in syntax" if constrained else "value",
                "where": label_for(p - 1, feats),
                "carriers": v["carriers"],
            }
        )
    rows.sort(key=lambda r: -(r["phylop"] if r["phylop"] is not None else -99))
    known = []
    for rsid, kchrom, kpos, note in KNOWN_VALUES.get(gene.upper(), []):
        if kchrom != chrom:
            continue
        v = variants.get(kpos)
        known.append(
            {
                "id": rsid,
                "pos": kpos,
                "inside_span": start < kpos <= end,
                "note": note,
                "carriers": v["carriers"] if v else {},
                "reference_allele_in_everyone_else": v is None,
            }
        )
    bases = end - start
    return {
        "gene": g.symbol,
        "chrom": chrom,
        "strand": "+" if g.locus.strand is Strand.PLUS else "-",
        "span": [start, end],
        "flank": flank,
        "bases": bases,
        "people": people,
        "threshold": threshold,
        "syntax_bases": span.above,
        "syntax_fraction": round(span.above / span.bases, 4) if span.bases else None,
        "phylop_mean": round(span.mean, 3) if span.bases else None,
        "variant_positions": len(variants),
        "snv_positions": len(snv_pos),
        "indel_positions": len(variants) - len(snv_pos),
        "values_in_syntax": in_syntax,
        "values_per_kb": round(len(snv_pos) / bases * 1000, 2) if bases else None,
        "syntax_hit_rate": round(in_syntax / len(snv_pos), 4) if snv_pos else None,
        "expected_if_random": round(span.above / span.bases, 4) if span.bases else None,
        "rows": rows,
        "known_values": known,
        "features": {
            "exons": sum(1 for f in feats if f[2].startswith("exon")),
            "ccres": sum(1 for f in feats if f[2].startswith("cCRE")),
        },
        "cost": {"seconds": round(time.time() - t0, 1), "mb_fetched": round(mb, 2)},
        "evidence": {
            "syntax": "curated: Zoonomia phyloP over 241 placental mammals, constrained at >= 2.27 (5% FDR)",
            "values": "measured: each person's own variant calls (GIAB benchmarks for HG002, HG003, HG004)",
            "features": "curated: GENCODE canonical transcript; ENCODE cCRE registry",
        },
        "reading": (
            "syntax is what selection held still; a value is where people differ; a value in syntax is a "
            "change on protected sequence, the positions worth asking about first"
        ),
    }


OPEN_CONSENT = {"HG002", "HG003", "HG004"}  # GIAB samples released for open use; anyone else stays private


def save_path(gene: str, people: list[str], root: Path = Path("data/results")) -> Path:
    """Under data/results only when every person read is an open-consent GIAB sample; otherwise the
    result carries someone's genotypes and lives under that person's own git-ignored directory."""
    private = [n for n in people if n not in OPEN_CONSENT]
    if private:
        from genomeos.genome.individuals import ROOT

        return ROOT / private[0] / f"syntax_{gene.upper()}.json"
    return root / f"syntax_{gene.upper()}.json"
