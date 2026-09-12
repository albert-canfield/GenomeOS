# SPDX-License-Identifier: AGPL-3.0-or-later
"""One locus across species: does the module reconstruct in mouse, chicken and fish?

The last step of ROADMAP area J is the experiment the conversation proposed: take one
developmental locus with every layer available, read it in human and in other vertebrates,
and see whether the same grammar is there. Two loci are the first candidates. The ZRS
(VISTA hs2496, chr7:156,791,086-156,791,875) is the limb enhancer of SHH, a megabase away
inside an intron of LMBR1, conserved to fish and the textbook long-range pointer. The
HERC2 intron-86 enhancer of OCA2 carries rs12913832, Albert's eye-colour value slot.

For each locus and each species, Ensembl Compara's pairwise LASTZ alignment of the human
region is fetched (one call, under a second), giving the other species' locus, the
identity over aligned columns and the share of the human region that aligns at all. The
JASPAR scan (genome/motifs.py) then runs over both sides of each alignment, and a factor's
motif counts as conserved only where it hits the *same aligned site* in both: identity alone
makes nearly every motif "present" in a close species (a first version called 170 factors
conserved across the ZRS), so each human hit is mapped through the alignment columns. The
factors holding their site in every species are the locus's grammar. Both
constraint axes (Zoonomia phyloP, gnomAD Gnocchi) are read over the human region as in
`genomeos variation`. Everything fetched is cached under data/knowledge/across.

Evidence: the alignment is `curated` (Ensembl LASTZ_NET), a motif hit `predicted`, a
"conserved grammar" call `inferred` from the two. What the alignment cannot say is
whether the element is active in the other species; VISTA's mouse assay says it for the
ZRS (positive, limb), and for the rest that is the reader's job.
"""

from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

from genomeos.results import RESULTS_DIR, save_result

CACHE = Path("data/knowledge/across")
REST = (
    "https://rest.ensembl.org/alignment/region/homo_sapiens/{region}"
    "?method=LASTZ_NET&species_set=homo_sapiens&species_set={species}"
)
SPECIES = ("mus_musculus", "monodelphis_domestica", "gallus_gallus", "xenopus_tropicalis", "danio_rerio")
LOCI = {
    "ZRS": {
        "chrom": "chr7",
        "start": 156_791_086,
        "end": 156_791_875,
        "gene": "SHH",
        "note": (
            "VISTA hs2496, limb positive; the SHH limb enhancer inside LMBR1 intron 5, 1 Mb from its target"
        ),
    },
    "HERC2_OCA2": {
        "chrom": "chr15",
        "start": 28_120_106,
        "end": 28_121_063,
        "gene": "OCA2",
        "note": (
            "the HERC2 intron-86 enhancer of OCA2 around rs12913832 (chr15:28,120,472), eye colour's slot"
        ),
    },
}
EVIDENCE = {
    "alignment": "curated: Ensembl Compara LASTZ_NET pairwise alignment of the human region",
    "motif": "predicted: JASPAR 2026 log-odds hit at >= 85% of the matrix range",
    "grammar": "inferred: a factor whose motif holds the same aligned site in human and every species",
}


def _get(url: str, cache_key: str, timeout: int = 180) -> Any:
    CACHE.mkdir(parents=True, exist_ok=True)
    p = CACHE / f"{cache_key}.json"
    if p.exists():
        return json.loads(p.read_text())
    req = urllib.request.Request(
        url, headers={"Content-Type": "application/json", "User-Agent": "GenomeOS/0.9 (across)"}
    )
    with urllib.request.urlopen(req, timeout=timeout) as r:  # noqa: S310
        d = json.load(r)
    p.write_text(json.dumps(d))
    return d


def fetch_alignment(chrom: str, start: int, end: int, species: str) -> list[dict]:
    """The LASTZ_NET blocks of the human region against one species (1-based inclusive region)."""
    region = f"{chrom.removeprefix('chr')}:{start + 1}-{end}"
    url = REST.format(region=urllib.parse.quote(region), species=species)
    d = _get(url, f"{chrom}_{start}_{end}_{species}")
    return d if isinstance(d, list) else [d]


def summarise_alignment(blocks: list[dict], species: str, length: int) -> dict[str, Any]:
    """Identity over aligned columns, share of the human region covered, the other locus, the sequences."""
    human_cols = other_cols = matches = 0
    human_seq_parts: list[str] = []
    other_seq_parts: list[str] = []
    loci: list[str] = []
    for b in blocks:
        h = next((a for a in b.get("alignments", []) if a["species"] == "homo_sapiens"), None)
        o = next((a for a in b.get("alignments", []) if a["species"] == species), None)
        if not h or not o:
            continue
        for x, y in zip(h["seq"], o["seq"], strict=False):
            if x != "-":
                human_cols += 1
            if y != "-":
                other_cols += 1
            if x != "-" and y != "-" and x.upper() == y.upper():
                matches += 1
        # gapped, column for column: the same-site test maps positions through these columns
        human_seq_parts.append(h["seq"])
        other_seq_parts.append(o["seq"])
        loci.append(f"{o['seq_region']}:{o['start']}-{o['end']}({'+' if o['strand'] == 1 else '-'})")
    aligned = min(human_cols, other_cols)
    return {
        "species": species,
        "blocks": len(loci),
        "locus": ", ".join(loci) or None,
        "human_bases_aligned": human_cols,
        "coverage": round(human_cols / length, 3) if length else None,
        "identity": round(matches / human_cols, 3) if human_cols else None,
        "other_sequence": "".join(other_seq_parts),
        "human_aligned_sequence": "".join(human_seq_parts),
        "aligned_columns": aligned,
    }


SAME_SITE_WINDOW = 15  # a motif counts as the same site when its start maps within this many bases


def column_map(gapped: str) -> list[int]:
    """For each ungapped base of a gapped sequence, its column index."""
    return [i for i, c in enumerate(gapped) if c != "-"]


def base_at_column(gapped: str) -> dict[int, int]:
    """Column index to ungapped base index."""
    return {i: n for n, i in enumerate(column_map(gapped))}


def conserved_grammar(
    pairs: dict[str, tuple[str, str]], motifs, window: int = SAME_SITE_WINDOW
) -> dict[str, Any]:
    """Factors hitting the human element, and in which species they hit *the same aligned site*.

    Sequence identity alone makes almost every motif "present" in a close species: the honest
    question is whether the hit sits at the aligned position, so each hit's start is mapped through
    the alignment columns and must land within `window` bases of the other species' hit.
    """
    from genomeos.genome.motifs import scan

    per_species: dict[str, dict[str, int]] = {}
    human_hits: dict[str, tuple[float, int]] = {}
    for sp, (human_gapped, other_gapped) in pairs.items():
        h_seq, o_seq = human_gapped.replace("-", ""), other_gapped.replace("-", "")
        if not h_seq or not o_seq:
            continue
        hits = scan(motifs, [h_seq, o_seq])
        h_cols, o_base = column_map(human_gapped), base_at_column(other_gapped)
        for tf, per_seq in hits.items():
            if 0 not in per_seq:
                continue
            score, pos, _ = per_seq[0]
            if tf not in human_hits or score > human_hits[tf][0]:
                human_hits[tf] = (score, pos)
            if 1 not in per_seq:
                continue
            col = h_cols[pos] if pos < len(h_cols) else None
            if col is None:
                continue
            # the nearest column that exists in the other sequence, then its base index
            near = min(o_base, key=lambda c: abs(c - col)) if o_base else None
            expected = o_base.get(near) if near is not None else None
            if expected is None:
                continue
            if abs(per_seq[1][1] - expected) <= window:
                per_species.setdefault(tf, {})[sp] = per_seq[1][1]
    species = [sp for sp, (h, o) in pairs.items() if h.replace("-", "") and o.replace("-", "")]
    rows = []
    for tf, (score, pos) in human_hits.items():
        same = sorted(per_species.get(tf, {}))
        rows.append(
            {
                "factor": tf,
                "human_score": score,
                "human_position": pos,
                "same_site_in": same,
                "everywhere": bool(species) and len(same) == len(species),
            }
        )
    rows.sort(key=lambda r: (-len(r["same_site_in"]), -r["human_score"]))
    return {
        "factors_in_human": len(rows),
        "species_compared": species,
        "same_site_window": window,
        "everywhere": [r["factor"] for r in rows if r["everywhere"]],
        "rows": rows,  # every factor that hit human: a capped list made absent and unranked look alike
    }


def constraint(chrom: str, start: int, end: int) -> dict[str, Any]:
    from genomeos.attribution.constraint import PHYLOP_THRESHOLD, phylop_over_blocks
    from genomeos.attribution.variation import gnocchi_over, stats_dict

    ps, _ = phylop_over_blocks(chrom, [(start, end)], PHYLOP_THRESHOLD)
    gs, _ = gnocchi_over(chrom, [(start, end)])
    p = ps[0]
    return {
        "phylop": {
            "bases": p.bases,
            "mean": round(p.mean, 3) if p.mean is not None else None,
            "maximum": round(p.maximum, 3) if p.bases else None,
            "fraction_above": round(p.fraction_above, 3) if p.fraction_above is not None else None,
        },
        "gnocchi": stats_dict(gs[0]),
    }


def human_sequence(chrom: str, start: int, end: int) -> str | None:
    from genomeos.coords import Locus
    from genomeos.genome.fetch import REFERENCE
    from genomeos.genome.genome import Genome

    fa = REFERENCE / f"{chrom}.fa.gz"
    if not fa.exists():
        return None
    return str(Genome.from_fasta(fa).fetch(Locus(chrom, start, end)))


def build(name: str, species: tuple[str, ...] = SPECIES, motifs=None, progress=None) -> dict[str, Any]:
    if name not in LOCI:
        raise KeyError(f"unknown locus {name}; known: {', '.join(LOCI)}")
    loc = LOCI[name]
    chrom, start, end = loc["chrom"], loc["start"], loc["end"]
    t0 = time.time()
    length = end - start
    per_species: dict[str, dict] = {}
    for sp in species:
        try:
            blocks = fetch_alignment(chrom, start, end, sp)
            per_species[sp] = summarise_alignment(blocks, sp, length)
        except OSError as e:
            per_species[sp] = {"species": sp, "error": str(e)[:120], "coverage": 0.0, "identity": None}
        if progress:
            v = per_species[sp]
            progress(f"{name} against {sp}: coverage {v.get('coverage')} identity {v.get('identity')}")
    human = human_sequence(chrom, start, end)
    grammar: dict[str, Any] = {}
    if human:  # the local reference confirms the locus; the scan reads the aligned sequences
        if motifs is None:
            from genomeos.genome.motifs import load_motifs

            motifs = load_motifs()
        pairs = {
            sp: (v["human_aligned_sequence"], v["other_sequence"])
            for sp, v in per_species.items()
            if v.get("other_sequence") and v.get("human_aligned_sequence")
        }
        grammar = conserved_grammar(pairs, motifs)
    out = {
        "locus": name,
        **loc,
        "length": length,
        "species": {
            sp: {k: v for k, v in d.items() if k not in ("other_sequence", "human_aligned_sequence")}
            for sp, d in per_species.items()
        },
        "aligned_in": [sp for sp, v in per_species.items() if (v.get("coverage") or 0) > 0],
        "grammar": grammar,
        "constraint": constraint(chrom, start, end),
        "evidence": EVIDENCE,
        "cost": {"seconds": round(time.time() - t0, 1)},
    }
    return out


def run_and_save(name: str, results_dir: Path = RESULTS_DIR, progress=None) -> dict:
    out = build(name, progress=progress)
    save_result(f"across_{name}", out, results_dir)
    return out
