# SPDX-License-Identifier: AGPL-3.0-or-later
"""Exploratory, not pre-registered (docs/DATA.md, "Phasing the trio candidates on the Q100 assembly"):
how often the trio's supported de novo SNV candidates, and a control of one-parent inherited SNVs,
are carried by any of the 88 HPRC haplotypes in the local human panel (data/knowledge/human_panel,
UCSC's hprc cactus90way, which does not include HG002). A missed inherited polymorphism should be seen
there about as often as the control; a new mutation should not.

Writes data/results/trio_q100_population.json: counts per group, never a position, allele or genotype.
"""

import gzip
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.genome import individuals as ind  # noqa: E402
from genomeos.results import save_result  # noqa: E402
from scripts.trio_q100_phase import manifest  # noqa: E402

C = Path("data/cache/q100")
VCF = C / "GRCh38_HG2-T2TQ100-V1.1_dipcall-z2k.dip.vcf.gz"
DIP = C / "GRCh38_HG2-T2TQ100-V1.1_dipcall-z2k.dip.bed"
MOS = C / "mergedQ100denovomosaicexclusions_240506_GRCh38_withoverlaprepeats.bed.gz"
chroms = [f"chr{i}" for i in range(1, 23)]


def bed(p):
    iv = {}
    with ind.open_variants(p) as fh:
        for l in fh:
            f = l.split("\t")
            iv.setdefault(f[0], []).append((int(f[1]), int(f[2])))
    return {c: (sorted(v), [a for a, _ in sorted(v)]) for c, v in iv.items()}


dip, mos = bed(DIP), bed(MOS)
calls = ind.read_assembly_calls(VCF, set(chroms), lambda c: ind._reference_base_reader(c)[0])
tally = {}  # group -> [n, seen_in_panel, not_covered]


def add(group, seen):
    t = tally.setdefault(group, [0, 0, 0])
    t[0] += 1
    if seen is None:
        t[2] += 1
    elif seen:
        t[1] += 1


for chrom in chroms:
    pc, pf, pm = (ind.vcf_path(n, chrom) for n in ("HG002", "HG003", "HG004"))
    b, g = ind._reference_base_reader(chrom)
    gc, gf, gm = ind._genotypes(pc, b), ind._genotypes(pf, b), ind._genotypes(pm, b)
    g.close()
    rf, rm, rc = (ind.load_regions(n, chrom) for n in ("HG003", "HG004", "HG002"))
    sf, sm, sc = [x for x, _ in rf], [x for x, _ in rm], [x for x, _ in rc]
    di, ds = dip[chrom]
    mi, ms = mos.get(chrom, ([], []))
    want = {}
    k = 0
    for key, zyg in gc.items():
        if len(key[1]) != 1 or len(key[2]) != 1:
            continue
        p0 = key[0] - 1
        if not (
            ind._inside(rf, sf, p0)
            and ind._inside(rm, sm, p0)
            and ind._inside(rc, sc, p0)
            and ind._inside(di, ds, p0)
        ):
            continue
        fh_, mh = key in gf, key in gm
        if not fh_ and not mh:
            cls = ind.assembly_class(key, calls.get(chrom))
            cls = {"hap1": "father", "hap2": "mother"}.get(cls, cls)
            inm = "in_mosaic" if ind._inside(mi, ms, p0) else "out_mosaic"
            want[p0] = (key[2], f"cand_{cls}_{inm}")
        elif fh_ != mh and zyg == "het":
            k += 1
            if k % 100 == 0:
                want[p0] = (key[2], "control_one_parent")
    meta = json.loads(Path(f"data/knowledge/human_panel/{chrom}/meta.json").read_text())
    keep = [i for i, a in enumerate(meta["assemblies"]) if a != "hs1"]
    found = set()
    with gzip.open(f"data/knowledge/human_panel/{chrom}/sites.tsv.gz", "rt") as fh:
        for l in fh:
            p, r, a = l.rstrip("\n").split("\t")
            p = int(p)
            if p in want:
                alt, grp = want[p]
                found.add(p)
                add(grp, any(a[i] == alt for i in keep if i < len(a)))
    for p, (_alt, grp) in want.items():
        if p not in found:
            add(grp, False)  # no panel site: no assembly differs from the reference here
    print(chrom, file=sys.stderr)
out = {
    "groups": {g: {"snvs": t[0], "seen_in_panel": t[1]} for g, t in sorted(tally.items())},
    "panel": "88 HPRC year-1 haplotypes (hs1 left out), UCSC hg38 hprc cactus90way, as distilled in "
    "data/knowledge/human_panel; seen = at least one haplotype carries the same alternative base",
    "control": "every 100th child heterozygous SNV carried by exactly one parent, same regions",
    "registered": False,
    "date": time.strftime("%Y-%m-%d"),
}
panel = [Path(f"data/knowledge/human_panel/{c}/{f}") for c in chroms for f in ("sites.tsv.gz", "meta.json")]
save_result(
    "trio_q100_population",
    out,
    manifest=manifest(
        {
            "control_every_nth_one_parent_het_snv": 100,
            "panel_assemblies_left_out": ["hs1"],
            "chromosomes": chroms,
        },
        extra_inputs=[
            mf.files_entry("data/knowledge/human_panel chr1-chr22 sites.tsv.gz and meta.json", panel)
        ],
        extra_sources=[
            {
                "accession": "UCSC hg38 hprc cactus90way (88 HPRC year-1 haplotypes), as distilled in "
                "data/knowledge/human_panel",
                "version": "hprc cactus90way",
            }
        ],
    ),
)
print(json.dumps(out["groups"], indent=1))
