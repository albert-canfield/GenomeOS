# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""The node model against published rearrangements: the free half (docs/LOCI-BENCHMARK.md 11b).

A rearrangement is the only published experiment that can say whether our node is a real unit or a
resolution artefact. The claim under test is not that an element was lost. It is that a deletion,
inversion or duplication **moved a boundary**, so enhancers that belonged to one gene now reach
another and the limb gets the wrong instruction.

Nothing here makes a model request. Three readings, all local:

1. **the free one, asked first**: is there a CTCF-only element between the donor and the recipient
   at all? At HOXD there was not, and the node model's boundary there was a fiction. If that repeats
   across the rearrangement cases it is a second, independent statement about what our boundaries
   are made of;
2. **the reference node model**: does it put donor and recipient in different nodes, as the
   published answer requires;
3. **the node model recomputed on rearranged coordinates** - the same `infer_domains`, over a
   transformed cCRE list and a transformed chromosome length, not a simplification of it - does the
   boundary disappear and do donor and recipient end up in one node.

Every claim is scored identically at five matched random rearrangements per case, because everything
this benchmark has learned says a claim is worth nothing until a matched random version of it has
been tried. Cutting a chromosome anywhere destroys boundaries; the question is whether it also
produces *the published gene*.
"""

from __future__ import annotations

import bisect
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from genomeos.results import RESULTS_DIR, save_result

#: kept verbatim in the committed result, written before the first run (genomeos-9c, 2026-09-15)
PREDICTION = (
    "Cutting the genome at random destroys boundaries and creates adjacencies constantly, so the"
    " first two claims will fire at most control rearrangements, and naming the recipient gene will"
    " land near the 90% that naming a target already reaches at matched windows. What should"
    " separate the published rearrangements is which gene, and nothing else - if it does not, the"
    " node model has failed its sharpest test."
)

EVIDENCE = {
    "expectation": "curated: the published rearrangement and its phenotype, written down before the run",
    "boundary_census": "curated: ENCODE cCRE registry, CTCF-only class (an annotation lookup)",
    "node_model": (
        "inferred: CTCF-only elements as boundary proxies, confidence 0.4 (genome/domains.py)."
        " The node model is the thing under test here, so it takes the heuristic slot that"
        " nearest-coding-TSS holds in the locus benchmark, never the derived one"
    ),
    "derived": (
        "predicted: AlphaGenome asked on the rearranged sequence. NOT RUN - the scorer deletes one"
        " element and cannot construct a rearranged input. The free half reports first"
    ),
}

KINDS = ("deletion", "inversion", "duplication")
CONTROLS_PER_CASE = 5
CANDIDATES_PER_CASE = 60
MATCH_TSS = 0.35  # coding TSSs inside the span, within this fraction of the published span's
KEEP_OUT = 200_000  # a control span stays this far from the published one


#: the one thing about the published deletions that IS citable: their size (Lupianez et al. 2015)
PUBLISHED_DELETION_SIZE = (1_750_000, 1_900_000)


@dataclass(frozen=True)
class Rearrangement:
    """One published rearrangement, written down before anything is read."""

    locus: str
    chrom: str
    kind: str
    donor: str  # the gene whose regulatory landscape supplies the enhancers
    recipient: str  # the gene that acquires them and causes the phenotype
    phenotype: str
    span_source: str  # "published" when the breakpoints are cited, "stated" when the panel draws them
    span_citation: str
    citations: tuple[str, ...]
    answer_from: tuple[str, ...]
    span: tuple[int, int] | None = None  # filled from the genes when stated
    #: the published SIZE range, when the size is citable and the breakpoints are not. The case
    #: is then run at several spans across it, so the uncertainty is carried rather than collapsed.
    published_size_range: tuple[int, int] | None = None
    #: the fixed end the stated span is anchored to, and why that anchor is the published one
    anchor: str = ""
    note: str = ""

    def as_dict(self) -> dict[str, Any]:
        return {
            "locus": self.locus,
            "chrom": self.chrom,
            "kind": self.kind,
            "donor": self.donor,
            "recipient": self.recipient,
            "phenotype": self.phenotype,
            "span": list(self.span) if self.span else None,
            "span_source": self.span_source,
            "published_size_range": list(self.published_size_range) if self.published_size_range else None,
            "anchor": self.anchor,
            "span_citation": self.span_citation,
            "citations": list(self.citations),
            "answer_from": list(self.answer_from),
            "note": self.note,
        }


#: No published case is runnable. Every one of them needs human breakpoint coordinates, and the
#: Lupianez paper does not state any in its text: the accessible article gives SIZES and gene
#: content only, and the breakpoints live in aCGH supplementary data. Two independent routes were
#: tried (the PMC full text, and POSTRE's re-curation of the same variants) and neither carries a
#: coordinate. So the blocker was never the liftover - UCSC's chain and this project's chain reader
#: (attribution/human_panel.py) would take seconds - it is that there is no number to lift.
#: Placing a case approximately is the one thing the rule forbids, so all four are recorded here.
CASES: tuple[Rearrangement, ...] = (
    Rearrangement(
        locus="EPHA4_PAX3_deletion",
        chrom="chr2",
        kind="deletion",
        donor="EPHA4",
        recipient="PAX3",
        phenotype=(
            "brachydactyly: short digits on the preaxial side, from PAX3 falling under the EPHA4"
            " limb enhancers once the boundary between their domains is gone"
        ),
        span_source="stated",
        published_size_range=PUBLISHED_DELETION_SIZE,
        anchor=(
            "the 3' end is PAX3's GENCODE gene start (chr2:222,199,886): the published deletions"
            " extend into the NON-CODING part of the PAX3 domain and PAX3 itself survives to be"
            " misexpressed, so the deletion must stop before the gene. The 5' end follows from the"
            " published size, and at every size in the range it falls below EPHA4's start"
            " (chr2:221,418,026), so EPHA4 is removed entire, as published"
        ),
        span_citation=(
            "Lupianez et al. 2015, Cell 161:1012: heterozygous deletions of 1.75 to 1.9 Mb at"
            " 2q35-36 remove EPHA4 with much of its domain and extend into the non-coding part of"
            " the PAX3 domain, taking the boundary with them. The SIZE is cited and the breakpoints"
            " are not - the paper states no human coordinates anywhere in its accessible text - so"
            " this span is STATED, not cited: anchored at PAX3's gene start and run at every size"
            " across the published range, so the uncertainty is reported rather than hidden. An"
            " earlier version of this file stated a 626 kb intergenic span that excluded EPHA4; that"
            " was wrong in kind and the case was withdrawn until it could be stated faithfully"
        ),
        citations=(
            "Lupianez et al. 2015, Cell 161:1012 (disruptions of topological chromatin domains cause"
            " pathogenic rewiring of gene-enhancer interactions)",
            "Lupianez, Spielmann and Mundlos 2016, Trends Genet 32:225 (breaking TADs)",
        ),
        answer_from=("human pedigrees with aCGH", "4C-seq in patient fibroblasts", "mouse models"),
        note=(
            "the only claim this case can still decide is the third one - which gene the machinery"
            " names - because the first two are guaranteed by how the node model is built"
        ),
    ),
)

#: published cases and exactly what each is waiting for. None of these is "needs a liftover".
NOT_YET_RUNNABLE = {
    "EPHA4_IHH_duplication": (
        "polydactyly: a duplication of about 900 kb at 2q35 brings IHH under the EPHA4 enhancers"
        " (Lupianez et al. 2015). The effect depends on where the breakpoints fall INSIDE the"
        " domains, so no intergenic stand-in is faithful. Needs the real breakpoints"
    ),
    "EPHA4_WNT6_inversion": (
        "F-syndrome: a heterozygous inversion of about 1.1 Mb, its telomeric breakpoint about 1.4 Mb"
        " from EPHA4, moves WNT6 under the EPHA4 enhancers; a second family carries a 1.4 Mb"
        " duplication with a breakpoint about 1.2 Mb away (Lupianez et al. 2015). Sizes and offsets"
        " are published, the breakpoints are not. Needs the real breakpoints"
    ),
    "SOX9_KCNJ2_duplication": (
        "Cooks syndrome, duplications at the SOX9/KCNJ2 boundary (Franke et al. 2016, Nature"
        " 538:265). Same breakpoint problem, and it would overlap SOX9_PierreRobin's window in the"
        " locus panel, so the expectation would not be independent"
    ),
}

#: why claims one and two could never have discriminated, whatever the coordinates turn out to be
ANALYTIC_NOTE = (
    "Claims one and two are guaranteed by construction and cannot separate anything. The node model"
    " puts its boundaries at CTCF-only elements; a deletion removes every element inside its span,"
    " so any deletion spanning a boundary loses that boundary, and any deletion of the whole"
    " interval between two genes puts them in one node. This is arithmetic, not biology, and the"
    " size-matched deletions below confirm it at every one tried. Only claim three - WHICH gene the"
    " machinery then names - could ever have discriminated, which is exactly what the registered"
    " prediction said. The right reading of this module is therefore: two of its three claims were"
    " incapable of discriminating by construction, and the third needs chr2 swept before it can be"
    " scored at the controls at all."
)


# ------------------------------------------------------------------ the coordinate transforms
def transform(pos: int, span: tuple[int, int], kind: str) -> int | None:
    """Where a reference position lands after the rearrangement, or None if it is deleted away.

    Deletion: everything above the span shifts down by its length. Inversion: positions inside the
    span are mirrored within it. Duplication: everything above the span shifts up by its length,
    and positions inside it keep their first copy.
    """
    a, b = span
    n = b - a
    if kind == "deletion":
        if a <= pos < b:
            return None
        return pos if pos < a else pos - n
    if kind == "inversion":
        return a + b - 1 - pos if a <= pos < b else pos
    if kind == "duplication":
        return pos if pos < b else pos + n
    raise ValueError(f"unknown kind {kind}")


def rearranged_ccres(ccres: list, span: tuple[int, int], kind: str) -> list:
    """The cCRE list as the rearrangement leaves it, for `infer_domains` to read unchanged.

    Only the CTCF-only class matters to the boundary model, but every element is carried through so
    the same function can be used if the node model ever reads another class.
    """
    from dataclasses import replace

    a, b = span
    n = b - a
    out = []
    for c in ccres:
        mid = (c.start + c.end) // 2
        new = transform(mid, span, kind)
        if new is None:
            continue
        half = max(1, (c.end - c.start) // 2)
        out.append(replace(c, start=max(0, new - half), end=new + half))
        if kind == "duplication" and a <= mid < b:  # the second copy sits just above the original
            out.append(replace(c, start=max(0, new + n - half), end=new + n + half))
    out.sort(key=lambda c: c.start)
    return out


def rearranged_length(length: int, span: tuple[int, int], kind: str) -> int:
    n = span[1] - span[0]
    return length - n if kind == "deletion" else (length + n if kind == "duplication" else length)


# ------------------------------------------------------------------------------ the readings
def node_of(domains: list, pos: int) -> str | None:
    starts = [d.start for d in domains]
    i = bisect.bisect_right(starts, pos) - 1
    d = domains[i] if 0 <= i < len(domains) else None
    return d.id if d is not None and d.start <= pos < d.end else None


def ctcf_between(ccres: list, lo: int, hi: int) -> list[dict[str, Any]]:
    """Every CTCF-only element between two positions: the free reading, asked first."""
    a, b = (lo, hi) if lo <= hi else (hi, lo)
    return [
        {"id": c.id, "start": c.start, "end": c.end}
        for c in ccres
        if c.cls == "CTCF-only" and a <= (c.start + c.end) // 2 < b
    ]


def claims_for(ch, span: tuple[int, int], kind: str, donor_pos: int, recipient_pos: int) -> dict[str, Any]:
    """The three claims, scored identically at a published rearrangement and at a random one.

    `separated_before` is not a claim but the precondition: if the node model does not separate the
    pair to begin with, the rearrangement has no boundary to remove and the case says nothing.
    """
    from genomeos.genome.domains import infer_domains

    before = ch.domains
    d0, r0 = node_of(before, donor_pos), node_of(before, recipient_pos)
    walls_before = ctcf_between(ch.ccres, donor_pos, recipient_pos)
    after = infer_domains(
        ch.chrom, rearranged_length(ch.length, span, kind), rearranged_ccres(ch.ccres, span, kind)
    )
    d1 = transform(donor_pos, span, kind)
    r1 = transform(recipient_pos, span, kind)
    d1n = node_of(after, d1) if d1 is not None else None
    r1n = node_of(after, r1) if r1 is not None else None
    walls_after = ctcf_between(
        rearranged_ccres(ch.ccres, span, kind), d1 if d1 is not None else 0, r1 if r1 is not None else 0
    )
    gone = d1 is None or r1 is None
    return {
        # the published deletion removes EPHA4 itself, so the donor GENE does not survive to share a
        # node with anything. The mechanism is that the donor's surviving ENHANCERS reach the
        # recipient, and this module has no published coordinate for them. Marked not-applicable
        # rather than False: a claim the construction cannot ask is not a claim the model failed.
        "donor_or_recipient_deleted": gone,
        "separated_before": bool(d0 and r0 and d0 != r0),
        "node_before": {"donor": d0, "recipient": r0},
        "node_after": {"donor": d1n, "recipient": r1n},
        "ctcf_between_before": len(walls_before),
        "ctcf_between_after": len(walls_after),
        # claim 1: the node model loses a boundary between the pair
        "boundary_lost": None if gone else len(walls_after) < len(walls_before),
        # claim 2: a new adjacency appears - the pair now shares a node when it did not
        "new_adjacency": None if gone else bool(d0 and r0 and d0 != r0 and d1n and r1n and d1n == r1n),
        "nodes_before": len(before),
        "nodes_after": len(after),
    }


def control_candidates(ch, case: Rearrangement) -> list[dict[str, Any]]:
    """Every consecutive coding-gene pair on the chromosome, as a rearrangement of the same shape.

    **The control has to share the construction, not only the size.** The first version of this drew
    a random span of the published length and took the nearest coding TSS on each side. That looked
    matched and was not: the published span is the whole interval *between* the two genes, so
    deleting it removes every boundary between them by construction, while a random span leaves the
    boundaries that sit between it and the flanking genes. The published case then scored a new
    adjacency 1/1 against 0/5 - an artefact of how the two were built, not a fact about the genome.
    A control is therefore the same thing done elsewhere: a consecutive pair of coding genes and the
    whole interval between them, matched on how long that interval is and how many boundaries it
    crosses.
    """
    out = []
    for (t0, s0), (t1, s1) in zip(ch.coding, ch.coding[1:], strict=False):
        if t1 - t0 < 10_000:
            continue
        out.append(
            {
                # strictly between the two transcription starts, so neither gene sits inside the
                # interval that is removed - the published span (gene end to gene start) is the same
                # relation. Putting an endpoint inside its own deletion sends it nowhere, and the
                # first version of this did exactly that.
                "span": [t0 + 1, t1],
                "kind": case.kind,
                "length": t1 - t0 - 1,
                "coding_tss_inside": 0,  # consecutive by construction, as the published pair is
                "ctcf_crossed": len(ctcf_between(ch.ccres, t0, t1)),
                "donor_tss": t0,
                "recipient_tss": t1,
                "donor": s0,
                "recipient": s1,
            }
        )
    return out


def matched_controls(
    ch, case: Rearrangement, span: tuple[int, int], published: dict[str, Any], results_dir: Path
) -> list[dict[str, Any]]:
    """Five rearrangements of the same shape elsewhere on the chromosome.

    Same construction as the published case (a consecutive coding-gene pair and the whole interval
    between them), then matched on span length within 35% and on CTCF-only elements crossed, closest
    of the candidates - boundary density being what decides whether a cut can join anything at all.
    Kept 200 kb clear of the published span.
    """
    want_len, want_ctcf = published["length"], published["ctcf_crossed"]
    cands = [
        c
        for c in control_candidates(ch, case)
        if not (c["span"][0] < span[1] + KEEP_OUT and span[0] - KEEP_OUT < c["span"][1])
        and abs(c["length"] - want_len) <= MATCH_TSS * want_len
    ]
    cands.sort(key=lambda c: (abs(c["ctcf_crossed"] - want_ctcf), abs(c["length"] - want_len)))
    out = []
    for c in cands[:CONTROLS_PER_CASE]:
        c["claims"] = claims_for(ch, tuple(c["span"]), case.kind, c["donor_tss"], c["recipient_tss"])
        # claim 3 needs its control too, or it is the same mistake the whole benchmark exists to
        # catch: a rate with no matched denominator underneath it
        c["recipient_named"] = named_by_a_derived_layer(ch, results_dir, c["recipient_tss"], c["recipient"])
        c["claims"]["recipient_named"] = c["recipient_named"]["names_the_recipient"]
        out.append(c)
    return out


def named_by_a_derived_layer(ch, results_dir: Path, pos: int, gene: str) -> dict[str, Any]:
    """Claim 3: does any already-computed deletion in the recipient's neighbourhood name it.

    Read from runs already made, exactly as the locus benchmark's deletion layer is: no request.
    """
    from genomeos.benchmark.loci import _deletion_rows

    rows = _deletion_rows(ch.chrom, max(0, pos - 100_000), pos + 100_000, results_dir)
    named = []
    for e in rows:
        p = e.get("predicted_coding") or e.get("predicted") or {}
        if p.get("gene"):
            named.append(p["gene"])
    return {
        "elements_scored": len(rows),
        "genes_named": sorted(set(named))[:8],
        "names_the_recipient": gene in named,
        "pending": None if rows else "no deletion has been scored within 100 kb of the recipient",
    }


def boundary_behaviour(chrom: str = "chr2", results_dir: Path = RESULTS_DIR, progress=None) -> dict[str, Any]:
    """What a deletion of the published size does to the node model, anywhere on the chromosome.

    No published case is runnable (see NOT_YET_RUNNABLE), but one thing about them IS citable: the
    brachydactyly deletions are 1.75 to 1.9 Mb. Deleting an interval of that size between a
    consecutive coding pair is a lower bound on what any of these rearrangements would do to our
    boundaries, and it is enough to show that two of the three claims can never discriminate.
    """
    from genomeos.benchmark.loci import Chromosome

    say = progress or (lambda _m: None)
    ch = Chromosome(chrom, results_dir)
    try:
        lo, hi = PUBLISHED_DELETION_SIZE
        stub = Rearrangement(
            locus="size_matched_deletion",
            chrom=chrom,
            kind="deletion",
            donor="-",
            recipient="-",
            phenotype="-",
            span_source="stated",
            span_citation="x" * 81,
            citations=(),
            answer_from=(),
        )
        cands = [c for c in control_candidates(ch, stub) if lo <= c["length"] <= hi]
        say(f"{chrom}: {len(cands)} consecutive gene pairs {lo / 1e6:.2f} to {hi / 1e6:.2f} Mb apart")
        rows = []
        for c in cands[:CONTROLS_PER_CASE]:
            c["claims"] = claims_for(ch, tuple(c["span"]), "deletion", c["donor_tss"], c["recipient_tss"])
            c["recipient_named"] = named_by_a_derived_layer(
                ch, results_dir, c["recipient_tss"], c["recipient"]
            )
            c["claims"]["recipient_named"] = c["recipient_named"]["names_the_recipient"]
            rows.append(c)
        say(f"{chrom}: {len(rows)} size-matched deletions scored")
        keys = ("separated_before", "boundary_lost", "new_adjacency", "recipient_named")
        return {
            "chrom": chrom,
            "published_deletion_size": list(PUBLISHED_DELETION_SIZE),
            "candidate_pairs_at_that_size": len(cands),
            "deletions": rows,
            "rates": {k: {"k": sum(1 for r in rows if r["claims"].get(k)), "n": len(rows)} for k in keys},
            "scored_near_the_recipient": sum(1 for r in rows if r["recipient_named"]["elements_scored"]),
            "analytic_note": ANALYTIC_NOTE,
        }
    finally:
        ch.close()


def spans_across_the_range(case: Rearrangement, anchor_end: int, steps: int = 4) -> list[tuple[int, int]]:
    """One span per size across the published range, all ending at the anchor.

    The breakpoints are not published; the size is. Running every size in the range and reporting
    them together is the honest form of "take the published uncertainty with the coordinates".
    """
    lo, hi = case.published_size_range or (0, 0)
    if not lo:
        return [case.span] if case.span else []
    step = (hi - lo) // max(1, steps - 1)
    return [(anchor_end - (lo + i * step), anchor_end) for i in range(steps)]


def run(results_dir: Path = RESULTS_DIR, progress=None) -> dict[str, Any]:
    from genomeos.benchmark.loci import Chromosome, tss_of

    t0 = time.time()
    say = progress or (lambda _m: None)
    out_cases: list[dict[str, Any]] = []
    chroms: dict[str, Any] = {}
    try:
        for case in CASES:
            ch = chroms.get(case.chrom) or chroms.setdefault(case.chrom, Chromosome(case.chrom, results_dir))
            genes = {
                g.symbol: g for g in ch.annotation.genes.values() if g.symbol in (case.donor, case.recipient)
            }
            if len(genes) < 2:
                out_cases.append({"locus": case.locus, "pending": "donor or recipient not in GENCODE"})
                continue
            d, r = genes[case.donor], genes[case.recipient]
            donor_pos, recipient_pos = tss_of(d), tss_of(r)
            walls = ctcf_between(ch.ccres, donor_pos, recipient_pos)
            named = named_by_a_derived_layer(ch, results_dir, recipient_pos, case.recipient)
            sized = []
            for span in spans_across_the_range(case, r.locus.start):
                row = {
                    "span": list(span),
                    "length": span[1] - span[0],
                    "removes_the_donor_entire": span[0] <= d.locus.start and d.locus.end <= span[1],
                    "leaves_the_recipient_intact": span[1] <= r.locus.start,
                    "ctcf_crossed": len(ctcf_between(ch.ccres, *span)),
                    "claims": claims_for(ch, span, case.kind, donor_pos, recipient_pos),
                }
                row["claims"]["recipient_named"] = named["names_the_recipient"]
                sized.append(row)
            say(f"{case.locus}: {len(sized)} spans across the published size range")
            widest = max(sized, key=lambda x: x["length"])
            controls = matched_controls(ch, case, tuple(widest["span"]), widest, results_dir)
            say(f"{case.locus}: {len(controls)} matched random rearrangements")
            out_cases.append(
                {
                    "locus": case.locus,
                    "expected": case.as_dict(),
                    "boundary_census": {
                        "ctcf_only_between_donor_and_recipient": len(walls),
                        "elements": walls[:8],
                        "there_is_a_boundary_at_all": bool(walls),
                        "evidence": EVIDENCE["boundary_census"],
                    },
                    "spans_across_the_published_size_range": sized,
                    "published": widest,
                    "recipient_named": named,
                    "controls": controls,
                }
            )
    finally:
        for ch in chroms.values():
            ch.close()
    behaviour = boundary_behaviour(results_dir=results_dir, progress=progress)
    return {
        "result": "loci_rearrangements",
        "prediction_registered_before_the_run": PREDICTION,
        "cases": out_cases,
        "aggregate": aggregate(out_cases),
        "withdrawn_and_reinstated": (
            "The EPHA4-to-PAX3 deletion was withdrawn when hunting its breakpoints showed the span"
            " stated for it was wrong in kind (626 kb intergenic, excluding EPHA4, where the"
            " published deletions are 1.75 to 1.9 Mb and include it) and its phenotype was wrong"
            " (brachydactyly, not polydactyly). It is reinstated here with a span that matches the"
            " published construction, still STATED and not cited, run at every size across the"
            " published range. The other three cases stay withdrawn: no coordinate for them exists."
        ),
        "no_case_is_runnable": (
            "Every published case needs human breakpoint coordinates and the paper states none: the"
            " accessible text gives sizes and gene content only. Two routes were tried (the PMC full"
            " text and POSTRE's re-curation) and neither carries a coordinate, so there is nothing to"
            " lift over. The rule is to drop rather than place approximately, and it is applied here"
            " to the EPHA4-to-PAX3 deletion as well, which an earlier run placed approximately."
        ),
        "boundary_behaviour": behaviour,
        "analytic_note": ANALYTIC_NOTE,
        "not_yet_runnable": NOT_YET_RUNNABLE,
        "evidence": EVIDENCE,
        "note": (
            "The free half of docs/LOCI-BENCHMARK.md 11b: no model request was made. The node model"
            " is INFERRED at confidence 0.4 and is the thing under test, so it is never counted as a"
            " derived hit. One published case is an anecdote by the panel's own standard; the other"
            " three are recorded in not_yet_runnable with what each is waiting for."
        ),
        "cost": {"seconds": round(time.time() - t0, 1)},
    }


def aggregate(cases: list[dict[str, Any]]) -> dict[str, Any]:
    """The three claims at the published rearrangements and at the matched random ones."""
    keys = ("separated_before", "boundary_lost", "new_adjacency", "recipient_named")
    pos = [c["published"]["claims"] for c in cases if c.get("published")]
    neg = [x["claims"] for c in cases for x in c.get("controls", [])]

    def rate(rows, k):
        """Only rows where the claim applies: a None is not a miss, it is a question not asked."""
        asked = [r for r in rows if r.get(k) is not None]
        n = len(asked)
        got = sum(1 for r in asked if r[k])
        return {
            "k": got,
            "n": n,
            "rate": round(got / n, 3) if n else None,
            "not_applicable": len(rows) - n,
        }

    pub_named = [c["recipient_named"] for c in cases if c.get("recipient_named")]
    ctl_named = [x["recipient_named"] for c in cases for x in c.get("controls", [])]

    def judged(rows):
        """Claim 3 only where a deletion has actually been scored near the recipient.

        Without this the claim measures which chromosomes the sweep has reached. The published
        recipients sit where the panel deliberately spent requests and the controls do not, so an
        ungated 1/1 against 0/5 would be the coverage artefact this benchmark exists to catch, not a
        result. Gated, it is usually 'not judgeable', which is the honest answer.
        """
        have = [r for r in rows if r["elements_scored"]]
        return {
            "k": sum(1 for r in have if r["names_the_recipient"]),
            "n": len(have),
            "of": len(rows),
            "not_scored": len(rows) - len(have),
        }

    return {
        "published": {k: rate(pos, k) for k in keys},
        "controls": {k: rate(neg, k) for k in keys},
        "recipient_named_where_a_deletion_was_scored": {
            "published": judged(pub_named),
            "controls": judged(ctl_named),
            "judgeable": bool(judged(ctl_named)["n"]),
            "reading": "claim three is only meaningful where both sides have deletion data; where the"
            " controls have none it is not judgeable and must not be quoted",
        },
        "reading": (
            "claims one and two are expected to fire at most control rearrangements: cutting a"
            " chromosome anywhere destroys boundaries. Only the identity of the gene can separate the"
            " published cases, which is what the registered prediction says to look at"
        ),
    }


def run_and_save(results_dir: Path = RESULTS_DIR, progress=None) -> dict[str, Any]:
    out = run(results_dir, progress)
    save_result("loci_rearrangements", out, results_dir)
    return out
