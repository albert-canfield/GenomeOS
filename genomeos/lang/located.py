# SPDX-License-Identifier: Apache-2.0
# Part of the BioLang engine (language, IR, VM, standard library); see LICENSING.md.
"""Where everything is: the location analysis of a located BioLang program (v0.4 §4).

A program that declares one `compartment` is *located*, and then three facts hold rather than
being assumed:

1. Every gene and protein has a place. A missing location is an error, never a default.
2. A gene is transcribed in its own compartment and translated where its mRNA meets ribosomes,
   which for a nuclear gene means a transport (the nuclear pore) must carry the mRNA out.
3. A rule acts at the compartment of its target and can only read a source that is there, or
   sits in a membrane facing it. Anything else is a compile error rather than a silent effect.
4. A rule whose threshold is a concentration acts only where a compartment declares an absolute
   volume, because a concentration without a volume is not a number (§4.1, §10 decision 2).

This module is pure analysis over the IR: the parser turns what it returns into compile errors,
`bio check` prints it, and the located runtime uses the same routes to place its species.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass, field

from genomeos.ir import UNKNOWN, Compartment, Gene, Module, Protein, Transport

MRNA = "mRNA"


@dataclass(slots=True)
class Places:
    """The containment tree: who contains whom, and therefore what is adjacent to what."""

    by_id: dict[str, Compartment] = field(default_factory=dict)
    children: dict[str, list[str]] = field(default_factory=dict)
    roots: list[str] = field(default_factory=list)

    def neighbours(self, cid: str) -> list[str]:
        parent = self.by_id[cid].parent
        return ([parent] if parent else []) + list(self.children.get(cid, ()))

    def adjacent(self, a: str, b: str) -> bool:
        """Contained one in the other, or separated by exactly one membrane (which a transport crosses)."""
        if b in self.neighbours(a):
            return True
        return any(
            c.membrane and a in self.neighbours(c.id) and b in self.neighbours(c.id)
            for c in self.by_id.values()
        )

    def faces(self, cid: str) -> list[str]:
        """Where a species located in `cid` can be read: there, and both sides if it is a membrane."""
        return [cid, *self.neighbours(cid)] if self.by_id[cid].membrane else [cid]


def places(module: Module) -> Places:
    pl = Places()
    for c in module.compartments():
        pl.by_id[c.id] = c
    for c in pl.by_id.values():
        if c.parent:
            pl.children.setdefault(c.parent, []).append(c.id)
        else:
            pl.roots.append(c.id)
    return pl


@dataclass(slots=True)
class Layout:
    """The resolved places of a located program.

    gene_site: the compartment each gene is read in. mrna_sites: every compartment its mRNA reaches.
    synthesis: where each protein is born (per producing gene). reach: every compartment a protein can
    get to from where it is born, following the transports that recognise it. errors: what makes the
    program invalid; notes: what is legal but worth knowing (a declared place no route reaches)."""

    places: Places
    gene_site: dict[str, str] = field(default_factory=dict)
    mrna_sites: dict[str, list[str]] = field(default_factory=dict)
    synthesis: dict[str, list[str]] = field(default_factory=dict)
    reach: dict[str, set[str]] = field(default_factory=dict)
    errors: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)


def _walk(start: list[str], moves: list[Transport], allowed) -> list[str]:
    """Breadth-first over transports whose cargo test passes; returns compartments in visiting order."""
    seen, order, queue = set(start), list(start), deque(start)
    while queue:
        here = queue.popleft()
        for t in moves:
            if t.from_compartment == here and t.to_compartment not in seen and allowed(t):
                seen.add(t.to_compartment)
                order.append(t.to_compartment)
                queue.append(t.to_compartment)
    return order


def _gene_site(gene: Gene, pl: Places, errors: list[str]) -> str:
    chrom = gene.locus.chrom if gene.locus is not None else ""
    readers = [c.id for c in pl.by_id.values() if chrom and c.reads(chrom)]
    if len(gene.location) > 1:
        errors.append(f"gene {gene.id!r} is read in one compartment, not {gene.location}")
        return ""
    if gene.location:
        site = gene.location[0]
        if site not in pl.by_id:
            errors.append(f"gene {gene.id!r} is located in undeclared compartment {site!r}")
            return ""
        if chrom and not pl.by_id[site].reads(chrom):
            where = f"; {chrom} is read in {readers}" if readers else ""
            errors.append(
                f"gene {gene.id!r} on {chrom} is declared in {site!r}, which does not read {chrom}{where}"
            )
            return ""
        return site
    if len(readers) == 1:
        return readers[0]
    if len(readers) > 1:
        errors.append(f"gene {gene.id!r} on {chrom} could be read in {readers}; state its location")
    else:
        errors.append(f"gene {gene.id!r} has no location (and no compartment reads its locus)")
    return ""


def _check_tree(pl: Places, errors: list[str]) -> None:
    if len(pl.roots) != 1:
        errors.append(f"compartments form one tree with one root; roots: {sorted(pl.roots)}")
    for c in pl.by_id.values():
        if c.parent and c.parent not in pl.by_id:
            errors.append(f"compartment {c.id!r} has undeclared parent {c.parent!r}")
    for c in pl.by_id.values():  # a cycle never reaches a root
        seen, here = {c.id}, c.parent
        while here and here in pl.by_id:
            if here in seen:
                errors.append(f"compartment {c.id!r} is inside itself")
                break
            seen.add(here)
            here = pl.by_id[here].parent


def _check_transports(module: Module, pl: Places, errors: list[str]) -> None:
    proteins = {p.id for p in module.proteins()}
    for t in module.transports():
        ends = (t.from_compartment, t.to_compartment)
        missing = [e for e in ends if e not in pl.by_id]
        if missing:
            errors.append(f"transport {t.id!r} names undeclared compartments {missing}")
        elif not pl.adjacent(*ends):
            errors.append(f"transport {t.id!r} joins {ends[0]!r} and {ends[1]!r}, which are not adjacent")
        if t.via and t.via not in proteins:
            errors.append(f"transport {t.id!r} is gated by undeclared protein {t.via!r}")


def layout(module: Module) -> Layout:
    """Resolve every place in a located program; `errors` is empty for a valid one."""
    pl = places(module)
    out = Layout(places=pl)
    _check_tree(pl, out.errors)
    _check_transports(module, pl, out.errors)
    moves = module.transports()
    proteins = {p.id: p for p in module.proteins()}
    for g in module.genes():
        site = _gene_site(g, pl, out.errors)
        if site:
            out.gene_site[g.id] = site
            out.mrna_sites[g.id] = _walk([site], moves, lambda t, gid=g.id: t.carries(gid, mrna=True))
    for p in proteins.values():
        if not p.location:
            out.errors.append(f"protein {p.id!r} has no location")
        for c in p.location:
            if c not in pl.by_id:
                out.errors.append(f"protein {p.id!r} is located in undeclared compartment {c!r}")
    for r in module.rules:
        if r.action.value != "produces" or r.source not in out.gene_site or r.target not in proteins:
            continue
        sites = [c for c in out.mrna_sites[r.source] if pl.by_id[c].translation]
        if not sites:
            out.errors.append(
                f"{r.source!r} produces {r.target!r} but its mRNA reaches no compartment with "
                f"translation (from {out.gene_site[r.source]!r}; is an mRNA transport missing?)"
            )
            continue
        born = out.synthesis.setdefault(r.target, [])
        if sites[0] not in born:
            born.append(sites[0])  # the nearest ribosomes on the mRNA's route
    for p in proteins.values():
        start = out.synthesis.get(p.id) or list(p.location)  # no producer: it is there from the start
        got = _walk(start, moves, lambda t, pid=p.id, sig=tuple(p.signals): t.carries(pid, signals=sig))
        out.reach[p.id] = set(got)
        for c in p.location:
            if c in pl.by_id and c not in out.reach[p.id]:
                out.notes.append(
                    f"protein {p.id!r} is declared in {c!r} but no transport carries it there from "
                    f"{start}: it will be stranded (signals: {p.signals or 'none'})"
                )
    out.errors += _check_rules(module, out, proteins)
    return out


def rule_sites(module: Module, out: Layout, rule, proteins: dict[str, Protein]) -> list[str]:
    """The compartments where a rule acts: its target's, where the source can be read."""
    if rule.target in out.gene_site:
        targets = [out.gene_site[rule.target]]
    elif rule.target in proteins:
        targets = list(proteins[rule.target].location)
    else:
        return []
    src = proteins.get(rule.source)
    if src is None:
        readable = {out.gene_site[rule.source]} if rule.source in out.gene_site else set()
    else:
        readable = {f for c in src.location if c in out.places.by_id for f in out.places.faces(c)}
    return [c for c in targets if c in readable]


def _check_rules(module: Module, out: Layout, proteins: dict[str, Protein]) -> list[str]:
    errors = []
    for r in module.rules:
        if r.action.value == "produces":
            continue
        # a concentration threshold needs something to divide by, and a fraction of a cell is not it
        if r.threshold_unit:
            for c in rule_sites(module, out, r, proteins):
                if out.places.by_id[c].absolute_volume_fl is UNKNOWN:
                    errors.append(
                        f"rule {r.id!r} states its threshold as a concentration "
                        f"({r.threshold} {r.threshold_unit}) but compartment {c!r}, where it acts, "
                        "declares no absolute_volume; state one with its citation, or state the "
                        "threshold as an amount (BIOLANG-v0.4-ECONOMY.md §4.1)"
                    )
        known = (r.source in proteins or r.source in out.gene_site) and (
            r.target in proteins or r.target in out.gene_site
        )
        if not known or rule_sites(module, out, r, proteins):
            continue
        where = proteins[r.source].location if r.source in proteins else [out.gene_site.get(r.source, "?")]
        target_at = (
            proteins[r.target].location if r.target in proteins else [out.gene_site.get(r.target, "?")]
        )
        errors.append(
            f"rule {r.id!r} acts across compartments: {r.source!r} is in {where}, {r.target!r} in "
            f"{target_at}; they share no compartment and no membrane faces one"
        )
    return errors


def check_locations(module: Module) -> list[str]:
    """Compile errors of a located program; an empty list for an unlocated or a valid one."""
    return layout(module).errors if module.located else []
