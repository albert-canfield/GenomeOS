# SPDX-License-Identifier: AGPL-3.0-or-later
"""A species-tree null for area J: what shared ancestry alone gives a library pair, or a library.

Area J compares the libraries across the 355 species of `origin_presence_genome_wide` (one bit per
gene per species: an Ensembl Compara orthologue of the human gene is there or not). Two of its
readings count agreement between presence profiles species by species, and species are not
independent observations: a clade of sixty fish that lost a gene in one ancestor is one event and
sixty columns. GRAMMAR-BY-COMPARISON §11, §14 and §19 each record that the project "holds no tree"
and that a null along one is what would let a pair be called a dependency. This module holds the
tree and the two nulls.

THE TREE. Ensembl Compara's own species trees, the ones its protein trees (and so the orthologue
calls the presence matrix is made of) were built on: the pan-taxonomic tree of Ensembl Genomes
release 63 (175 species: bacteria, archaea, protists, plants, fungi, invertebrates and a few
vertebrates) and the vertebrate tree of Ensembl release 116 (200 species). The vertebrate tree's
Vertebrata clade replaces the pan tree's, and twelve species whose taxonomy was renamed since the
production name are mapped by an explicit table (`RENAMES`). Both files are topology only (the pan
one has no branch lengths), so branch lengths are a modelling choice stated in advance: Grafen's
(1989, rho 1: a node's height is its number of descendant tips minus one) as the registered one,
every branch of length one as the stated sensitivity.

NULL A, LIBRARY PAIRS (§14). §14 correlates two libraries' residual profiles (share of members
present per species, minus the genome's share) on the members exclusive to each, and compares the
correlation with pairs of random gene sets matching each side's origin-stratum make-up. That null
holds age. It does not hold the tree: a random set's profile is close to its age class's average, so
the null's spread is small, and a library's own departure from its age expectation -- which is one
trait evolving along the tree, autocorrelated within clades -- is compared as if each species were a
fresh observation. Null A splits each side's profile into the age expectation E (the exact mean of
§14's own matched null) and the library's deviation D = profile - E; it keeps E, replaces D with an
independent Brownian motion along the tree whose rate is D's own (the REML estimate under the tree's
covariance), and reads §14's statistic on the result. Two deviations that each sit in one deep clade
correlate easily under this null, which is the point: shared ancestry alone produces them.

NULL B, LIBRARY COHERENCE (§19). §19 scores a library by the mean Jaccard of its member pairs'
presence sets and compares it with random genes matched on prevalence. Null B puts each gene on the
tree: it is born at the most recent common ancestor of human and the species carrying it (Dollo: one
gain), and lost along the branches below independently, with a per-branch loss probability read from
every gene's Dollo reconstruction (the share of genes present at the branch's parent that are lost on
it, smoothed by half a gene) raised to a per-gene power fitted so that the gene's expected number of
species equals its observed one. Losses on the path to human are impossible by construction, since
every gene is present in human. The library's statistic is its observed mean pair Jaccard minus the
mean over simulated draws of its members (what each member's own birth and loss rate predict along
the tree); the null is the same statistic for random gene sets matched member by member on birth node.
The model removes what birth clade and loss rate explain; the real-gene sets remove whatever the
model misses for every gene alike, which §19 showed is the failure that matters ("a null that holds
the obvious covariate fixed can still be beaten by an artefact").

Nothing here reads a result under either null: the module computes, `scripts/species_null_register.py`
registers what will be read and how, and `scripts/species_null_run.py` reads it once.
"""

from __future__ import annotations

import hashlib
import math
import re
import urllib.request
from dataclasses import dataclass, field
from itertools import combinations
from pathlib import Path
from typing import Any

import numpy as np

TREE_DIR = Path("data/knowledge/species_tree")
PAN_TREE_URL = (
    "https://ftp.ensemblgenomes.ebi.ac.uk/pub/pan_ensembl/release-63/compara/species_trees/"
    "default_protein-trees_default.nh"
)
VERT_TREE_URL = (
    "https://ftp.ensembl.org/pub/release-116/compara/species_trees/default_protein-trees_default.nh"
)
PAN_TREE = TREE_DIR / "pan_default_protein-trees_default.nh"
VERT_TREE = TREE_DIR / "vert_default_protein-trees_default.nh"
#: sha256 of the two files as downloaded on 2026-10-10; their BSD `sum` matches the CHECKSUMS file
#: Ensembl publishes beside each (59496 8 and 61092 6).
TREE_SHA256 = {
    PAN_TREE.name: "e3f0237f4360a09d01c2a754412563e135d00d45c2c3dedb3cc5f7f4496a63b7",
    VERT_TREE.name: "2033c3eb2a6ab20f530737ffc66babb8e5e39edbb396ec6d315960db65142e66",
}
HUMAN = "homo_sapiens"
#: The two tips whose common ancestor is Vertebrata in the vertebrate tree.
VERTEBRATA_ANCHORS = ("homo_sapiens", "petromyzon_marinus")
#: Production names (assembly suffix stripped) whose species the trees spell under a newer name.
RENAMES = {
    "ralstonia_solanacearum_gmi1000": "ralstonia_pseudosolanacearum_gmi1000",
    "rhizobium_leguminosarum_bv_viciae_3841": "rhizobium_johnstonii_3841",
    "terrapene_carolina_triunguis": "terrapene_triunguis",
    "neovison_vison": "neogale_vison",
    "physeter_catodon": "physeter_macrocephalus",
    "desulfovibrio_vulgaris_str_hildenborough": "nitratidesulfovibrio_vulgaris_str_hildenborough",
    "giardia_intestinalis": "giardia_duodenalis_strain_wb_c6",
    "lactobacillus_plantarum_wcfs1": "lactiplantibacillus_plantarum_wcfs1",
    "mycoplasma_pneumoniae_m129": "mycoplasmoides_pneumoniae_m129",
    "rhodobacter_sphaeroides_2_4_1": "cereibacter_sphaeroides_2_4_1",
    "thermosynechococcus_elongatus_bp_1": "thermosynechococcus_vestitus_bp_1",
    "ustilago_maydis": "mycosarcoma_maydis_strain_521",
}

MIN_MEMBERS = 8  # §14's and §19's own floor, kept so the readings are the same readings
# null A
A_SIMS = 50_000  # Brownian draws per bank; the smallest attainable p is 1 / (A_SIMS + 1)
A_SEED = 20261010
A_CAL_SEED = 20261011
A_HALVES_SEED = 20261012
Q_ALPHA = 0.05
CAL_MAX_SHARE = 0.10  # random age-matched pairs at p <= 0.05 above this share: the null is not calibrated
HALVES_MIN_SHARE = 0.50  # libraries seen against their own halves below this share: low power
# null B
B_SIMS = 20  # model draws per gene
B_SETS = 200  # random birth-matched gene sets per library; §19 drew 20, and its classes are kept below
B_ALPHA = 0.05
B_SEED = 20261013
B_SIM_SEED = 20261014
LOSS_PSEUDO = (0.5, 1.0)  # p = (losses + 0.5) / (at risk + 1): no branch at exactly 0 or 1
B_TOLERANCE = 0.005  # observed mean pair Jaccard against §19's stored (sampled for large libraries)


# ------------------------------------------------------------------------------------------- tree


@dataclass(eq=False)
class Node:
    name: str = ""
    length: float | None = None
    children: list[Node] = field(default_factory=list)
    parent: Node | None = field(default=None, repr=False)


def parse_newick(text: str) -> Node:
    """A Newick string as a tree of Nodes. Names may hold any character but `(),:;`; quotes are dropped."""
    s = text.strip()
    if s.endswith(";"):
        s = s[:-1]
    root = Node()
    node = root
    i, n = 0, len(s)
    while i < n:
        c = s[i]
        if c == "(":
            child = Node(parent=node)
            node.children.append(child)
            node = child
            i += 1
        elif c == ",":
            if node.parent is None:
                raise ValueError(f"',' outside any parenthesis at {i}")
            child = Node(parent=node.parent)
            node.parent.children.append(child)
            node = child
            i += 1
        elif c == ")":
            if node.parent is None:
                raise ValueError(f"unbalanced ')' at {i}")
            node = node.parent
            i += 1
        else:
            j = i
            while j < n and s[j] not in "(),:;":
                j += 1
            label = s[i:j].strip().strip("'")
            if label:
                node.name = label
            i = j
            if i < n and s[i] == ":":
                j = i + 1
                while j < n and s[j] not in "(),;":
                    j += 1
                node.length = float(s[i + 1 : j])
                i = j
    if node is not root:
        raise ValueError("unbalanced '(': the string ends inside a clade")
    return root


def norm(name: str) -> str:
    """Lower case, every run of other characters one underscore: `Mus_musculus_reference_CL57BL6_strain`."""
    return re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")


def preorder(root: Node) -> list[Node]:
    out, stack = [], [root]
    while stack:
        v = stack.pop()
        out.append(v)
        stack.extend(reversed(v.children))
    return out


def leaves(root: Node) -> list[Node]:
    return [v for v in preorder(root) if not v.children]


def mrca(nodes: list[Node]) -> Node:
    """The most recent common ancestor of the given nodes."""
    if not nodes:
        raise ValueError("no nodes")

    def path(v: Node | None) -> list[Node]:
        out = []
        while v is not None:
            out.append(v)
            v = v.parent
        return out[::-1]

    paths = [path(v) for v in nodes]
    common = paths[0][0]
    for k in range(min(len(p) for p in paths)):
        if all(p[k] is paths[0][k] for p in paths):
            common = paths[0][k]
        else:
            break
    return common


def copy_tree(root: Node, parent: Node | None = None) -> Node:
    new = Node(name=root.name, length=root.length, parent=parent)
    new.children = [copy_tree(c, new) for c in root.children]
    return new


def merge_trees(pan: Node, vert: Node, anchors: tuple[str, str] = VERTEBRATA_ANCHORS) -> Node:
    """The pan-taxonomic tree with its vertebrate clade replaced by the vertebrate tree's Vertebrata.

    The vertebrate tree also carries outgroups (yeast, fly, worm, the two Ciona). One that the pan
    tree lacks (Ciona intestinalis) is placed where the vertebrate tree places it: as sister to the
    pan tree's clade of the tips that share its parent in the vertebrate tree, and only when that
    clade in the pan tree holds nothing else; otherwise the merge raises, because placing it would be
    a guess this module does not make. Lengths are dropped: the two files have none in common, and
    branch lengths are chosen by `grafen_lengths` or `equal_lengths` on the merged topology.
    """
    pan, vert = copy_tree(pan), copy_tree(vert)
    vtips = {norm(t.name): t for t in leaves(vert)}
    missing = [a for a in anchors if a not in vtips]
    if missing:
        raise ValueError(f"the vertebrate tree has no tip for {missing}")
    v_vert = mrca([vtips[a] for a in anchors])
    vert_names = {norm(t.name) for t in leaves(v_vert)}
    pan_tips = leaves(pan)
    inside = [t for t in pan_tips if norm(t.name) in vert_names]
    if not inside:
        raise ValueError("the pan tree holds no vertebrate of the vertebrate tree")
    p_vert = mrca(inside)
    stray = sorted(norm(t.name) for t in leaves(p_vert) if norm(t.name) not in vert_names)
    if stray:
        raise ValueError(f"the pan tree's vertebrate clade holds non-vertebrates {stray}")
    outgroups = {norm(t.name) for t in leaves(vert)} - vert_names
    pan_names = {norm(t.name) for t in pan_tips}
    for name in sorted(outgroups - pan_names):
        tip = vtips[name]
        clade = {norm(t.name) for t in leaves(tip.parent)} - {name} if tip.parent else set()
        here = [t for t in leaves(pan) if norm(t.name) in clade]
        if not here:
            raise ValueError(f"vertebrate-tree outgroup {name} has no relative in the pan tree")
        anchor = mrca(here)
        if {norm(t.name) for t in leaves(anchor)} - clade:
            raise ValueError(f"vertebrate-tree outgroup {name}: its clade is not a clade of the pan tree")
        joint = Node(name=tip.parent.name if tip.parent else "", parent=anchor.parent)
        if anchor.parent is not None:
            sib = anchor.parent.children
            sib[sib.index(anchor)] = joint
        joint.children = [anchor, Node(name=tip.name, parent=joint)]
        anchor.parent = joint
        if joint.parent is None:
            pan = joint
    graft = copy_tree(v_vert)
    if p_vert.parent is None:
        out = graft
    else:
        siblings = p_vert.parent.children
        siblings[siblings.index(p_vert)] = graft
        graft.parent = p_vert.parent
        out = pan
    for v in preorder(out):
        v.length = None
    return out


def match_species(species: list[str], tip_names: list[str]) -> dict[str, str]:
    """Each production name to exactly one tree tip (normalised), or raise naming what failed.

    In order: the name itself; the name with its assembly suffix (`_gca...`) stripped; the RENAMES
    table; a tip that extends the stripped name with a strain (`gallus_gallus` ->
    `gallus_gallus_reference_breed`) or that the stripped name extends; and only when none does, a
    tip that extends the genus and species (`aedes_aegypti_lvpagwg` ->
    `aedes_aegypti_strain_lvp_agwg`). Two candidates is an error, never a choice, and so is two
    species on one tip.
    """
    tips = [norm(t) for t in tip_names]
    tipset = set(tips)
    out: dict[str, str] = {}
    problems = []
    for sp in species:
        base = re.sub(r"_gca_?\d.*$", "", sp)
        if sp in tipset:
            out[sp] = sp
            continue
        if base in tipset:
            out[sp] = base
            continue
        if base in RENAMES and RENAMES[base] in tipset:
            out[sp] = RENAMES[base]
            continue
        cands = [t for t in tips if t.startswith(base + "_") or base.startswith(t + "_")]
        if not cands:  # a strain the tree spells differently: genus and species, then anything
            binomial = "_".join(base.split("_")[:2])
            cands = [t for t in tips if t.startswith(binomial + "_")]
        if len(cands) == 1:
            out[sp] = cands[0]
        else:
            problems.append((sp, cands))
    if problems:
        raise ValueError(f"species without exactly one tree tip: {problems[:10]} ({len(problems)})")
    seen: dict[str, str] = {}
    for sp, t in out.items():
        if t in seen:
            raise ValueError(f"two species on one tip {t}: {seen[t]} and {sp}")
        seen[t] = sp
    return out


def prune(root: Node, keep: set[str]) -> Node:
    """A copy holding only the tips whose normalised name is in `keep`, unary nodes collapsed."""

    def rec(v: Node) -> Node | None:
        if not v.children:
            return Node(name=norm(v.name)) if norm(v.name) in keep else None
        kids = [k for k in (rec(c) for c in v.children) if k is not None]
        if not kids:
            return None
        if len(kids) == 1:
            return kids[0]
        new = Node(name=v.name, children=kids)
        for k in kids:
            k.parent = new
        return new

    out = rec(root)
    if out is None:
        raise ValueError("no tip kept")
    out.parent = None
    return out


def to_newick(root: Node) -> str:
    """Topology with tip names (internal names dropped); deterministic for a given tree."""

    def rec(v: Node) -> str:
        if not v.children:
            return v.name
        return "(" + ",".join(rec(c) for c in v.children) + ")"

    return rec(root) + ";"


def topology_digest(root: Node) -> str:
    return hashlib.sha256(to_newick(root).encode()).hexdigest()


@dataclass
class TreeArrays:
    """A tree as arrays, nodes in preorder (a subtree is a contiguous range), tips in column order."""

    names: list[str]  # node names, preorder
    parent: np.ndarray  # parent index per node, -1 at the root
    size: np.ndarray  # nodes in each node's subtree, itself included
    tip_node: np.ndarray  # node index of each column's tip
    under: np.ndarray  # nodes x columns bool: the column's tip is in the node's subtree
    columns: list[str]  # the tip name of each column


def arrays(root: Node, columns: list[str]) -> TreeArrays:
    """The tree's arrays with tip columns in the given order (the tips must be exactly `columns`)."""
    order = preorder(root)
    idx = {id(v): i for i, v in enumerate(order)}
    parent = np.array([idx[id(v.parent)] if v.parent is not None else -1 for v in order], dtype=np.int64)
    size = np.ones(len(order), dtype=np.int64)
    for i in range(len(order) - 1, 0, -1):
        size[parent[i]] += size[i]
    col = {c: k for k, c in enumerate(columns)}
    tips = [(i, v.name) for i, v in enumerate(order) if not v.children]
    if sorted(n for _, n in tips) != sorted(columns):
        raise ValueError("the tree's tips are not the columns")
    tip_node = np.empty(len(columns), dtype=np.int64)
    for i, name in tips:
        tip_node[col[name]] = i
    under = np.zeros((len(order), len(columns)), dtype=bool)
    for k, t in enumerate(tip_node):
        v = int(t)
        while v >= 0:
            under[v, k] = True
            v = int(parent[v])
    return TreeArrays([v.name for v in order], parent, size, tip_node, under, list(columns))


def grafen_lengths(t: TreeArrays) -> np.ndarray:
    """Grafen's branch lengths, rho 1: node height = descendant tips - 1, root at height 1."""
    tips = t.under.sum(axis=1).astype(float)
    height = (tips - 1) / (tips[0] - 1)
    out = np.zeros(len(t.parent))
    out[1:] = height[t.parent[1:]] - height[1:]
    return out


def newick_lengths(root: Node) -> np.ndarray:
    """The tree's own branch lengths in `arrays` order (preorder), 0 where a branch has none."""
    out = np.array([v.length or 0.0 for v in preorder(root)], dtype=float)
    out[0] = 0.0
    return out


def equal_lengths(t: TreeArrays) -> np.ndarray:
    out = np.ones(len(t.parent))
    out[0] = 0.0
    return out


def covariance(t: TreeArrays, lengths: np.ndarray) -> np.ndarray:
    """Brownian covariance between tips: the length of the path they share from the root."""
    u = t.under.astype(float)
    return u.T @ (lengths[:, None] * u)


# ------------------------------------------------------------------------------------------- data


def presence_bits(presence: dict) -> tuple[list[str], list[str], np.ndarray]:
    """(species, genes, genes x species bool) from `origin_presence_genome_wide`.

    Bit i from the left of a gene's hex is species[i], as `homology.presence_matrix` writes it and
    `profiling.unpack` reads it: the last n bits of the zero-filled binary string.
    """
    species: list[str] = presence["species"]
    n = len(species)
    genes = sorted(g for g, h in (presence.get("genes") or {}).items() if h)
    hexes = [presence["genes"][g] for g in genes]
    width = max(len(h) for h in hexes)
    nbytes = (width * 4 + 7) // 8
    buf = b"".join(int(h, 16).to_bytes(nbytes, "big") for h in hexes)
    arr = np.frombuffer(buf, dtype=np.uint8).reshape(len(hexes), nbytes)
    return species, genes, np.unpackbits(arr, axis=1)[:, -n:].astype(bool)


def pearson(a: np.ndarray, b: np.ndarray) -> float | None:
    a = a - a.mean()
    b = b - b.mean()
    den = math.sqrt(float(a @ a) * float(b @ b))
    return float(a @ b) / den if den else None


def bh(pvalues: list[float]) -> list[float]:
    """Benjamini-Hochberg q, monotone from the largest p down, in input order (as §14 computes it)."""
    m = len(pvalues)
    order = sorted(range(m), key=lambda i: pvalues[i])
    q = [1.0] * m
    prev = 1.0
    for k in range(m - 1, -1, -1):
        i = order[k]
        prev = min(prev, pvalues[i] * m / (k + 1))
        q[i] = prev
    return q


# ------------------------------------------------------------------------------- null A: pairs


class Brownian:
    """Two independent banks of Brownian draws along one tree, rows centred, for one-pass Pearson.

    A row of a bank is one realisation of a unit-rate Brownian motion at the tips (covariance C). Each
    tested pair uses bank one for its first library and bank two for its second, so the two
    deviations are independent; every pair reuses the same banks, which leaves each pair's p exact and
    makes the p values of different pairs dependent through the shared draws.
    """

    def __init__(self, cov: np.ndarray, sims: int = A_SIMS, seed: int = A_SEED):
        self.cov = cov
        self.n = cov.shape[0]
        chol = np.linalg.cholesky(cov)
        rng = np.random.default_rng(seed)
        banks = []
        for _ in range(2):
            w = rng.standard_normal((sims, self.n)) @ chol.T
            w -= w.mean(axis=1, keepdims=True)
            banks.append(w)
        self.w1, self.w2 = banks
        self.n11 = np.einsum("ij,ij->i", self.w1, self.w1)
        self.n22 = np.einsum("ij,ij->i", self.w2, self.w2)
        self.n12 = np.einsum("ij,ij->i", self.w1, self.w2)
        self.inv = np.linalg.inv(cov)
        self.inv1 = self.inv @ np.ones(self.n)
        self.one_inv_one = float(self.inv1.sum())

    def rate(self, d: np.ndarray) -> tuple[float, float]:
        """GLS root mean and REML Brownian rate of a trait at the tips."""
        mu = float(self.inv1 @ d) / self.one_inv_one
        r = d - mu
        return mu, float(r @ self.inv @ r) / (self.n - 1)

    def phylo_correlation(self, da: np.ndarray, db: np.ndarray) -> float | None:
        """The correlation of two traits once the tree's covariance is removed (GLS; descriptive)."""
        ra = da - float(self.inv1 @ da) / self.one_inv_one
        rb = db - float(self.inv1 @ db) / self.one_inv_one
        den = math.sqrt(float(ra @ self.inv @ ra) * float(rb @ self.inv @ rb))
        return float(ra @ self.inv @ rb) / den if den else None

    def null_r(self, fa: np.ndarray, fb: np.ndarray, sa: float, sb: float) -> np.ndarray:
        """Pearson r of (fa + sa w1) and (fb + sb w2) for every row: the statistic under the null."""
        a = fa - fa.mean()
        b = fb - fb.mean()
        w1a, w1b = self.w1 @ a, self.w1 @ b
        w2a, w2b = self.w2 @ a, self.w2 @ b
        num = float(a @ b) + sb * w2a + sa * w1b + sa * sb * self.n12
        va = float(a @ a) + 2 * sa * w1a + sa * sa * self.n11
        vb = float(b @ b) + 2 * sb * w2b + sb * sb * self.n22
        return num / np.sqrt(va * vb)


@dataclass
class PairInputs:
    """What null A needs from the data: profiles per species, origins, the genome's profile."""

    species: list[str]
    genes: list[str]
    bits: np.ndarray  # genes x species
    row: dict[str, int]
    origin_of_row: np.ndarray  # stratum name per gene row ('' when the gene has no origin)
    genome: np.ndarray
    stratum_mean: dict[str, np.ndarray]
    pools: dict[str, np.ndarray]


def pair_inputs(presence: dict, origin: dict) -> PairInputs:
    species, genes, bits = presence_bits(presence)
    row = {g: i for i, g in enumerate(genes)}
    origins = {g: v.get("origin") for g, v in (origin.get("genes") or {}).items()}
    origin_of_row = np.array([origins.get(g) or "" for g in genes], dtype=object)
    f = bits.astype(float)
    pools: dict[str, np.ndarray] = {}
    means: dict[str, np.ndarray] = {}
    for s in sorted({o for o in origin_of_row if o}):
        rows = np.flatnonzero(origin_of_row == s)
        pools[s] = rows
        means[s] = f[rows].mean(axis=0)
    return PairInputs(species, genes, bits, row, origin_of_row, f.mean(axis=0), means, pools)


def profile(d: PairInputs, rows: np.ndarray) -> np.ndarray:
    return d.bits[rows].mean(axis=0) if len(rows) else np.zeros(len(d.species))


def age_expectation(d: PairInputs, rows: np.ndarray) -> np.ndarray:
    """The exact mean profile of §14's matched null: a stratum-for-stratum draw of the same make-up."""
    strata = d.origin_of_row[rows]
    strata = strata[strata != ""]
    if not len(strata):
        return np.zeros(len(d.species))
    names, counts = np.unique(strata.astype(str), return_counts=True)
    total = np.zeros(len(d.species))
    for s, c in zip(names, counts, strict=True):
        total += c * d.stratum_mean[str(s)]
    return total / counts.sum()


def matched_rows(d: PairInputs, rows: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """A random gene set of the same origin-stratum make-up (§14's draw, without replacement per stratum)."""
    strata = d.origin_of_row[rows]
    out = []
    for s in sorted(set(strata.tolist()) - {""}):
        k = int((strata == s).sum())
        pool = d.pools[s]
        out.append(rng.choice(pool, size=min(k, len(pool)), replace=False))
    return np.sort(np.concatenate(out)) if out else np.array([], dtype=np.int64)


def pair_test(d: PairInputs, bm: Brownian, rows_a: np.ndarray, rows_b: np.ndarray) -> dict[str, Any]:
    """§14's statistic for one pair, and where it falls under the Brownian deviation null."""
    pa, pb = profile(d, rows_a), profile(d, rows_b)
    ea, eb = age_expectation(d, rows_a), age_expectation(d, rows_b)
    r_obs = pearson(pa - d.genome, pb - d.genome)
    if r_obs is None:
        return {"r_exclusive": None}
    da, db = pa - ea, pb - eb
    _, va = bm.rate(da)
    _, vb = bm.rate(db)
    null = bm.null_r(ea - d.genome, eb - d.genome, math.sqrt(va), math.sqrt(vb))
    exceed = int((null >= r_obs).sum())
    return {
        "r_exclusive": r_obs,
        "r_age_expected": pearson(ea - d.genome, eb - d.genome),
        "brownian_rate": [va, vb],
        "tree_null_mean": float(null.mean()),
        "tree_null_sd": float(null.std(ddof=1)),
        "tree_draws_at_or_above": exceed,
        "p_tree": (exceed + 1) / (len(null) + 1),
        "deviation_phylo_correlation": bm.phylo_correlation(da, db),
    }


def library_sets(
    members: dict[str, list[str]], d: PairInputs, min_members: int = MIN_MEMBERS
) -> dict[str, set[str]]:
    sets = {lib: {s for s in syms if s in d.row} for lib, syms in members.items()}
    return {lib: s for lib, s in sorted(sets.items()) if len(s) >= min_members}


def testable_pairs(sets: dict[str, set[str]], min_members: int = MIN_MEMBERS) -> list[tuple]:
    """Every pair §14 tested: at least `min_members` members exclusive to each side."""
    out = []
    for a, b in combinations(sorted(sets), 2):
        only_a, only_b = sets[a] - sets[b], sets[b] - sets[a]
        if len(only_a) >= min_members and len(only_b) >= min_members:
            out.append((a, b, sorted(only_a), sorted(only_b)))
    return out


def rows_of(d: PairInputs, symbols) -> np.ndarray:
    return np.array(sorted(d.row[s] for s in symbols), dtype=np.int64)


def null_a(
    d: PairInputs,
    bm: Brownian,
    sets: dict[str, set[str]],
    calibrate: bool = True,
    halves: bool = True,
    progress=None,
) -> dict[str, Any]:
    """Every testable pair under the tree null, plus the calibration pairs and the library halves."""
    pairs = testable_pairs(sets)
    rows: list[dict[str, Any]] = []
    for k, (a, b, oa, ob) in enumerate(pairs):
        row: dict[str, Any] = {"libraries": [a, b], "exclusive_members": [len(oa), len(ob)]}
        row.update(pair_test(d, bm, rows_of(d, oa), rows_of(d, ob)))
        rows.append(row)
        if progress and k % 100 == 0:
            progress("pairs", k, len(pairs))
    tested = [r for r in rows if r.get("p_tree") is not None]
    for r, q in zip(tested, bh([r["p_tree"] for r in tested]), strict=True):
        r["q_tree"] = q
    out: dict[str, Any] = {"pairs": rows, "tested": len(tested)}
    if calibrate:
        rng = np.random.default_rng(A_CAL_SEED)
        cal = []
        for k, (a, b, oa, ob) in enumerate(pairs):
            ra = matched_rows(d, rows_of(d, oa), rng)
            rb = matched_rows(d, rows_of(d, ob), rng)
            t = pair_test(d, bm, ra, rb)
            cal.append({"for": [a, b], "p_tree": t.get("p_tree"), "r": t.get("r_exclusive")})
            if progress and k % 100 == 0:
                progress("calibration", k, len(pairs))
        ps = [c["p_tree"] for c in cal if c["p_tree"] is not None]
        hits = sum(1 for p in ps if p <= 0.05)
        out["calibration"] = {
            "pairs": len(ps),
            "at_or_below_005": hits,
            "share_at_or_below_005": hits / len(ps) if ps else None,
            "q_at_or_below_005": sum(1 for q in bh(ps) if q <= Q_ALPHA) if ps else 0,
        }
    if halves:
        rng = np.random.default_rng(A_HALVES_SEED)
        hv = []
        for lib in sorted(sets):
            r = rows_of(d, sets[lib])
            perm = rng.permutation(r)
            h1, h2 = np.sort(perm[: len(r) // 2]), np.sort(perm[len(r) // 2 :])
            t = pair_test(d, bm, h1, h2)
            hv.append({"library": lib, "p_tree": t.get("p_tree"), "r": t.get("r_exclusive")})
        ps = [h["p_tree"] for h in hv if h["p_tree"] is not None]
        hits = sum(1 for p in ps if p <= 0.05)
        out["halves"] = {
            "libraries": len(ps),
            "at_or_below_005": hits,
            "share_at_or_below_005": hits / len(ps) if ps else None,
            "rows": hv,
        }
    return out


# ---------------------------------------------------------------------------- null B: coherence


class Dollo:
    """Each gene born once on the tree and lost independently below, at its own rate.

    Built from the real genes' presence (genes x species, human excluded) and a tree whose columns are
    the species plus human. The birth node is the most recent common ancestor of human and the species
    carrying the gene; a node holds the gene when it is under the birth node and has a carrier (human
    included) below it; a loss is a branch whose parent holds the gene and whose child does not.
    """

    def __init__(self, t: TreeArrays, bits: np.ndarray, human_column: int):
        self.t = t
        n_nodes, n_cols = t.under.shape
        self.human_column = human_column
        self.species_columns = np.array([k for k in range(n_cols) if k != human_column], dtype=np.int64)
        if bits.shape[1] != len(self.species_columns):
            raise ValueError("presence columns are not the tree's species columns")
        path = []
        v = int(t.tip_node[human_column])
        while v >= 0:
            path.append(v)
            v = int(t.parent[v])
        self.path = np.array(path[::-1], dtype=np.int64)  # root first
        self.on_path = np.zeros(n_nodes, dtype=bool)
        self.on_path[self.path] = True
        g = np.zeros((bits.shape[0], n_cols), dtype=np.float32)
        g[:, self.species_columns] = bits
        g[:, human_column] = 1.0
        carriers = g.sum(axis=1)
        below = g @ t.under.astype(np.float32).T  # genes x nodes: carriers below each node
        holds_all = below[:, self.path] >= carriers[:, None] - 0.5
        # the deepest path node still holding every carrier; the root always does
        self.birth = (holds_all.shape[1] - 1 - np.argmax(holds_all[:, ::-1], axis=1)).astype(np.int64)
        self.carriers = bits.sum(axis=1).astype(np.int64)
        present = below > 0.5
        for j, u in enumerate(self.path):
            rows = np.flatnonzero(self.birth == j)
            if len(rows):
                outside = np.ones(n_nodes, dtype=bool)
                outside[int(u) : int(u + t.size[u])] = False
                present[np.ix_(rows, np.flatnonzero(outside))] = False
        par = t.parent.copy()
        par[0] = 0
        at_risk = present[:, par]
        at_risk[:, 0] = False
        lost = at_risk & ~present
        self.losses = lost.sum(axis=0)
        self.at_risk = at_risk.sum(axis=0)
        p = (self.losses + LOSS_PSEUDO[0]) / (self.at_risk + LOSS_PSEUDO[1])
        p[self.on_path] = 0.0
        self.p = p
        self.theta = self._fit_theta()

    def edges_below(self, j: int) -> np.ndarray:
        """Off-path nodes strictly under path node j: the branches a gene born there can be lost on."""
        u = int(self.path[j])
        nodes = np.arange(u + 1, u + int(self.t.size[u]))
        return nodes[~self.on_path[nodes]]

    def log_survival(self, j: int) -> tuple[np.ndarray, np.ndarray]:
        """For a gene born at path node j: the species columns under it, and the sum of log(1 - p) to each."""
        u = int(self.path[j])
        cols = self.species_columns[self.t.under[u, self.species_columns]]
        e = self.edges_below(j)
        logq = np.log1p(-self.p[e])
        ls = self.t.under[np.ix_(e, cols)].astype(float).T @ logq
        return cols, ls

    def _fit_theta(self) -> np.ndarray:
        theta = np.zeros(len(self.birth))
        for j in range(len(self.path)):
            rows = np.flatnonzero(self.birth == j)
            if not len(rows):
                continue
            cols, ls = self.log_survival(j)
            if not len(cols):
                continue
            k = self.carriers[rows].astype(float)
            lo = np.full(len(rows), -20.0)  # bisection on log theta
            hi = np.full(len(rows), 12.0)
            for _ in range(80):
                mid = (lo + hi) / 2
                f = np.exp(np.exp(mid)[:, None] * ls[None, :]).sum(axis=1)
                too_many = f > k
                lo = np.where(too_many, mid, lo)
                hi = np.where(too_many, hi, mid)
            th = np.exp((lo + hi) / 2)
            th[k >= len(cols)] = 0.0
            theta[rows] = th
        return theta

    def expected_carriers(self, row: int) -> float:
        _, ls = self.log_survival(int(self.birth[row]))
        return float(np.exp(self.theta[row] * ls).sum())

    def simulate(self, draws: int, seed: int, chunk: int = 400) -> np.ndarray:
        """genes x draws x species bool: every gene's presence drawn `draws` times along the tree."""
        rng = np.random.default_rng(seed)
        n_sp = len(self.species_columns)
        out = np.zeros((len(self.birth), draws, n_sp), dtype=bool)
        col_of = {int(c): k for k, c in enumerate(self.species_columns)}
        for j in range(len(self.path)):
            rows = np.flatnonzero(self.birth == j)
            if not len(rows):
                continue
            u = int(self.path[j])
            cols = self.species_columns[self.t.under[u, self.species_columns]]
            if not len(cols):
                continue
            e = self.edges_below(j)
            inc = self.t.under[np.ix_(e, cols)].astype(np.float32)
            logq = np.log1p(-self.p[e])
            local = np.array([col_of[int(c)] for c in cols], dtype=np.int64)
            for s in range(0, len(rows), chunk):
                r = rows[s : s + chunk]
                surv = np.exp(self.theta[r][:, None] * logq[None, :])  # genes x edges
                dead = rng.random((len(r), draws, len(e))) >= surv[:, None, :]
                hits = dead.reshape(-1, len(e)).astype(np.float32) @ inc
                block = np.zeros((len(r), draws, n_sp), dtype=bool)
                block[:, :, local] = (hits < 0.5).reshape(len(r), draws, len(cols))
                out[r] = block
        return out


def mean_pair_jaccard(x: np.ndarray) -> float | None:
    """Exact mean Jaccard over all pairs of rows of a bool matrix (0 for a pair with an empty union)."""
    n = x.shape[0]
    if n < 2:
        return None
    f = x.astype(np.float32)
    inter = f @ f.T
    c = f.sum(axis=1)
    union = c[:, None] + c[None, :] - inter
    j = np.divide(inter, union, out=np.zeros_like(inter), where=union > 0)
    total = float(j.sum(dtype=np.float64)) - float(np.trace(j, dtype=np.float64))
    return total / (n * (n - 1))


def tree_excess(rows: np.ndarray, bits: np.ndarray, sims: np.ndarray) -> tuple[float, float, list[float]]:
    """(observed mean pair Jaccard, its mean over the model's draws, the draws) for one gene set."""
    obs = mean_pair_jaccard(bits[rows])
    draws = [mean_pair_jaccard(sims[rows, d]) for d in range(sims.shape[1])]
    return float(obs), float(np.mean(draws)), [float(x) for x in draws]


def birth_matched(dollo: Dollo, rows: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """A random gene set matched member for member on birth node, drawn without replacement per node."""
    out = []
    births = dollo.birth[rows]
    for j in sorted(set(births.tolist())):
        k = int((births == j).sum())
        pool = np.flatnonzero(dollo.birth == j)
        out.append(rng.choice(pool, size=min(k, len(pool)), replace=False))
    return np.sort(np.concatenate(out)) if out else np.array([], dtype=np.int64)


def classify(at_or_above: int, sets: int = B_SETS) -> str:
    """§19's three classes, by what they meant: clear, marginal, at or below.

    Clear is an empirical one-sided p, (k + 1) / (sets + 1), at or below 0.05; at or below is half or
    more of the random sets at or above the library (the library at or below their median). With
    §19's 20 draws these are exactly §19's "0 draws at or above" and "10 or more"; with 200 they are
    k <= 9 and k >= 100.
    """
    if (at_or_above + 1) / (sets + 1) <= B_ALPHA:
        return "clear"
    if 2 * at_or_above >= sets:
        return "at_or_below"
    return "marginal"


def chance_ceiling(n: int, alpha: float = B_ALPHA, level: float = 0.95) -> int:
    """The smallest k such that n independent libraries with no effect give k or fewer clear
    `level` of the time."""
    cdf = 0.0
    for k in range(n + 1):
        cdf += math.comb(n, k) * alpha**k * (1 - alpha) ** (n - k)
        if cdf >= level:
            return k
    return n


def null_b(
    dollo: Dollo,
    bits: np.ndarray,
    sims: np.ndarray,
    library_rows: dict[str, np.ndarray],
    sets: int = B_SETS,
    seed: int = B_SEED,
    progress=None,
) -> list[dict[str, Any]]:
    """Each library's excess over its own tree model, against birth-matched real gene sets' excess."""
    rng = np.random.default_rng(seed)
    out = []
    for k, (lib, rows) in enumerate(sorted(library_rows.items())):
        obs, model, _ = tree_excess(rows, bits, sims)
        t_lib = obs - model
        null_t, misfit = [], 0
        for _ in range(sets):
            r = birth_matched(dollo, rows, rng)
            o, m, draws = tree_excess(r, bits, sims)
            null_t.append(o - m)
            misfit += int(all(o > x for x in draws))
        above = sum(1 for x in null_t if x >= t_lib)
        out.append(
            {
                "library": lib,
                "members_with_a_profile": len(rows),
                "mean_pair_jaccard": obs,
                "tree_model_mean": model,
                "excess_over_tree_model": t_lib,
                "random_sets_excess_mean": float(np.mean(null_t)),
                "random_sets_excess_sd": float(np.std(null_t, ddof=1)),
                "random_sets_at_or_above": above,
                "class": classify(above, sets),
                "random_sets_above_every_model_draw": misfit,
            }
        )
        if progress:
            progress("libraries", k, len(library_rows))
    return out


# ----------------------------------------------------------------------------------------- fetch


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fetch_trees(tree_dir: Path = TREE_DIR) -> dict[str, str]:
    """Download the two species trees once (about 13 KB together) and check them against the pins."""
    tree_dir.mkdir(parents=True, exist_ok=True)
    out = {}
    for url, name in ((PAN_TREE_URL, PAN_TREE.name), (VERT_TREE_URL, VERT_TREE.name)):
        path = tree_dir / name
        if not path.exists():
            with urllib.request.urlopen(url, timeout=60) as fh:  # noqa: S310 - a pinned https URL
                path.write_bytes(fh.read())
        digest = sha256(path)
        if digest != TREE_SHA256[name]:
            raise ValueError(f"{path} hashes to {digest}, not the pinned {TREE_SHA256[name]}")
        out[name] = digest
    return out


def species_tree(species: list[str], with_human: bool, tree_dir: Path = TREE_DIR) -> tuple[Node, dict]:
    """The merged tree pruned to the presence species (and human), tips renamed to the species names."""
    pan = parse_newick((tree_dir / PAN_TREE.name).read_text())
    vert = parse_newick((tree_dir / VERT_TREE.name).read_text())
    merged = merge_trees(pan, vert)
    mapping = match_species(species, [t.name for t in leaves(merged)])
    keep = set(mapping.values()) | ({HUMAN} if with_human else set())
    pruned = prune(merged, keep)
    back = {v: k for k, v in mapping.items()}
    for t in leaves(pruned):
        t.name = back.get(t.name, t.name)
    return pruned, mapping


# -------------------------------------------------------------------------------------- outcomes
# The exact wording of every outcome, fixed before the run. The run fills the numbers into the
# sentence its decision rule selects and changes no other word. Each states its denominators.

A_OUTCOMES = {
    "not_reproduced": "No reading: §14's observed statistic did not reproduce. {n_mismatch} of "
    "{n_tested} testable library pairs differ from the stored r_exclusive by more than 0.0001, or the "
    "testable pairs are not §14's {n_tested_14}, so the pairs this null would judge are not the pairs "
    "§14 judged.",
    "null_not_calibrated": "No reading: the species-tree null is not calibrated on this data. Of "
    "{n_cal} pairs of random gene sets drawn to the tested pairs' origin-stratum make-up, {cal_hits} "
    "({cal_share}) have a tree-null p at or below 0.05, above the registered ceiling of 10%, so the "
    "null finds random genes alike more often than it should and a library pair's p under it is not "
    "evidence. §14's {n_14} of {n_tested} pairs stand as written, with their caveat.",
    "none_survive": "None of the {n_tested} testable library pairs is more alike than age and shared "
    "ancestry predict: under a Brownian null along the Ensembl Compara species tree (Grafen branch "
    "lengths), no pair has a Benjamini-Hochberg q at or below 0.05, and that includes all {n_14} "
    "pairs §14 called more alike than age at the same false discovery rate. The null is calibrated "
    "({cal_hits} of {n_cal} random age-matched pairs at p at or below 0.05) and sees {halves_hits} "
    "of {n_halves} libraries against their own random halves at p at or below 0.05. On this model "
    "the {n_14} pairs are what age and shared ancestry give, and no library pair is called a "
    "dependency.",
    "none_survive_low_power": "None of the {n_tested} testable library pairs is more alike than age "
    "and shared ancestry predict under a Brownian null along the Ensembl Compara species tree "
    "(Grafen branch lengths): no pair has a Benjamini-Hochberg q at or below 0.05, including the "
    "{n_14} pairs §14 called more alike than age. This is weak evidence against them: the null is "
    "calibrated ({cal_hits} of {n_cal} random age-matched pairs at p at or below 0.05) but sees only "
    "{halves_hits} of {n_halves} libraries against their own random halves at p at or below 0.05, "
    "below the registered half, so on this tree the test has little power and a pair's failure does "
    "not show that it is inertia. No library pair is called a dependency.",
    "some_survive": "{n_surv} of {n_tested} testable library pairs are more alike than age and "
    "shared ancestry predict under a Brownian null along the Ensembl Compara species tree (Grafen "
    "branch lengths; Benjamini-Hochberg q at or below 0.05): {k_14} of the {n_14} pairs §14 called "
    "more alike than age, and {k_new} that §14 did not call. The other {lost_14} of §14's {n_14} do "
    "not survive the tree. The null is calibrated ({cal_hits} of {n_cal} random age-matched pairs at "
    "p at or below 0.05) and sees {halves_hits} of {n_halves} libraries against their own random "
    "halves at p at or below 0.05. The surviving pairs meet, on this model and these branch lengths, "
    "the condition the area set before a pair could be called a dependency; agreement of presence "
    "along the tree does not show that one library acts on the other, and an orthologue's presence "
    "is not function.",
}
A_SENSITIVITY = (
    " Under equal branch lengths instead of Grafen's, {n_surv_eq} of {n_tested} pairs have q at or "
    "below 0.05, {k_14_eq} of them among §14's {n_14}, and {cal_hits_eq} of {n_cal} random "
    "age-matched pairs have p at or below 0.05; the registered reading is the Grafen one."
)
A_RULE = {
    "1_reproduce": "if the testable pairs are not §14's (same pairs, same count) or any pair's "
    "recomputed r_exclusive differs from the stored one by more than 0.0001, select not_reproduced",
    "2_calibration": "else if more than 10% of the calibration pairs (one random pair of gene sets "
    "per testable pair, each set drawn to its side's origin-stratum make-up) have p_tree at or below "
    "0.05, select null_not_calibrated",
    "3_none": "else if no library pair has q_tree at or below 0.05: select none_survive_low_power "
    "when fewer than half of the libraries' random halves have p_tree at or below 0.05, otherwise "
    "none_survive",
    "4_some": "else select some_survive",
    "ties": "'at or below' includes equality; 'more than 10%' and 'fewer than half' are strict",
    "always": "the equal-length sensitivity sentence is appended to whichever outcome is selected, "
    "except not_reproduced, and decides nothing",
}

B_OUTCOMES = {
    "not_reproduced": "No reading: §19's observed statistic did not reproduce. {n_mismatch} of "
    "{n_lib} libraries differ from §19's member count or from its mean pair Jaccard by more than "
    "0.005, so the libraries this null would judge are not the libraries §19 judged.",
    "none_clear": "No library's members co-occur beyond what the species tree and real genes "
    "predict. Each library's mean pair Jaccard, less what a stochastic Dollo model of every member's "
    "own birth and losses along the Ensembl Compara species tree predicts, was set against the same "
    "quantity for {n_sets} random gene sets matched member for member on birth node: 0 of {n_lib} "
    "libraries are clear (empirical one-sided p at or below 0.05), {n_marg} are marginal and "
    "{n_below} sit at or below the random sets' median. None of the {n_19} libraries §19 called "
    "clear of prevalence-matched genes stays clear. On this null their coherence is what birth on "
    "the tree and loss rate give; nothing is called a dependency.",
    "at_chance": "{n_clear} of {n_lib} libraries' members co-occur beyond what the species tree and "
    "real genes predict, which is no more than chance gives. Each library's mean pair Jaccard, less "
    "what a stochastic Dollo model of every member's own birth and losses along the Ensembl Compara "
    "species tree predicts, was set against the same quantity for {n_sets} random gene sets matched "
    "member for member on birth node: {n_clear} are clear (empirical one-sided p at or below 0.05), "
    "{n_marg} marginal and {n_below} at or below the random sets' median. With no effect anywhere "
    "about {expected_chance} of {n_lib} would be clear, and {ceiling} or fewer 95 times in 100 if the "
    "libraries were independent. {k_19} of the {n_19} libraries §19 called clear are among the "
    "{n_clear}. No library's coherence is read as beyond the tree; nothing is called a dependency.",
    "some_clear": "{n_clear} of {n_lib} libraries' members co-occur beyond what the species tree and "
    "real genes predict. Each library's mean pair Jaccard, less what a stochastic Dollo model of "
    "every member's own birth and losses along the Ensembl Compara species tree predicts, was set "
    "against the same quantity for {n_sets} random gene sets matched member for member on birth "
    "node: {n_clear} are clear (empirical one-sided p at or below 0.05), {n_marg} marginal and "
    "{n_below} at or below the random sets' median; with no effect anywhere about {expected_chance} "
    "would be clear, and more than {ceiling} would arise by chance at most 5 times in 100 if the "
    "libraries were independent. {k_19} of the {n_19} libraries §19 called clear stay clear, and "
    "{k_new} that §19 "
    "did not call clear are clear here. This is agreement of presence across species beyond birth "
    "clade and loss rate; it names no pair, calls nothing a dependency, and an orthologue's presence "
    "is not function.",
}
B_RULE = {
    "1_reproduce": "if any library's member count differs from §19's members_with_a_profile, or its "
    f"exact mean pair Jaccard differs from §19's stored value by more than {B_TOLERANCE}, or the "
    "libraries are not §19's 42, select not_reproduced",
    "2_none": "else if no library is clear, select none_clear",
    "3_chance": "else if the number of clear libraries is at or below chance_ceiling(n_lib) (the "
    "smallest k with binomial(n_lib, 0.05) CDF at or above 0.95), select at_chance",
    "4_some": "else select some_clear",
    "classes": "clear: (k + 1) / (B_SETS + 1) at or below 0.05, k the random sets whose excess over "
    "their own tree model is at or above the library's; at_or_below: 2k at or above B_SETS; "
    "marginal otherwise",
}


def select_a(
    grafen: dict[str, Any],
    equal: dict[str, Any] | None,
    stored_14: dict[str, Any],
) -> tuple[str, str, dict[str, Any]]:
    """Apply A_RULE to null A's output; return (outcome, filled wording, the figures used)."""
    stored = {tuple(r["libraries"]): r for r in stored_14["pairs"] if r.get("r_exclusive") is not None}
    called = {tuple(r["libraries"]) for r in stored_14["significant_q05"]}
    rows = [r for r in grafen["pairs"] if r.get("p_tree") is not None]
    mism = sum(
        1
        for r in rows
        if tuple(r["libraries"]) not in stored
        or abs(r["r_exclusive"] - stored[tuple(r["libraries"])]["r_exclusive"]) > 1e-4
    )
    same_pairs = {tuple(r["libraries"]) for r in rows} == set(stored)
    surv = {tuple(r["libraries"]) for r in rows if r["q_tree"] <= Q_ALPHA}
    cal = grafen.get("calibration") or {}
    hv = grafen.get("halves") or {}
    f: dict[str, Any] = {
        "n_tested": len(rows),
        "n_tested_14": len(stored),
        "n_mismatch": mism,
        "n_14": len(called),
        "n_surv": len(surv),
        "k_14": len(surv & called),
        "k_new": len(surv - called),
        "lost_14": len(called - surv),
        "cal_hits": cal.get("at_or_below_005"),
        "n_cal": cal.get("pairs"),
        "cal_share": f"{cal.get('share_at_or_below_005', 0):.1%}" if cal.get("pairs") else None,
        "halves_hits": hv.get("at_or_below_005"),
        "n_halves": hv.get("libraries"),
    }
    if mism or not same_pairs:
        key = "not_reproduced"
    elif (cal.get("share_at_or_below_005") or 0) > CAL_MAX_SHARE:
        key = "null_not_calibrated"
    elif not surv:
        low = (hv.get("share_at_or_below_005") or 0) < HALVES_MIN_SHARE
        key = "none_survive_low_power" if low else "none_survive"
    else:
        key = "some_survive"
    text = A_OUTCOMES[key].format(**f)
    if equal is not None and key != "not_reproduced":
        rows_eq = [r for r in equal["pairs"] if r.get("p_tree") is not None]
        surv_eq = {tuple(r["libraries"]) for r in rows_eq if r["q_tree"] <= Q_ALPHA}
        f["n_surv_eq"] = len(surv_eq)
        f["k_14_eq"] = len(surv_eq & called)
        f["cal_hits_eq"] = (equal.get("calibration") or {}).get("at_or_below_005")
        text += A_SENSITIVITY.format(**f)
    f["survivors"] = sorted(list(p) for p in surv)
    f["survivors_among_14"] = sorted(list(p) for p in surv & called)
    return key, text, f


def select_b(rows: list[dict[str, Any]], stored_19: dict[str, Any]) -> tuple[str, str, dict[str, Any]]:
    """Apply B_RULE to null B's output; return (outcome, filled wording, the figures used)."""
    stored = {r["library"]: r for r in stored_19["libraries"]}
    called = {r["library"] for r in stored_19["libraries"] if r["matched_draws_at_or_above_observed"] == 0}
    mism = sum(
        1
        for r in rows
        if r["library"] not in stored
        or r["members_with_a_profile"] != stored[r["library"]]["members_with_a_profile"]
        or abs(r["mean_pair_jaccard"] - stored[r["library"]]["mean_pair_jaccard"]) > B_TOLERANCE
    )
    same = {r["library"] for r in rows} == set(stored)
    clear = {r["library"] for r in rows if r["class"] == "clear"}
    n = len(rows)
    f: dict[str, Any] = {
        "n_lib": n,
        "n_mismatch": mism,
        "n_sets": B_SETS,
        "n_clear": len(clear),
        "n_marg": sum(1 for r in rows if r["class"] == "marginal"),
        "n_below": sum(1 for r in rows if r["class"] == "at_or_below"),
        "n_19": len(called),
        "k_19": len(clear & called),
        "k_new": len(clear - called),
        "expected_chance": f"{n * B_ALPHA:.1f}",
        "ceiling": chance_ceiling(n),
    }
    if mism or not same:
        key = "not_reproduced"
    elif not clear:
        key = "none_clear"
    elif len(clear) <= f["ceiling"]:
        key = "at_chance"
    else:
        key = "some_clear"
    f["clear"] = sorted(clear)
    f["clear_among_19"] = sorted(clear & called)
    return key, B_OUTCOMES[key].format(**f), f
