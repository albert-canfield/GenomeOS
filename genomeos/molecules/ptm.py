# SPDX-License-Identifier: AGPL-3.0-or-later
"""Post-translational state: the modified residues UniProt records, and who writes them.

The protein model separates what a protein is (ProteinDefinition) from one protein in one place at
one time (ProteinState), and until now nothing populated the modifications of a state. UniProt's
"Modified residue", "Glycosylation", "Lipidation", "Disulfide bond" and "Cross-link" features carry
the site, the chemistry and, often, the writer ("Phosphoserine; by CDK5, PRPK, AMPK, NUAK1 and ATM").
This module reads them out of the compiled definitions: per protein, the sites by class with their
writers; genome-wide, the writer → substrate edges (a kinase's substrates, an acetyltransferase's) as
a curated layer of the knowledge graph, and the counts as a committed summary.

Curated evidence (UniProt, largely from the literature it cites); a modification is a site that can be
modified, not a measurement that it is modified in a given cell.
"""

from __future__ import annotations

import json
import re
from collections import Counter
from pathlib import Path
from typing import Any

CACHE = Path("data/knowledge/proteins")
INDEX = CACHE / "_ptm_index.json"
EVIDENCE = (
    "curated: UniProt sequence features (modified residue, glycosylation, lipidation, cross-link, disulfide)"
)
CLASSES = (
    ("phospho", ("phospho",)),
    ("acetyl", ("acetyl",)),
    ("methyl", ("methyl",)),
    ("ubiquitin", ("ubiquitin", "glycyl lysine isopeptide")),
    ("sumo", ("sumo",)),
    ("glyco", ("glycos", "glycan", "glcnac", "fucos", "mannos", "galactos", "xylos")),
    ("lipid", ("myristoyl", "palmitoyl", "farnesyl", "geranyl", "gpi-anchor", "lipoyl", "octanoyl")),
    ("disulfide", ("disulfide",)),
    ("hydroxy", ("hydroxy",)),
    ("nitro", ("nitro",)),
    ("adp-ribosyl", ("adp-ribosyl",)),
    ("citrulline", ("citrulline",)),
    ("sulfo", ("sulfo",)),
)
_BY = re.compile(r";\s*by\s+([^;.]+)")
_AUTO = re.compile(r"\bautocatal", re.I)


def classify(description: str, feature_type: str) -> str:
    d = description.lower()
    if feature_type == "Disulfide bond":
        return "disulfide"
    if feature_type == "Glycosylation":
        return "glyco"
    if feature_type == "Lipidation":
        return "lipid"
    for cls, keys in CLASSES:
        if any(k in d for k in keys):
            return cls
    return "other"


def writers_of(description: str) -> list[str]:
    """The enzymes named after "by": 'Phosphoserine; by CDK5, PRPK, AMPK, NUAK1 and ATM' → 5 symbols."""
    m = _BY.search(description)
    if not m:
        return []
    text = m.group(1)
    if _AUTO.search(text) or "autocatalysis" in text.lower():
        return ["(self)"]
    parts = re.split(r",|\band\b|/|\bor\b", text)
    out = []
    for p in parts:
        p = p.strip().strip(".")
        p = re.sub(r"\s*\(.*?\)", "", p)
        p = re.sub(r"^(isoform|in vitro|in )\b.*", "", p).strip()
        if p and len(p) <= 15 and re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9\-]*", p):
            out.append(p.upper())
    return out


def sites(defn: dict[str, Any]) -> list[dict[str, Any]]:
    """The modifiable sites of one compiled definition, with class and writers."""
    out = []
    for f in (defn["sections"].get("modifications") or {}).get("items") or []:
        desc = f.get("description") or ""
        out.append(
            {
                "type": f.get("type"),
                "start": f.get("start"),
                "end": f.get("end"),
                "class": classify(desc, f.get("type") or ""),
                "description": desc.split(";")[0].strip(),
                "writers": writers_of(desc),
            }
        )
    return out


def states(defn: dict[str, Any]) -> list[dict[str, Any]]:
    """ProteinState-shaped records: one per modifiable site, the writer as the condition."""
    pid = defn.get("id") or f"gene:{defn['gene']}"
    out = []
    for s in sites(defn):
        out.append(
            {
                "protein": pid,
                "modifications": [f"{s['description']} at {s['start']}"],
                "written_by": s["writers"],
                "class": s["class"],
                "evidence": EVIDENCE,
                "confidence": 0.8,
            }
        )
    return out


def build_index(cache_dir: Path = CACHE, out: Path | None = None) -> dict[str, Any]:
    """Genome-wide: sites per class, writers and their substrates, from every cached definition."""
    per_class: Counter = Counter()
    writers: dict[str, Counter] = {}
    proteins_with = 0
    proteins = 0
    site_total = 0
    per_protein: dict[str, dict[str, int]] = {}
    for p in sorted(cache_dir.glob("*.json")):
        if p.name.startswith("_"):
            continue
        try:
            d = json.loads(p.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        if "sections" not in d:
            continue
        proteins += 1
        ss = sites(d)
        if not ss:
            continue
        proteins_with += 1
        site_total += len(ss)
        counts = Counter(s["class"] for s in ss)
        per_class.update(counts)
        per_protein[d["gene"]] = dict(counts)
        for s in ss:
            for w in s["writers"]:
                if w != "(self)":
                    writers.setdefault(w, Counter())[d["gene"]] += 1
    index = {
        "proteins": proteins,
        "proteins_with_sites": proteins_with,
        "sites": site_total,
        "by_class": dict(per_class.most_common()),
        "writers": {w: dict(c) for w, c in writers.items()},
        "per_protein": per_protein,
        "evidence": EVIDENCE,
    }
    out = out or INDEX
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(index))
    return index


def load_index(path: Path = INDEX) -> dict[str, Any] | None:
    try:
        return json.loads(path.read_text()) if path.exists() else None
    except (OSError, json.JSONDecodeError):
        return None


def summary(index: dict[str, Any], top: int = 20) -> dict[str, Any]:
    """The committed summary: counts, and the writers with the most substrates."""
    writers = index.get("writers", {})
    ranked = sorted(writers.items(), key=lambda kv: -len(kv[1]))
    return {
        "proteins": index["proteins"],
        "proteins_with_sites": index["proteins_with_sites"],
        "sites": index["sites"],
        "by_class": index["by_class"],
        "writers": len(writers),
        "writer_edges": sum(len(c) for c in writers.values()),
        "top_writers": [
            {"writer": w, "substrates": len(c), "sites": sum(c.values())} for w, c in ranked[:top]
        ],
        "evidence": index["evidence"],
    }


# Observation, not occupancy -----------------------------------------------------------------------
#
# Ochoa et al. 2020 (Nat. Biotechnol. 38:365) reanalysed 6,801 public phosphoproteomics raw files from
# 104 cell types or tissues; the authors' funscoR package (LGPL) ships, per site, the number of distinct
# cell lines or tissues in which it was identified at 1% site-level FDR and its spectral count. Joined to
# the curated sites above, that turns "can be modified" into "identified phosphorylated in N cell types or
# tissues". It is not occupancy (the fraction of the protein modified), not a per-cell state, and a curated
# site missing from the reference is not evidence that it is never phosphorylated.

OBSERVATION = Path("data/results/phosphosite_observation.json")
OBSERVATION_MEANING = (
    "observed_in_cell_types_or_tissues: distinct cell lines or tissues whose public mass-spectrometry data "
    "identified this phosphosite (Ochoa et al. 2020, 1% site-level FDR); not occupancy, "
    "not a per-cell state, and absence is not evidence of no phosphorylation"
)
PHOSPHO_ACCEPTORS = frozenset("STY")


def accession_of(defn: dict[str, Any]) -> str | None:
    ident = ((defn.get("sections") or {}).get("identity") or {}).get("items") or {}
    acc = ident.get("accession")
    if acc:
        return acc
    pid = defn.get("id") or ""
    return pid.split(":", 1)[1] if pid.startswith("UniProt:") else None


def sequence_of(defn: dict[str, Any]) -> str:
    ident = ((defn.get("sections") or {}).get("identity") or {}).get("items") or {}
    return ident.get("sequence") or ""


def join_observations(
    reference: dict[tuple[str, int], tuple[str, int, int]],
    definitions: Any,
) -> dict[str, Any]:
    """Join curated phospho sites to a reference keyed by (accession, position).

    ``reference`` maps (accession, position) to (residue, cell types or tissues, spectral count).
    ``definitions`` yields compiled protein definitions. Only curated sites of class ``phospho`` are
    joined; a match counts only when the reference residue equals the residue at that position in the
    current UniProt sequence, and disagreements are counted, not joined.
    """
    curated = phospho = matched_position = disagree = 0
    other_class_on_reference = Counter()
    off_acceptor: list[str] = []
    joined: dict[str, list[list[Any]]] = {}
    for d in definitions:
        acc = accession_of(d)
        seq = sequence_of(d)
        for s in sites(d):
            curated += 1
            pos = s.get("start")
            if acc is None or not isinstance(pos, int):
                continue
            ref = reference.get((acc, pos))
            if s["class"] != "phospho":
                if ref is not None:
                    other_class_on_reference[s["class"]] += 1
                continue
            phospho += 1
            if ref is None:
                continue
            matched_position += 1
            here = seq[pos - 1] if 0 < pos <= len(seq) else ""
            if here != ref[0]:
                disagree += 1
                continue
            if here not in PHOSPHO_ACCEPTORS:
                off_acceptor.append(f"{acc}:{here}{pos}")
            joined.setdefault(acc, []).append([pos, here, ref[1], ref[2]])
    n_joined = sum(len(v) for v in joined.values())
    return {
        "curated_sites": curated,
        "curated_phospho_sites": phospho,
        "accession_position_matches": matched_position,
        "residue_disagreements": disagree,
        "joined": n_joined,
        "off_acceptor": off_acceptor,
        "other_class_on_reference_position": dict(other_class_on_reference.most_common()),
        "sites": {a: sorted(v) for a, v in sorted(joined.items())},
    }


def load_observation(path: Path = OBSERVATION) -> dict[str, Any] | None:
    return load_index(path)


def observed_in(accession: str, position: int, observation: dict[str, Any] | None) -> dict[str, Any] | None:
    """The observation record of one curated phosphosite, or None when the reference lacks it."""
    if not observation:
        return None
    for pos, residue, n, psm in observation.get("sites", {}).get(accession, []):
        if pos == position:
            return {
                "residue": residue,
                "observed_in_cell_types_or_tissues": n,
                "spectral_count": psm,
                "meaning": OBSERVATION_MEANING,
            }
    return None


# ---- relative abundance in tumours: CPTAC phosphoproteomics through cBioPortal (2026-09-28) -----------
#
# The CPTAC phosphoproteome profiles on cBioPortal name a site by gene symbol and a position on a RefSeq
# protein (sometimes with the accession, sometimes without), and the values are log2 abundance ratios to
# a pooled reference, per tumour. The functions below parse an entity into one site, carry a RefSeq
# position over to the UniProt canonical sequence of the compiled definition, and read the committed
# per-site summary. The registration the join is scored against is in docs/PROTEIN.md,
# "Relative abundance in tumours (2026-09-28)", and the constants that follow it.

CPTAC_RESULT = Path("data/results/phosphosite_tumour_abundance.json")
CPTAC_FIELD = "relative_abundance_in_tumours"
CPTAC_MEANING = (
    "relative_abundance_in_tumours: per CPTAC study on cBioPortal, the number of tumours with a value and "
    "the median of the per-tumour log2 ratio of this phosphosite's abundance to a pooled reference of "
    "tumours from the same study; not occupancy (the fraction of the protein modified), not comparable "
    "across studies as an absolute level, not per patient, and absence is not evidence of no phosphorylation"
)
# Registered 2026-09-28, after the mapping stage and before any curated site was joined or any value
# fetched (docs/PROTEIN.md, "Relative abundance in tumours (2026-09-28)"). Positions are UniProt
# canonical numbering, so EGFR's legacy Y1068 (mature-protein numbering) is Y1092 here.
CPTAC_REGISTRATION: dict[str, Any] = {
    "join": "UniProt accession + position + residue letter equal to a curated site of class phospho (S/T/Y)",
    "committed_tier": "refseq",  # entities that name their RefSeq protein; gene-keyed tiers are counted only
    "expected_refseq_keyed": (0.35, 0.55),  # share of the curated phospho sites gaining a value
    "expected_any_profile": (0.45, 0.65),  # counted, not committed
    "check_present_lung": ("P00533", 1092, "Y"),  # EGFR pY1068 (legacy), in lung adenocarcinoma
    "check_present_any": ("P06748", 125, "S"),  # NPM1 S125, the constitutive CK2 site
    "check_absent": ("P04637", 1),  # TP53 M1: not an acceptor
    "check_reported_only": ("P31749", 473, "S"),  # AKT1 S473: reported, no expectation (poor tryptic peptide)
    "refseq_keyed_mismatch_bound": 0.02,  # residue letter disagrees at the RefSeq end
    "gene_keyed_mismatch_bound": 0.02,  # a gene-keyed tier would be committed only under the same bound
    "median_log2_ratio_within": (-0.5, 0.5),  # pooled-reference ratios should centre near zero
}
CPTAC_LICENCE = (
    "Values: CPTAC data through cBioPortal (https://www.cbioportal.org), Open Database License 1.0 "
    "(https://opendatacommons.org/licenses/odbl/1-0/); this per-site summary is a derived database and is "
    "itself offered under ODbL 1.0 (share-alike), attribution: CPTAC and the cBioPortal study of each "
    "profile named in `studies`. Residue mapping: UniProt (CC BY 4.0, https://www.uniprot.org) and NCBI "
    "RefSeq (public, https://www.ncbi.nlm.nih.gov/refseq/). The repository keeps its own licences."
)
_SITE = re.compile(r"^([STY])(\d+)$")
_REFSEQ = re.compile(r"^(?:NP|XP|YP)_\d+\.\d+$")


def parse_cptac_entity(profile: str, entity: dict[str, Any]) -> dict[str, Any]:
    """One cBioPortal generic-assay entity of a CPTAC phospho profile, as a single site or a drop reason.

    Returns ``{"refseq": accession or None, "gene", "residue", "position"}`` for a single, localised
    phosphosite, else ``{"drop": reason}``. ``refseq`` is None when the profile names the gene only.
    """
    sid = entity.get("stableId") or ""
    props = entity.get("genericEntityMetaProperties") or {}
    gene = props.get("GENE_SYMBOL") or ""
    if "acetylprotein" in sid:
        return {"drop": "gene-level aggregate, not a site"}
    if profile.startswith("luad_cptac_2020"):
        m = re.match(r"^((?:NP|XP|YP)_\d+\.\d+)_(\d+)_(\d+)_(\d+)_(\d+)$", sid)
        if not m:
            return {"drop": "not a RefSeq protein"}
        acc, n, loc, a, b = m.group(1), *map(int, m.groups()[1:])
        name_sites = re.findall(r"_([STY])(\d+)[sty]", props.get("NAME") or "")
        if n != 1 or len(name_sites) != 1:
            return {"drop": "several sites in one entity"}
        if loc != 1 or a != b:
            return {"drop": "site not localised"}
        res, pos = name_sites[0][0], int(name_sites[0][1])
        if pos != a:
            return {"drop": "entity name and position disagree"}
        return {"refseq": acc, "gene": gene, "residue": res, "position": pos}
    if profile.startswith("brca_cptac_2020"):
        m = re.match(r"^((?:[STY]\d+[sty])+)_(\d+)_(\d+)_(\d+)_(\d+)$", props.get("PHOSPHOSITES") or "")
        if not m:
            return {"drop": "unparsed"}
        toks = re.findall(r"([STY])(\d+)[sty]", m.group(1))
        n, loc, a, b = map(int, m.groups()[1:])
        if n != 1 or len(toks) != 1:
            return {"drop": "several sites in one entity"}
        if loc != 1 or a != b:
            return {"drop": "site not localised"}
        return {"refseq": None, "gene": gene, "residue": toks[0][0], "position": int(toks[0][1])}
    if profile.startswith("brain_cptac_2020"):
        parts = (props.get("NAME") or "").split()
        acc = props.get("DESCRIPTION") or ""
        if len(parts) != 4 or not _REFSEQ.match(acc):
            return {"drop": "unparsed"}
        toks = re.findall(r"([STY])(\d+)", parts[1])
        if len(toks) != 1 or parts[3] != "1_1":
            return {"drop": "several sites in one entity"}
        return {"refseq": acc, "gene": parts[0], "residue": toks[0][0], "position": int(toks[0][1])}
    if profile.startswith(("gbm_cptac_2021", "paad_cptac_2021")):
        if ":" not in sid:
            return {"drop": "unparsed"}
        left, acc = sid.rsplit(":", 1)
        if not _REFSEQ.match(acc):
            return {"drop": "not a RefSeq protein"}
        toks = re.findall(r"(?:^|_)([STY])(\d+)(\.\d+)?(?=_|$)", left)
        if len(toks) != 1:
            return {"drop": "several sites in one entity"} if len(toks) > 1 else {"drop": "unparsed"}
        if toks[0][2]:
            return {"drop": "second entity of one site (.N suffix)"}
        return {"refseq": acc, "gene": gene, "residue": toks[0][0], "position": int(toks[0][1])}
    if profile.startswith("ucec_cptac_2020"):
        toks = sid.split("_")
        trailing = [t for t in toks[1:] if _SITE.match(t)]
        if len(trailing) != 1 or not _SITE.match(toks[-1]):
            return {"drop": "several sites in one entity"} if len(trailing) > 1 else {"drop": "unparsed"}
        m = _SITE.match(toks[-1])
        return {"refseq": None, "gene": gene, "residue": m.group(1), "position": int(m.group(2))}
    site = props.get("PHOSPHOSITE")
    if site is not None:
        m = re.match(r"^p?([STY])(\d+)$", site)
        if not m:
            return {"drop": "several sites in one entity"} if "_" in site else {"drop": "unparsed"}
        return {"refseq": None, "gene": gene, "residue": m.group(1), "position": int(m.group(2))}
    return {"drop": "unparsed"}


def transfer_position(source: str, target: str, position: int, flank: int = 7) -> tuple[int | None, str]:
    """Carry a 1-based position on ``source`` over to ``target``.

    Identical sequences carry the position unchanged. Otherwise the window of ``flank`` residues either
    side of the site (clipped at the ends) must occur exactly once in ``target``; the site then sits at the
    same offset in that occurrence. Returns ``(position, how)`` or ``(None, reason)``. The caller checks
    the residue letter.
    """
    if not 0 < position <= len(source):
        return None, "position beyond the RefSeq sequence"
    if source == target:
        return position, "identical"
    lo = max(0, position - 1 - flank)
    hi = min(len(source), position + flank)
    window = source[lo:hi]
    first = target.find(window)
    if first < 0:
        return None, "site window not found in the UniProt sequence"
    if target.find(window, first + 1) >= 0:
        return None, "site window occurs more than once in the UniProt sequence"
    return first + (position - 1 - lo) + 1, "window"


def load_tumour_abundance(path: Path = CPTAC_RESULT) -> dict[str, Any] | None:
    return load_index(path)


def relative_abundance_in_tumours(
    accession: str, position: int, result: dict[str, Any] | None
) -> dict[str, Any] | None:
    """The per-study summary of one curated phosphosite in the CPTAC tumour profiles, or None."""
    if not result:
        return None
    for row in result.get("sites", {}).get(accession, []):
        if row[0] == position:
            studies = result.get("studies") or []
            return {
                "residue": row[1],
                CPTAC_FIELD: [
                    {"study": studies[i], "tumours_with_value": n, "median_log2_ratio": med}
                    for i, n, med in row[2]
                ],
                "meaning": CPTAC_MEANING,
            }
    return None
