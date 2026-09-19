"""The federated protein compiler.

One call compiles everything the public databases know about a protein into
one stable object, a ProteinDefinition, keyed by its UniProt accession. The
sources are adapters, each returning a normalised section with its own
evidence and confidence, and each allowed to fail without breaking the rest:

    Ensembl   gene → transcripts → protein products (which isoform each makes)
    UniProt   accession, name, sequence, function, location, isoforms, domains
              (InterPro), modifications, structures (PDB, with method), pathways
              (Reactome), disease associations, existence level
    AlphaFold predicted structure (kept separate from experimental ones)
    STRING    protein associations with the channel scores that support them
    HPA       where it is expressed (tissue, cell type) and its main location

A definition is what the protein *is*; a ProteinState is one protein in one
place at one time (tissue, cell type, level, localisation, modifications).
The two are kept apart on purpose. Compiled definitions are cached under
data/knowledge/proteins/ so simulations never call the databases; the
compiler is the importer, not the runtime.
"""

from __future__ import annotations

import json
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

UA = {"Accept": "application/json", "User-Agent": "GenomeOS/0.1 (protein compiler)"}
CACHE = Path("data/knowledge/proteins")

SOURCES = {
    "ensembl": ("Ensembl REST (rest.ensembl.org)", "curated", 0.95),
    "gencode": ("GENCODE gene models (local)", "curated", 0.95),
    "uniprot": ("UniProtKB/Swiss-Prot reviewed (rest.uniprot.org)", "curated", 0.95),
    "interpro": ("InterPro via UniProt cross-references", "curated", 0.9),
    "pdb": ("PDB via UniProt cross-references (method and resolution as deposited)", "experimental", 0.95),
    "alphafold": ("AlphaFold DB", "predicted", 0.7),
    "reactome": ("Reactome via UniProt cross-references", "curated", 0.85),
    "string": ("STRING v12 network API, combined score ≥ 0.7", "predicted", 0.5),
    "hpa": ("Human Protein Atlas entry JSON", "experimental", 0.8),
}


SECTION_SOURCE = {
    "identity": "uniprot",
    "genomic_origin": "ensembl",
    "interactions": "string",
    "expression": "hpa",
}


def _get(url: str, timeout: int = 60, retries: int = 3) -> Any:
    """GET JSON with a short back-off on server errors and rate limits (Ensembl, STRING)."""
    delay = 1.5
    for attempt in range(retries + 1):
        req = urllib.request.Request(url, headers=UA)
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:  # noqa: S310
                return json.load(r)
        except urllib.error.HTTPError as e:
            if attempt == retries or e.code not in (429, 500, 502, 503, 504):
                raise
            time.sleep(float(e.headers.get("Retry-After") or delay))
        except (urllib.error.URLError, TimeoutError):
            if attempt == retries:
                raise
            time.sleep(delay)
        delay *= 2
    raise RuntimeError("unreachable")


def _section(source: str, items: Any, **extra: Any) -> dict[str, Any]:
    desc, kind, conf = SOURCES[source]
    return {"items": items, "source": desc, "evidence": kind, "confidence": conf, **extra}


def _unavailable(source: str, err: Exception) -> dict[str, Any]:
    desc, _, _ = SOURCES[source]
    return {"items": None, "source": desc, "evidence": "none", "confidence": 0.0, "error": str(err)[:200]}


# ---- adapters ----------------------------------------------------------------------


def ensembl_gene(symbol: str) -> dict[str, Any]:
    """Gene → transcripts → protein products, from Ensembl's gene model."""
    d = _get(f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{urllib.parse.quote(symbol)}?expand=1")
    txs = []
    for t in d.get("Transcript", []):
        tr = t.get("Translation") or {}
        txs.append(
            {
                "transcript": t["id"],
                "name": t.get("display_name"),
                "biotype": t.get("biotype"),
                "canonical": bool(t.get("is_canonical")),
                "protein": tr.get("id"),
                "protein_length": tr.get("length"),
            }
        )
    txs.sort(key=lambda x: (not x["canonical"], x["name"] or ""))
    return {
        "gene_id": d["id"],
        "biotype": d.get("biotype"),
        "locus": f"chr{d.get('seq_region_name')}:{d.get('start')}-{d.get('end')}"
        f"({'+' if d.get('strand', 1) > 0 else '-'})",
        "description": d.get("description"),
        "transcripts": txs,
        "protein_products": sum(1 for x in txs if x["protein"]),
    }


UNIPROT_FIELDS = ",".join(
    [
        "accession",
        "protein_name",
        "gene_primary",
        "length",
        "sequence",
        "protein_existence",
        "ft_domain",
        "ft_region",
        "ft_topo_dom",
        "ft_transmem",
        "ft_act_site",
        "ft_binding",
        "ft_site",
        "ft_mod_res",
        "ft_lipid",
        "ft_carbohyd",
        "ft_disulfid",
        "ft_signal",
        "ft_propep",
        "ft_chain",
        "cc_function",
        "cc_subcellular_location",
        "cc_alternative_products",
        "cc_disease",
        "cc_catalytic_activity",
        "xref_pdb",
        "xref_interpro",
        "xref_reactome",
        "xref_alphafolddb",
        "xref_ensembl",
        "keyword",
    ]
)


def uniprot_raw(symbol: str, organism: int = 9606) -> dict[str, Any] | None:
    """The reviewed entry whose primary gene name is the symbol; gene_exact also matches synonyms,
    so the first hit can be another protein (MIF returned a 560-residue entry once)."""
    q = f"gene_exact:{symbol} AND organism_id:{organism} AND reviewed:true"
    url = f"https://rest.uniprot.org/uniprotkb/search?query={urllib.parse.quote(q)}&format=json&size=25"
    d = _get(url + f"&fields={UNIPROT_FIELDS}")
    hits = d.get("results") or []
    same = [
        h
        for h in hits
        if (((h.get("genes") or [{}])[0].get("geneName") or {}).get("value", "")).upper() == symbol.upper()
    ]
    if same:
        # a gene can own several reviewed entries (MIEF1: the 463-residue protein and a 70-residue
        # microprotein from an upstream ORF); the main product is the longest
        return max(same, key=lambda h: h.get("sequence", {}).get("length", 0))
    # no entry names this symbol as its primary gene. UniProt's parser treats some symbols as stop
    # words (WAS returned every reviewed human entry); a plain gene: query still resolves them.
    q2 = f"gene:{symbol} AND organism_id:{organism} AND reviewed:true"
    url2 = f"https://rest.uniprot.org/uniprotkb/search?query={urllib.parse.quote(q2)}&format=json&size=25"
    hits2 = _get(url2 + f"&fields={UNIPROT_FIELDS}").get("results") or []
    same2 = [
        h
        for h in hits2
        if (((h.get("genes") or [{}])[0].get("geneName") or {}).get("value", "")).upper() == symbol.upper()
    ]
    if same2:
        return max(same2, key=lambda h: h.get("sequence", {}).get("length", 0))
    # still nothing: the first hit is a synonym match (CAPS → CADPS); keep it, marked, so the
    # mismatch stays visible in the identity section
    if hits:
        hits[0]["_symbol_match"] = False
        return hits[0]
    return None


def normalise_uniprot(d: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """Split one UniProt entry into the compiler's sections."""
    desc = d.get("proteinDescription", {})
    name = (
        (desc.get("recommendedName") or desc.get("submissionNames", [{}])[0]).get("fullName", {}).get("value")
    )
    comments = d.get("comments", [])

    def texts(kind: str) -> list[str]:
        return [t["value"] for c in comments if c.get("commentType") == kind for t in c.get("texts", [])]

    location = [
        x.get("location", {}).get("value", "")
        for c in comments
        if c.get("commentType") == "SUBCELLULAR LOCATION"
        for x in c.get("subcellularLocations", [])
    ]
    isoforms = [
        {
            "id": i["isoformIds"][0],
            "name": i.get("name", {}).get("value"),
            "synonyms": [s["value"] for s in i.get("synonyms", [])],
            "status": i.get("isoformSequenceStatus"),
        }
        for c in comments
        if c.get("commentType") == "ALTERNATIVE PRODUCTS"
        for i in c.get("isoforms", [])
    ]
    diseases = [
        {
            "name": c["disease"].get("diseaseId"),
            "acronym": c["disease"].get("acronym"),
            "mim": (c["disease"].get("diseaseCrossReference") or {}).get("id"),
            "description": c["disease"].get("description", "")[:300],
        }
        for c in comments
        if c.get("commentType") == "DISEASE" and c.get("disease")
    ]
    feats = [
        {
            "type": f["type"],
            "description": f.get("description", ""),
            "start": f["location"]["start"].get("value"),
            "end": f["location"]["end"].get("value"),
        }
        for f in d.get("features", [])
    ]
    mod_types = {"Modified residue", "Lipidation", "Glycosylation", "Disulfide bond", "Cross-link"}
    domain_types = {
        "Domain",
        "Region",
        "Topological domain",
        "Transmembrane",
        "Active site",
        "Binding site",
        "Site",
    }
    processing_types = {"Signal", "Propeptide", "Chain", "Peptide"}
    xrefs = d.get("uniProtKBCrossReferences", [])

    def props(x: dict) -> dict[str, str]:
        return {p["key"]: p["value"] for p in x.get("properties", [])}

    pdb = []
    for x in xrefs:
        if x["database"] != "PDB":
            continue
        p = props(x)
        res = p.get("Resolution", "-")
        pdb.append(
            {
                "source": "PDB",
                "id": x["id"],
                "method": p.get("Method", "?").replace("X-ray", "X_RAY").replace("EM", "CRYO_EM").upper(),
                "resolution_A": float(res.split()[0]) if res not in ("-", "") else None,
                "chains": p.get("Chains", ""),
            }
        )
    interpro = [
        {"id": x["id"], "name": props(x).get("EntryName")} for x in xrefs if x["database"] == "InterPro"
    ]
    reactome = [
        {"id": x["id"], "name": props(x).get("PathwayName")} for x in xrefs if x["database"] == "Reactome"
    ]
    ensembl = [
        {"transcript": x["id"], "protein": props(x).get("ProteinId"), "isoform": x.get("isoformId")}
        for x in xrefs
        if x["database"] == "Ensembl"
    ]
    af = [x["id"] for x in xrefs if x["database"] == "AlphaFoldDB"]
    identity = {
        "symbol_match": d.get("_symbol_match", True),
        "accession": d["primaryAccession"],
        "name": name,
        "gene": (d.get("genes") or [{}])[0].get("geneName", {}).get("value"),
        "length": d["sequence"]["length"],
        "sequence": d["sequence"]["value"],
        "existence": d.get("proteinExistence"),
        "keywords": [k["name"] for k in d.get("keywords", [])][:30],
    }
    return {
        "identity": _section("uniprot", identity),
        "function": _section(
            "uniprot",
            {
                "summary": texts("FUNCTION"),
                "catalytic": texts("CATALYTIC ACTIVITY")[:10],
                "location": [x for x in location if x],
            },
        ),
        "isoforms": _section("uniprot", isoforms, ensembl_products=ensembl),
        "domains": _section(
            "interpro", {"interpro": interpro, "features": [f for f in feats if f["type"] in domain_types]}
        ),
        "modifications": _section("uniprot", [f for f in feats if f["type"] in mod_types]),
        "processing": _section("uniprot", [f for f in feats if f["type"] in processing_types]),
        "structures_experimental": _section("pdb", pdb, count=len(pdb)),
        "pathways": _section("reactome", reactome),
        "diseases": _section("uniprot", diseases),
        "alphafold_ids": af,
    }


def string_network(symbol: str, min_score: float = 0.7, limit: int = 50) -> list[dict[str, Any]]:
    url = (
        f"https://string-db.org/api/json/network?identifiers={urllib.parse.quote(symbol)}"
        f"&species=9606&required_score={int(min_score * 1000)}&limit={limit}"
    )
    rows = _get(url)
    out = []
    for r in rows:
        a, b = r["preferredName_A"], r["preferredName_B"]
        partner = b if a.upper() == symbol.upper() else a
        if partner.upper() == symbol.upper():
            continue
        out.append(
            {
                "partner": partner,
                "score": r["score"],
                "experimental": r.get("escore", 0),
                "database": r.get("dscore", 0),
                "textmining": r.get("tscore", 0),
                "coexpression": r.get("ascore", 0),
                "physical_evidence": r.get("escore", 0) >= 0.4,
            }
        )
    best: dict[str, dict[str, Any]] = {}
    for x in out:
        if x["partner"] not in best or x["score"] > best[x["partner"]]["score"]:
            best[x["partner"]] = x
    return sorted(best.values(), key=lambda x: -x["score"])


def _numeric(d: Any) -> Any:
    """HPA gives per-tissue values as strings; keep them as numbers."""
    if not isinstance(d, dict):
        return d
    out = {}
    for k, v in d.items():
        try:
            out[k] = float(v)
        except (TypeError, ValueError):
            out[k] = v
    return out


def hpa_entry(ensembl_gene_id: str) -> dict[str, Any]:
    d = _get(f"https://www.proteinatlas.org/{ensembl_gene_id}.json")
    return {
        "protein_class": d.get("Protein class"),
        "protein_evidence": d.get("Evidence"),
        "subcellular_main": d.get("Subcellular main location"),
        "subcellular_additional": d.get("Subcellular additional location"),
        "tissue_specificity": d.get("RNA tissue specificity"),
        "tissue_distribution": d.get("RNA tissue distribution"),
        "tissue_ntpm": _numeric(d.get("RNA tissue specific nTPM")),
        "cell_type_specificity": d.get("RNA single cell type specificity"),
        "cell_type_ncpm": _numeric(d.get("RNA single cell type specific nCPM")),
        "blood_cell_specificity": d.get("RNA blood cell specificity"),
        "cancer_specificity": d.get("RNA cancer specificity"),
        "disease_involvement": d.get("Disease involvement"),
    }


# ---- the compiled object -------------------------------------------------------------


@dataclass(slots=True)
class ProteinState:
    """One protein in one place at one time. Separate from the definition on purpose."""

    protein: str  # "UniProt:P04637"
    tissue: str | None = None
    cell_type: str | None = None
    level: str | float | None = None  # HPA category or a measured concentration
    localisation: str | None = None
    modifications: list[str] = field(default_factory=list)
    isoform: str | None = None
    evidence: str = "none"
    confidence: float = 0.0

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def compile_protein(
    symbol: str,
    cache_dir: Path = CACHE,
    refresh: bool = False,
    sources: set[str] | None = None,
    origin: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Compile one protein definition from the federated sources; cached by gene symbol.
    `origin` (gene id, locus, transcripts from local gene models) replaces the Ensembl call."""
    symbol = symbol.upper()
    want = sources or set(SOURCES)
    if origin is not None:
        want = want - {"ensembl"}
    cache_dir.mkdir(parents=True, exist_ok=True)
    cached = cache_dir / f"{symbol}.json"
    prior: dict[str, Any] | None = None
    if cached.exists() and not refresh:
        prior = json.loads(cached.read_text())
        failed = {k for k, v in prior["sections"].items() if v.get("error") and v.get("evidence") == "none"}
        # a section never attempted (its prerequisite failed) counts as failed too
        failed |= {
            k
            for k in SECTION_SOURCE
            if k not in prior["sections"] and SECTION_SOURCE[k] in (sources or SOURCES)
        }
        if not failed:
            return prior
        # a transient failure must not become permanent: refetch only what failed
        want = {SECTION_SOURCE.get(k, k) for k in failed}
    out: dict[str, Any] = prior or {"id": None, "gene": symbol, "sections": {}}
    sec = out["sections"]
    raw = None
    if prior is None or "uniprot" in want:
        try:
            raw = uniprot_raw(symbol)
        except (urllib.error.URLError, TimeoutError, OSError, ValueError) as e:
            sec["identity"] = _unavailable("uniprot", e)
    if raw:
        parts = normalise_uniprot(raw)
        af_ids = parts.pop("alphafold_ids")
        sec.update(parts)
        out["id"] = f"UniProt:{parts['identity']['items']['accession']}"
        if "alphafold" in want:
            sec["structures_predicted"] = _section(
                "alphafold",
                [
                    {"source": "AlphaFold", "id": a, "method": "PREDICTED", "confidence": "pLDDT per residue"}
                    for a in af_ids
                ],
            )
    elif "identity" not in sec and prior is None:
        sec["identity"] = {
            "items": None,
            "source": SOURCES["uniprot"][0],
            "evidence": "none",
            "confidence": 0.0,
            "error": "no reviewed entry",
        }
    ens = None
    if origin is not None:
        sec["genomic_origin"] = _section("gencode", origin)
        ens = origin
    elif prior is not None and "ensembl" not in want:
        ens = (sec.get("genomic_origin") or {}).get("items")
    if "ensembl" in want:
        try:
            ens = ensembl_gene(symbol)
            sec["genomic_origin"] = _section("ensembl", ens)
        except (urllib.error.URLError, TimeoutError, OSError, ValueError, KeyError) as e:
            sec["genomic_origin"] = _unavailable("ensembl", e)
    if "string" in want:
        try:
            sec["interactions"] = _section("string", string_network(symbol))
        except (urllib.error.URLError, TimeoutError, OSError, ValueError, KeyError) as e:
            sec["interactions"] = _unavailable("string", e)
    if "hpa" in want and ens:
        try:
            sec["expression"] = _section("hpa", hpa_entry(ens["gene_id"]))
        except (urllib.error.URLError, TimeoutError, OSError, ValueError) as e:
            sec["expression"] = _unavailable("hpa", e)
    out["coverage"] = coverage(out)
    cached.write_text(json.dumps(out))
    return out


def coverage(defn: dict[str, Any]) -> dict[str, bool]:
    """Which questions this compiled definition answers, in Albert's table order."""
    s = defn["sections"]

    def has(name: str, key: str | None = None) -> bool:
        items = (s.get(name) or {}).get("items")
        if items is None:
            return False
        if key:
            v = items.get(key) if isinstance(items, dict) else None
            return bool(v)
        return bool(items)

    return {
        "genomic_origin": has("genomic_origin"),
        "sequence": has("identity", "sequence"),
        "name": has("identity", "name"),
        "function": has("function", "summary"),
        "domains": has("domains", "interpro") or has("domains", "features"),
        "pathways": has("pathways"),
        "interactions": has("interactions"),
        "expression": has("expression", "tissue_specificity") or has("expression", "subcellular_main"),
        "structure_experimental": has("structures_experimental"),
        "structure_predicted": has("structures_predicted"),
        "disease": has("diseases"),
    }


def states_from_definition(defn: dict[str, Any]) -> list[ProteinState]:
    """Derive ProteinState records (where/at what level) from the HPA section of a definition."""
    exp = (defn["sections"].get("expression") or {}).get("items") or {}
    pid = defn.get("id") or f"gene:{defn['gene']}"
    loc = ", ".join(exp.get("subcellular_main") or []) or None
    out = []
    for tissue, ntpm in (exp.get("tissue_ntpm") or {}).items():
        out.append(
            ProteinState(
                pid,
                tissue=tissue,
                level=ntpm,
                localisation=loc,
                evidence="experimental: HPA RNA nTPM",
                confidence=0.7,
            )
        )
    for ct, ncpm in (exp.get("cell_type_ncpm") or {}).items():
        out.append(
            ProteinState(
                pid,
                cell_type=ct,
                level=ncpm,
                localisation=loc,
                evidence="experimental: HPA single-cell nCPM",
                confidence=0.7,
            )
        )
    if not out and loc:
        out.append(
            ProteinState(pid, localisation=loc, evidence="experimental: HPA subcellular", confidence=0.7)
        )
    return out


def to_biolang(defn: dict[str, Any], max_items: int = 12) -> str:
    """One BioLang `protein` block from a compiled definition (the importer's output)."""
    s = defn["sections"]
    ident = (s.get("identity") or {}).get("items") or {}
    if not ident:
        return f"# {defn['gene']}: no reviewed UniProt entry\n"
    dom = (s.get("domains") or {}).get("items") or {}
    names = [d["name"] for d in dom.get("interpro", []) if d.get("name")][:max_items]
    pw = [x["id"] for x in (s.get("pathways") or {}).get("items") or []][:max_items]
    inter = [x["partner"] for x in (s.get("interactions") or {}).get("items") or [] if x["physical_evidence"]]
    iso = [i["id"] for i in (s.get("isoforms") or {}).get("items") or []]
    pdb = [x["id"] for x in (s.get("structures_experimental") or {}).get("items") or []][:max_items]
    af = [f"AF-{x['id']}" for x in (s.get("structures_predicted") or {}).get("items") or []]
    fn = ((s.get("function") or {}).get("items") or {}).get("summary") or []
    sources = "UniProt, Ensembl, InterPro, PDB, AlphaFold, Reactome, STRING, HPA"
    lines = [f"# {ident.get('name')} — compiled from {sources}"]
    if fn:
        lines.append("# " + fn[0][:160].replace("\n", " "))
    props = [
        f"accession: {ident['accession']}",
        f"sequence: {ident['sequence']}",
    ]
    if iso:
        props.append(f"isoforms: {', '.join(iso)}")
    if names:
        props.append(f"domains: {', '.join(names)}")
    if pw:
        props.append(f"pathways: {', '.join(pw)}")
    if inter:
        props.append(f"interactions: {', '.join(inter[:max_items])}")
    if pdb or af:
        props.append(f"structures: {', '.join(pdb + af)}")
    props.append('evidence: curated "UniProtKB/Swiss-Prot; InterPro; Reactome; PDB; STRING physical channel"')
    props.append(f"confidence: {min(v['confidence'] for v in s.values() if v.get('items'))}")
    lines.append(f"protein {defn['gene']} {{")
    lines.extend(f"  {x};" for x in props)
    lines.append("}")
    lines.extend(writer_rules(defn["gene"], writer_counts(defn)))
    return "\n".join(lines) + "\n"


def writer_counts(defn: dict[str, Any]) -> dict[str, int]:
    """Writers UniProt names on the protein's modified residues (kinases, acetyltransferases, …) with
    how many sites each one writes; the post-translational layer's edges, per protein."""
    from genomeos.molecules.ptm import sites

    counts: dict[str, int] = {}
    for s in sites(defn):
        for w in s["writers"]:
            if w != "(self)" and w != defn["gene"]:
                counts[w] = counts.get(w, 0) + 1
    return counts


def writer_rules(gene: str, counts: dict[str, int], max_rules: int = 12) -> list[str]:
    """The writers as BioLang rules: `rule PKA modifies TP53 { … }`, curated from UniProt, so a program
    that imports the protein carries who acts on it. A rule may only name declared entities, so each
    writer is declared first as a bare protein (a program that also imports the writer itself should
    drop the stub, since an id is declared once)."""
    out = []
    ranked = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))[:max_rules]
    for w, _ in ranked:
        ev = f'evidence: curated "UniProt: named as a writer on {gene}"'
        out.append(f"protein {w} {{ {ev}; confidence: 0.6 }}")
    for w, n in ranked:
        sites_ = f"{n} modified residue{'s' if n != 1 else ''}"
        ev = f'evidence: curated "UniProt: {sites_} written by {w}"'
        out.append(f"rule {w} modifies {gene} {{ strength: 1.0; {ev}; confidence: 0.8 }}")
    return out
