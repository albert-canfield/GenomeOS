"""A gene dossier: every layer GenomeOS knows about one gene, in one document.

`genomeos report APP --chrom chr21` (or `/api/report`) assembles, from local
data and the local knowledge cache only:

    origin        locus, strand, size, transcripts, canonical (GENCODE)
    flow          spliced mRNA, CDS, our protein, identity with UniProt
    regulation    the node it sits in, promoter elements, reachable enhancers
    RNA           GTEx expression pattern and top tissues (if fetched)
    protein       function, domains, structures, pathways, partners, diseases
    pathways      reactions lost across its pathways when absent (if run)
    individual    the test human's variants inside the gene, coding ones traced

Every section names its evidence. Markdown out; the JSON behind it is the
same dict.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from genomeos.molecules.compiler import CACHE


def _load_json(p: Path) -> dict[str, Any] | None:
    try:
        return json.loads(p.read_text()) if p.exists() else None
    except json.JSONDecodeError:
        return None


def gene_report(symbol: str, chrom: str, annotation, genome, root: Path = Path(".")) -> dict[str, Any]:
    from genomeos.flow import trace_gene
    from genomeos.genome.regulation import regulation_of
    from genomeos.genome.regulatory import load_ccres

    symbol = symbol.upper()
    g = annotation.gene(symbol)
    module = annotation.to_module("report")
    rep: dict[str, Any] = {"gene": g.symbol, "chrom": chrom, "sections": {}}
    sec = rep["sections"]
    txs = module.entities[g.id].transcripts
    sec["origin"] = {
        "locus": str(g.locus),
        "strand": g.locus.strand.value,
        "length_bp": g.locus.length,
        "type": g.type,
        "transcripts": len(txs),
        "coding_transcripts": sum(1 for t in txs if t.cds_segments),
        "evidence": "curated: GENCODE",
    }
    tr = trace_gene(genome, g, txs)
    defn = _load_json(CACHE / f"{symbol}.json")
    ident = ((defn or {}).get("sections", {}).get("identity") or {}).get("items") or {}
    if tr is not None:
        uni = ident.get("sequence") or ""
        ident_frac = None
        if uni:
            n = max(len(tr.protein), len(uni))
            ident_frac = (
                round(sum(1 for a, b in zip(tr.protein, uni, strict=False) if a == b) / n, 4) if n else None
            )
        sec["flow"] = {
            "canonical": tr.transcript,
            "exons": len(tr.exons),
            "mrna_nt": len(tr.mrna),
            "utr5": tr.utr5,
            "cds_nt": len(tr.cds),
            "utr3": tr.utr3,
            "protein_aa": len(tr.protein),
            "codon_table": tr.table_name,
            "identity_with_uniprot": ident_frac,
            "uniprot_aa": len(uni) if uni else None,
            "evidence": "derived: splice and translation of the reference; UniProt for the comparison",
        }
    ccres = load_ccres(chrom)
    if ccres:
        r = regulation_of(symbol, chrom, ccres, annotation, genome.lengths[chrom])
        d = r["domain"]
        sec["regulation"] = {
            "node": d["id"] if d else None,
            "node_length": d["length"] if d else None,
            "genes_in_node": d["coding_genes"] if d else None,
            "promoter_elements": len(r["promoters"]),
            "enhancers_in_node": r["enhancers_in_domain"],
            "enhancers_nearest": r["enhancers_nearest_to_this_gene"],
            "enhancers_inside_gene": r["enhancers_inside_gene"],
            "insulators": len(r["insulators_bounding"]),
            "enhancers_predicted": r.get("enhancers_predicted", 0),
            "enhancers_predicted_this_gene": r.get("enhancers_predicted_this_gene", 0),
            "strongest_predicted": next(
                (
                    {
                        "element": e["id"],
                        "distance": e["distance"],
                        "log2_fold_change": e["predicted"]["log2_fold_change"],
                        "tissue": e["predicted"]["tissue"],
                        "action": e["predicted"]["action"],
                    }
                    for e in sorted(
                        (e for e in r["enhancers"] if e.get("predicted_this_gene")),
                        key=lambda e: -abs(e["predicted"]["log2_fold_change"]),
                    )
                ),
                None,
            ),
            "evidence": "curated: ENCODE elements; inferred: targets bounded by the CTCF domain"
            + (
                "; predicted: AlphaGenome deletion effects where scored"
                if r.get("enhancers_predicted")
                else ""
            ),
        }
    # the reader: which cell types read this gene (promoter open), from every reader result on this chromosome
    read_in, silent_in = [], []
    for f in sorted((root / "data" / "results").glob(f"reader_*_{chrom}.json")):
        rr = _load_json(f)
        if not rr or "silent_genes" not in rr:
            continue
        cell = rr.get("cell_type") or f.name[len("reader_") : -len(f"_{chrom}.json")]
        (silent_in if symbol in set(rr["silent_genes"]) else read_in).append(cell)
    if read_in or silent_in:
        sec["reader"] = {
            "read_in": read_in,
            "silent_in": silent_in,
            "evidence": "experimental: ENCODE DNase peaks; inferred: read = promoter (TSS ± 1 kb) open",
        }
    gtex = _load_json(root / "data" / "knowledge" / "expression" / f"gtex_{symbol}.json")
    if gtex and gtex.get("tissues"):
        sec["rna"] = {
            "pattern": gtex["pattern"],
            "median_tpm": gtex["median_tpm"],
            "top": gtex["top"][:6],
            "tissues_expressed": gtex["tissues_expressed"],
            "evidence": gtex["evidence"],
        }
    if defn:
        s = defn["sections"]
        fn = (s.get("function") or {}).get("items") or {}
        dom = (s.get("domains") or {}).get("items") or {}
        exp = s.get("structures_experimental") or {}
        inter = [x for x in ((s.get("interactions") or {}).get("items") or []) if x["physical_evidence"]]
        sec["protein"] = {
            "accession": ident.get("accession"),
            "name": ident.get("name"),
            "length": ident.get("length"),
            "existence": ident.get("existence"),
            "function": (fn.get("summary") or [""])[0][:600],
            "location": (fn.get("location") or [])[:4],
            "domains": [x["name"] for x in dom.get("interpro", [])][:10],
            "structures_experimental": exp.get("count", 0),
            "alphafold": bool((s.get("structures_predicted") or {}).get("items")),
            "pathways": len((s.get("pathways") or {}).get("items") or []),
            "pathway_names": [p["name"] for p in ((s.get("pathways") or {}).get("items") or [])[:6]],
            "partners_physical": list(dict.fromkeys(x["partner"] for x in inter))[:10],
            "diseases": [x["name"] for x in ((s.get("diseases") or {}).get("items") or []) if x.get("name")][
                :6
            ],
            "evidence": "curated: UniProt, InterPro, Reactome, PDB; predicted: AlphaFold, STRING",
        }
    ko = _load_json(root / "data" / "results" / f"knockout_{symbol}.json")
    if ko:
        rows = ko["pathways"]
        sec["knockout"] = {
            "pathways": len(rows),
            "reactions_lost": sum(r["lost"] for r in rows),
            "most_affected": [(r["name"], r["fraction_lost"]) for r in rows[:3]],
            "evidence": "curated: Reactome; inferred reachability",
        }
    # every local individual inside the gene: the test human and any imported genome
    from genomeos.genome.individuals import rows_in, sources

    people = []
    for name, vcf, evidence in sources(chrom):
        inside, coding = 0, []
        for f in rows_in(vcf, g.locus.start + 1, g.locus.end):
            inside += 1
            if tr is not None and len(f[3]) == 1 and len(f[4]) == 1 and len(coding) < 20:
                sub = tr.substitute(int(f[1]) - 1, f[3], f[4])
                if sub.get("region") == "CDS":
                    coding.append(
                        {
                            "pos": int(f[1]),
                            "change": f"{f[3]}>{f[4]}",
                            "genotype": f[9].split(":")[0] if len(f) > 9 else "",
                            "consequence": sub.get("consequence"),
                            "hgvs_p": sub.get("hgvs_p"),
                        }
                    )
        people.append(
            {
                "sample": name,
                "variants_in_gene": inside,
                "coding_snvs": coding,
                "evidence": f"{evidence}; consequence derived by the local trace",
            }
        )
    if people:
        sec["individual"] = people[0]  # the test human, as before
        sec["individuals"] = people
    return rep


def to_markdown(rep: dict[str, Any]) -> str:
    sec = rep["sections"]
    lines = [f"# {rep['gene']} ({rep['chrom']})", ""]

    def para(title: str, body: str, evidence: str) -> None:
        lines.extend([f"**{title}** — {body} _{evidence}_", ""])

    o = sec["origin"]
    para(
        "Origin",
        f"{o['locus']} ({o['strand']}), {o['length_bp']:,} bp, {o['type']}; {o['transcripts']} transcripts, "
        f"{o['coding_transcripts']} coding.",
        o["evidence"],
    )
    if "flow" in sec:
        f = sec["flow"]
        idn = ""
        if f["identity_with_uniprot"] is not None:
            idn = f" Identity with UniProt {f['identity_with_uniprot']:.1%} ({f['uniprot_aa']} aa)."
        para(
            "DNA → RNA → protein",
            f"canonical {f['canonical']}, {f['exons']} exons; mRNA {f['mrna_nt']:,} nt = 5'UTR {f['utr5']} + "
            f"CDS {f['cds_nt']} + 3'UTR {f['utr3']}; protein {f['protein_aa']} aa "
            f"({f['codon_table']} code).{idn}",
            f["evidence"],
        )
    if "regulation" in sec:
        r = sec["regulation"]
        predicted = ""
        if r.get("enhancers_predicted"):
            predicted = (
                f" AlphaGenome scored {r['enhancers_predicted']} of them, "
                f"{r['enhancers_predicted_this_gene']} name this gene"
            )
            sp = r.get("strongest_predicted")
            if sp:
                predicted += (
                    f", strongest {sp['element']} at {sp['distance']:,} bp "
                    f"({sp['log2_fold_change']:+.2f} in {sp['tissue']})"
                )
            predicted += "."
        para(
            "Regulation",
            f"node {r['node']} ({(r['node_length'] or 0) / 1000:.0f} kb, {r['genes_in_node']} coding genes); "
            f"{r['promoter_elements']} promoter elements; {r['enhancers_in_node']} enhancers can reach it "
            f"({r['enhancers_nearest']} nearest to it, {r['enhancers_inside_gene']} inside the gene); "
            f"{r['insulators']} bounding insulators." + predicted,
            r["evidence"],
        )
    if "reader" in sec:
        r = sec["reader"]
        n_all = len(r["read_in"]) + len(r["silent_in"])
        para(
            "Reader",
            f"promoter open (read) in {len(r['read_in'])} of {n_all} cell types: "
            f"{', '.join(r['read_in']) or 'none'}; silent in: {', '.join(r['silent_in']) or 'none'}.",
            r["evidence"],
        )
    if "rna" in sec:
        r = sec["rna"]
        top = ", ".join(f"{t} {v:.0f}" for t, v in r["top"])
        para(
            "Expression",
            f"{r['pattern']}; median {r['median_tpm']} TPM, {r['tissues_expressed']} tissues ≥ 1 TPM; "
            f"top: {top}.",
            r["evidence"],
        )
    if "protein" in sec:
        p = sec["protein"]
        lines.extend([f"**Protein** — {p['accession']} {p['name']}, {p['length']} aa, {p['existence']}.", ""])
        lines.extend([p["function"], ""])
        af = ", AlphaFold" if p["alphafold"] else ""
        para(
            "Protein facts",
            f"Location: {'; '.join(p['location']) or '?'}. Domains: {', '.join(p['domains']) or 'none'}. "
            f"Structures: {p['structures_experimental']} experimental{af}. "
            f"Pathways: {p['pathways']} ({'; '.join(p['pathway_names'])}). "
            f"Physical partners: {', '.join(p['partners_physical']) or 'none recorded'}. "
            f"Diseases: {'; '.join(p['diseases']) or 'none recorded'}.",
            p["evidence"],
        )
    if "knockout" in sec:
        k = sec["knockout"]
        worst = "; ".join(f"{n} ({fr:.0%})" for n, fr in k["most_affected"])
        para(
            "Without it",
            f"{k['reactions_lost']} reactions lost across {k['pathways']} pathways; most affected: {worst}.",
            k["evidence"],
        )
    for i in sec.get("individuals") or ([sec["individual"]] if "individual" in sec else []):
        cod = (
            "; ".join(f"{c['hgvs_p'] or c['consequence']} ({c['genotype']})" for c in i["coding_snvs"])
            or "none"
        )
        para(
            i["sample"],
            f"{i['variants_in_gene']:,} variants inside the gene; coding SNVs: {cod}.",
            i["evidence"],
        )
    return "\n".join(lines)
