"""The human body as counted populations: the numbers a population program needs.

Sender & Milo 2021 (Nat Med 27:45, "The distribution of cellular turnover in
the human body"; data github.com/milo-lab/cellular_turnover, Summary.xlsx):
per cell type the number of cells in a reference adult, the lifespan in
days, the turnover (cells per day) and the cell mass. Streamed once and
distilled to one small table; the body program under data/organisms/human
cites these numbers as curated counts, not mechanism.
"""

from __future__ import annotations

import io
import math
import urllib.request
import xml.etree.ElementTree as ET
import zipfile

MILO_URL = "https://raw.githubusercontent.com/milo-lab/cellular_turnover/master/Summary.xlsx"
_NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
_M = "{" + _NS["m"] + "}"
_R = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"


def read_xlsx(data: bytes) -> dict[str, list[list[str]]]:
    """Every sheet of an .xlsx as rows of strings, with the standard library only."""
    z = zipfile.ZipFile(io.BytesIO(data))
    names = z.namelist()
    shared: list[str] = []
    if "xl/sharedStrings.xml" in names:
        for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall("m:si", _NS):
            shared.append("".join(t.text or "" for t in si.iter(_M + "t")))
    rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
    target = {r.get("Id"): r.get("Target") for r in rels}
    wb = ET.fromstring(z.read("xl/workbook.xml"))

    def value(c: ET.Element) -> str:
        t = c.get("t")
        if t == "inlineStr":
            return "".join(x.text or "" for x in c.iter(_M + "t"))
        v = c.find("m:v", _NS)
        if v is None or v.text is None:
            return ""
        return shared[int(v.text)] if t == "s" else v.text

    out: dict[str, list[list[str]]] = {}
    for sh in wb.iter(_M + "sheet"):
        path = target[sh.get(_R + "id")].lstrip("/")
        path = path if path.startswith("xl/") else "xl/" + path
        root = ET.fromstring(z.read(path))
        rows = []
        for r in root.findall(".//m:sheetData/m:row", _NS):
            rows.append([value(c) for c in r.findall("m:c", _NS)])
        out[sh.get("name")] = rows
    return out


def fetch_milo(url: str = MILO_URL, timeout: float = 60.0) -> bytes:
    with urllib.request.urlopen(url, timeout=timeout) as r:  # noqa: S310 (fixed public URL)
        return r.read()


def _num(s: str) -> float | None:
    try:
        return float(s)
    except ValueError:
        return None


def distil(sheets: dict[str, list[list[str]]]) -> dict:
    colours = sheets.get("Colors", [])
    tissue_of = {r[0]: (r[1] if len(r) > 1 else "") for r in colours[1:] if r}
    types = []
    for name, rows in sheets.items():
        if name == "Colors" or not rows or not rows[0] or rows[0][0] != "Parameter":
            continue
        params: dict[str, dict] = {}
        for r in rows[1:]:
            if len(r) < 3 or not r[0]:
                continue
            params[r[0]] = {
                "value": _num(r[1]) if len(r) > 1 else None,
                "unit": r[2],
                "uncertainty": _num(r[3]) if len(r) > 3 else None,
            }
        number = params.get("number", {}).get("value")
        lifespan = params.get("lifespan", {}).get("value")
        turnover = params.get("cellular turnover rate", {}).get("value")
        if turnover is None:
            turnover = params.get("extrapolated cellular turnover rate", {}).get("value")
            extrapolated = turnover is not None
        else:
            extrapolated = False
        types.append(
            {
                "cell_type": name,
                "tissue": tissue_of.get(name, ""),
                "cells": number,
                "lifespan_days": lifespan,
                "turnover_per_day": turnover,
                "turnover_extrapolated_from_rodents": extrapolated,
                "cell_mass_pg": params.get("cell mass", {}).get("value"),
                "total_mass_g": params.get("total cellular mass", {}).get("value"),
            }
        )
    total_cells = sum(t["cells"] or 0 for t in types)
    total_turnover = sum(t["turnover_per_day"] or 0 for t in types)
    return {
        "source": "Sender & Milo 2021, Nat Med 27:45; github.com/milo-lab/cellular_turnover Summary.xlsx",
        "url": MILO_URL,
        "cell_types": len(types),
        "total_cells": total_cells,
        "total_turnover_per_day": total_turnover,
        "types": sorted(types, key=lambda t: -(t["cells"] or 0)),
    }


# Sender & Milo cell types grouped into the tissue populations of the body program, with the germ
# layer each derives from (curated developmental biology) and a Cell Ontology id for the group.
GROUPS: list[dict] = [
    {"id": "Erythrocyte", "layer": "Mesoderm", "cl": "CL:0000232", "types": ["Erythrocytes"]},
    {"id": "Neutrophil", "layer": "Mesoderm", "cl": "CL:0000775", "types": ["Neutrophils"]},
    {
        "id": "Lymphocyte",
        "layer": "Mesoderm",
        "cl": "CL:0000542",
        "types": [
            "Mature T cells",
            "Mature B cells",
            "B cells progenitors",
            "Transitional B cells",
            "Thymocytes",
        ],
    },
    {
        "id": "Monocyte",
        "layer": "Mesoderm",
        "cl": "CL:0000576",
        "types": ["Monocytes", "Alveolar macrophages", "Kupffer cells"],
    },
    {"id": "EndothelialCell", "layer": "Mesoderm", "cl": "CL:0000115", "types": ["Endothelial cells"]},
    {"id": "FibroblastDermal", "layer": "Mesoderm", "cl": "CL:0000057", "types": ["Dermal fibroblasts"]},
    {"id": "Adipocyte", "layer": "Mesoderm", "cl": "CL:0000136", "types": ["Adipocytes"]},
    {"id": "CardiomyocyteHuman", "layer": "Mesoderm", "cl": "CL:0000746", "types": ["Cardiomyocytes"]},
    {"id": "Myocyte", "layer": "Mesoderm", "cl": "CL:0000187", "types": ["Myocytes"]},
    {"id": "LungInterstitial", "layer": "Mesoderm", "cl": "CL:0000499", "types": ["Lung interstitial cells"]},
    {"id": "NeuronHuman", "layer": "Ectoderm", "cl": "CL:0000540", "types": ["Neurons"]},
    {"id": "Glia", "layer": "Ectoderm", "cl": "CL:0000125", "types": ["Glial cells"]},
    {"id": "Keratinocyte", "layer": "Ectoderm", "cl": "CL:0000312", "types": ["Epidermal cells"]},
    {"id": "Hepatocyte", "layer": "Endoderm", "cl": "CL:0000182", "types": ["Hepatocytes", "Stellate cells"]},
    {
        "id": "GutEpithelium",
        "layer": "Endoderm",
        "cl": "CL:0002563",
        "types": ["Small intestine epithelia", "Colon epithelia", "Stomach epithelia"],
    },
    {
        "id": "LungEpithelium",
        "layer": "Endoderm",
        "cl": "CL:0000082",
        "types": ["Alveolar epithelial cells", "Bronchial epithelial cells"],
    },
]


# Birth: about 1.8e12 cells (a 3.5 kg newborn at the adult cell density of Hatton 2023: 36e12 cells in
# 70 kg), which is also what the fetal doublings in body.bio produce. Each tissue's share at birth follows
# from the germ-layer shares in body.bio and the adult shares within each layer; its growth rate is then
# whatever takes it from that share to the adult count over childhood (inferred; the adult count caps the
# growth, so tissues over-represented at birth simply stop, and a rate aimed at 16.5 years leaves room).
BIRTH_TOTAL = 1.8e12
LAYER_SHARE = {"Ectoderm": 0.45, "Mesoderm": 0.35, "Endoderm": 0.20}
CHILDHOOD_DAYS = 6000.0


def to_bio_tissues(table: dict) -> str:
    """The adult tissue populations as BioLang: cell types, germ-layer shares at organogenesis,
    growth to the adult count, and daily turnover (loss and replacement) once born."""
    by_type = {t["cell_type"]: t for t in table["types"]}
    src = "Sender & Milo 2021, Nat Med 27:45"
    lines = [
        "# Generated by `genomeos data distil --only human_cell_turnover` from Sender & Milo 2021.",
        "# Adult counts, germ-layer shares (share of the layer's adult cell count), growth to the",
        "# adult count and daily turnover per tissue population. Counts, not mechanism.",
        "module organism.human.tissues",
        "import bio.std.development",
        "",
    ]
    for layer, cl in (("Ectoderm", "CL:0000221"), ("Mesoderm", "CL:0000222"), ("Endoderm", "CL:0000223")):
        lines.append(
            f"cell_type {layer} {{ parent: Blastomere; ontology: {cl}; "
            f'evidence: curated "Cell Ontology"; confidence: 0.9 }}'
        )
    lines.append("")
    groups = []
    for g in GROUPS:
        cells = sum(by_type[t]["cells"] or 0 for t in g["types"])
        turnover = sum(by_type[t]["turnover_per_day"] or 0 for t in g["types"])
        extrapolated = any(by_type[t]["turnover_extrapolated_from_rodents"] for t in g["types"])
        groups.append({**g, "cells": cells, "turnover": turnover, "extrapolated": extrapolated})
        lines.append(
            f"cell_type {g['id']} {{ parent: {g['layer']}; ontology: {g['cl']}; "
            f'evidence: curated "Cell Ontology; grouped from {", ".join(g["types"])}"; confidence: 0.8 }}'
        )
    lines.append("")
    layer_total = {layer: sum(g["cells"] for g in groups if g["layer"] == layer) for layer in LAYER_SHARE}
    for layer in LAYER_SHARE:
        members = sorted((g for g in groups if g["layer"] == layer), key=lambda g: -g["cells"])
        remaining = layer_total[layer]
        lines.append(f"# {layer}: {remaining:.3g} adult cells across {len(members)} populations")
        for i, g in enumerate(members):
            frac = 1.0 if i == len(members) - 1 else g["cells"] / remaining
            remaining -= g["cells"]
            lines.append(
                f"decision {layer.lower()}_to_{g['id']} {{ action: differentiate; when: cell_type = {layer}, "
                f"stage = Organogenesis; to: {g['id']}; fraction: {frac:.4f}; "
                f'evidence: inferred "share of the layer\'s adult cell count ({src})"; confidence: 0.4 }}'
            )
        lines.append("")
    lines.append("# growth stops at the adult count (quiescence re-evaluated at every growth step)")
    for g in groups:
        lines.append(
            f"decision {g['id']}_adult_count {{ action: quiesce; when: cell_type = {g['id']}, "
            f"count = >={g['cells']:.3g}, stage = Fetal|Childhood; "
            f'evidence: curated "{src}: {g["cells"]:.3g} cells in the reference adult"; confidence: 0.7 }}'
        )
    lines.append("")
    lines.append("# turnover: the daily share lost, the share produced to replace it, and in childhood the")
    lines.append(f"# extra growth that takes each tissue from its share of ~{BIRTH_TOTAL:.2g} cells at birth")
    lines.append("# to its adult count by 18 years (inferred rates; the adult count caps the growth)")
    for g in groups:
        cid = g["id"]
        f = g["turnover"] / g["cells"] if g["cells"] else 0.0
        r = f / (1.0 - f)
        kind = "inferred" if g["extrapolated"] else "experimental"
        note = "lifespan extrapolated from rodents" if g["extrapolated"] else "measured lifespan"
        conf = 0.4 if g["extrapolated"] else 0.7
        at_birth = BIRTH_TOTAL * LAYER_SHARE[g["layer"]] * g["cells"] / layer_total[g["layer"]]
        growth = max(0.0, math.log(g["cells"] / at_birth) / CHILDHOOD_DAYS) if at_birth > 0 else 0.0
        if f > 0:
            lines.append(
                f"decision {cid}_loss {{ action: die; when: cell_type = {cid}, stage = Childhood|Adult; "
                f"fraction: {f:.4g}; after: 1 d; "
                f'evidence: {kind} "{src}: {g["turnover"]:.3g} cells/day, {note}"; confidence: {conf} }}'
            )
        lines.append(
            f"decision {cid}_childhood {{ action: divide; when: cell_type = {cid}, stage = Childhood; "
            f"fraction: {r + growth:.4g}; after: 1 d; "
            f'evidence: inferred "replacement plus {growth:.2g}/day growth: ~{at_birth:.2g} at birth to '
            f'{g["cells"]:.3g} by 18 years"; confidence: 0.3 }}'
        )
        if f > 0:
            lines.append(
                f"decision {cid}_replace {{ action: divide; when: cell_type = {cid}, stage = Adult; "
                f"fraction: {r:.4g}; after: 1 d; "
                f'evidence: {kind} "{src}: production balances loss"; confidence: {conf} }}'
            )
    return "\n".join(lines) + "\n"
