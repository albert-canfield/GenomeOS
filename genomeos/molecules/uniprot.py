"""UniProt (public REST, no key): the curated protein record for a gene."""

from __future__ import annotations

import json
import urllib.parse
import urllib.request
from dataclasses import dataclass, field

BASE = "https://rest.uniprot.org/uniprotkb"
FIELDS = ",".join(
    [
        "accession",
        "protein_name",
        "length",
        "sequence",
        "ft_domain",
        "ft_region",
        "ft_topo_dom",
        "ft_transmem",
        "ft_act_site",
        "ft_binding",
        "cc_function",
        "cc_subcellular_location",
    ]
)
EVIDENCE = "UniProtKB/Swiss-Prot (reviewed) via rest.uniprot.org"


@dataclass(slots=True)
class UniProtEntry:
    accession: str
    name: str
    length: int
    sequence: str
    function: str = ""
    location: str = ""
    features: list[dict] = field(default_factory=list)  # {type, description, start, end}

    def to_dict(self) -> dict:
        return {
            "accession": self.accession,
            "name": self.name,
            "length": self.length,
            "sequence": self.sequence,
            "function": self.function,
            "location": self.location,
            "features": self.features,
            "evidence": EVIDENCE,
        }


def uniprot_entry(gene: str, organism: int = 9606, timeout: int = 60) -> UniProtEntry | None:
    q = f"gene_exact:{gene} AND organism_id:{organism} AND reviewed:true"
    url = f"{BASE}/search?query={urllib.parse.quote(q)}&format=json&fields={FIELDS}&size=1"
    req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "GenomeOS/0.1"})
    with urllib.request.urlopen(req, timeout=timeout) as r:  # noqa: S310
        data = json.load(r)
    if not data.get("results"):
        return None
    d = data["results"][0]
    desc = d.get("proteinDescription", {})
    name = (
        (desc.get("recommendedName") or desc.get("submissionNames", [{}])[0])
        .get("fullName", {})
        .get("value", gene)
    )
    function = location = ""
    for c in d.get("comments", []):
        if c.get("commentType") == "FUNCTION" and c.get("texts"):
            function = c["texts"][0]["value"]
        if c.get("commentType") == "SUBCELLULAR LOCATION":
            locs = [x.get("location", {}).get("value", "") for x in c.get("subcellularLocations", [])]
            location = "; ".join(x for x in locs if x)
    feats = [
        {
            "type": f["type"],
            "description": f.get("description", ""),
            "start": f["location"]["start"]["value"],
            "end": f["location"]["end"]["value"],
        }
        for f in d.get("features", [])
        if "value" in f["location"]["start"] and "value" in f["location"]["end"]
    ]
    return UniProtEntry(
        d["primaryAccession"],
        name,
        d["sequence"]["length"],
        d["sequence"]["value"],
        function,
        location,
        feats,
    )
