# SPDX-License-Identifier: AGPL-3.0-or-later
"""Origin per gene, age per library, paralogues: Ensembl Compara streamed once and distilled.

The conversation behind ROADMAP area J asked for the library ladder the genome carries
(`Life.Core`, `Eukaryote.Core`, `Animal.Core`, `Vertebrate.Core`, `Mammal.Core`, primate,
human-specific): for every gene, the deepest clade in which an orthologue exists is the
gene's *origin*, and a library's members' origins say how old the library is and whether
it is one library at all. Ensembl Compara's homology dump for human (one gzipped TSV,
109 MB, `Compara.<release>.protein_default.homologies.tsv.gz`) lists every orthologue and
paralogue of every human protein-coding gene with the partner species; streamed once over
HTTP and never stored, it gives per gene:

- origin: the deepest stratum of a ladder of clades shared with human (Eukaryota down to
  Homo) that has an orthologue; a gene with no orthologue in Ensembl's species set is
  `Homo` (human-specific as far as this set sees, which for Ensembl vertebrates means:
  absent from yeast, fly, worm, ciona, lamprey and every vertebrate genome it holds);
- how many species carry an orthologue and how many are one-to-one;
- its paralogues (within-species and older), which are the copy-and-paste edges the
  knowledge graph lacked.

Each species is placed on the ladder once through Ensembl's taxonomy classification
(`species_strata.json`, a cache). Evidence is `curated` (Ensembl Compara gene trees); an
origin is a lower bound on age, since a gene lost in every sampled outgroup reads younger
than it is, and Ensembl vertebrates samples three non-animal or invertebrate outgroups.
"""

from __future__ import annotations

import csv
import gzip
import io
import json
import time
import urllib.request
from collections import Counter
from pathlib import Path
from typing import Any

from genomeos.results import RESULTS_DIR, load_result, save_result

RELEASE = 116
DUMP_URL = (
    "https://ftp.ensembl.org/pub/release-{release}/tsv/ensembl-compara/homologies/homo_sapiens/"
    "Compara.{release}.protein_default.homologies.tsv.gz"
)
KNOWLEDGE = Path("data/knowledge/homology")
GENES_TSV = Path("data/cache/gencode_genes.tsv")
# deep to shallow: the clade ladder every species is placed on, relative to human
STRATA = (
    "Eukaryota",
    "Opisthokonta",
    "Metazoa",
    "Bilateria",
    "Chordata",
    "Vertebrata",
    "Gnathostomata",
    "Euteleostomi",
    "Sarcopterygii",
    "Tetrapoda",
    "Amniota",
    "Mammalia",
    "Theria",
    "Eutheria",
    "Boreoeutheria",
    "Euarchontoglires",
    "Primates",
    "Haplorrhini",
    "Simiiformes",
    "Catarrhini",
    "Hominoidea",
    "Hominidae",
    "Homininae",
    "Homo",
)
RANK = {s: i for i, s in enumerate(STRATA)}
# the conversation's library ladder, as the stratum each name begins at
LADDER = {
    "Eukaryota": "Life/Eukaryote core",
    "Opisthokonta": "Eukaryote core",
    "Metazoa": "Animal core",
    "Bilateria": "Animal core",
    "Chordata": "Chordate",
    "Vertebrata": "Vertebrate core",
    "Gnathostomata": "Vertebrate core",
    "Euteleostomi": "Vertebrate core",
    "Sarcopterygii": "Tetrapod",
    "Tetrapoda": "Tetrapod",
    "Amniota": "Amniote",
    "Mammalia": "Mammal core",
    "Theria": "Mammal core",
    "Eutheria": "Mammal core",
    "Boreoeutheria": "Mammal core",
    "Euarchontoglires": "Mammal core",
    "Primates": "Primate",
    "Haplorrhini": "Primate",
    "Simiiformes": "Primate",
    "Catarrhini": "Primate",
    "Hominoidea": "Ape",
    "Hominidae": "Ape",
    "Homininae": "Ape",
    "Homo": "Human-specific (in this species set)",
}
# Ensembl's classification lists only some nodes (no Amniota, Tetrapoda, Theria, Boreoeutheria ...), so
# these node names stand in for the stratum they imply relative to human
PROXIES = {
    "Fungi": "Opisthokonta",
    "Ascomycota": "Opisthokonta",
    "Saccharomycetes": "Opisthokonta",
    "Choanoflagellata": "Opisthokonta",
    "Protostomia": "Bilateria",
    "Ecdysozoa": "Bilateria",
    "Arthropoda": "Bilateria",
    "Nematoda": "Bilateria",
    "Lophotrochozoa": "Bilateria",
    "Tunicata": "Chordata",
    "Urochordata": "Chordata",
    "Cephalochordata": "Chordata",
    "Cyclostomata": "Vertebrata",
    "Myxini": "Vertebrata",
    "Hyperoartia": "Vertebrata",
    "Petromyzontiformes": "Vertebrata",
    "Chondrichthyes": "Gnathostomata",
    "Actinopterygii": "Euteleostomi",
    "Coelacanthiformes": "Sarcopterygii",
    "Coelacanthimorpha": "Sarcopterygii",
    "Dipnoi": "Sarcopterygii",
    "Dipnomorpha": "Sarcopterygii",
    "Amphibia": "Tetrapoda",
    "Batrachia": "Tetrapoda",
    "Sauropsida": "Amniota",
    "Aves": "Amniota",
    "Archelosauria": "Amniota",
    "Archosauria": "Amniota",
    "Testudines": "Amniota",
    "Crocodylia": "Amniota",
    "Lepidosauria": "Amniota",
    "Squamata": "Amniota",
    "Monotremata": "Mammalia",
    "Prototheria": "Mammalia",
    "Metatheria": "Theria",
    "Marsupialia": "Theria",
    "Didelphimorphia": "Theria",
    "Diprotodontia": "Theria",
    "Afrotheria": "Eutheria",
    "Xenarthra": "Eutheria",
    "Laurasiatheria": "Boreoeutheria",
    "Glires": "Euarchontoglires",
    "Rodentia": "Euarchontoglires",
    "Lagomorpha": "Euarchontoglires",
    "Scandentia": "Euarchontoglires",
    "Dermoptera": "Euarchontoglires",
    "Strepsirrhini": "Primates",
    "Tarsiiformes": "Primates",
    "Platyrrhini": "Simiiformes",
    "Cercopithecoidea": "Catarrhini",
    "Cercopithecidae": "Catarrhini",
    "Hylobatidae": "Hominoidea",
    "Ponginae": "Hominidae",
    "Pongo": "Hominidae",
    "Gorilla": "Homininae",
    "Pan": "Homininae",
}


def place(names: set[str]) -> str | None:
    """A species' stratum: the shallowest clade shared with human that its classification implies."""
    found = [RANK[n] for n in names if n in RANK and n != "Homo"]
    found += [RANK[PROXIES[n]] for n in names if n in PROXIES]
    return STRATA[max(found)] if found else None


def classify_species(
    species: list[str], cache: Path = KNOWLEDGE / "species_strata.json", progress=None
) -> dict:
    """Place every species on the ladder through Ensembl's taxonomy classification; resumable, cached."""
    done = json.loads(cache.read_text()) if cache.exists() else {}
    cache.parent.mkdir(parents=True, exist_ok=True)
    for i, sp in enumerate(sorted(species)):
        if sp in done and done[sp].get("stratum"):
            continue
        try:
            nodes = _rest(f"https://rest.ensembl.org/taxonomy/classification/{sp}")
            names = {x.get("scientific_name") for x in nodes} | {x.get("name") for x in nodes}
            done[sp] = {"stratum": place(names), "nodes": sorted(n for n in names if n)}
        except Exception as e:  # noqa: BLE001 - one species failing must not stop the pass
            done[sp] = {"stratum": None, "error": str(e)[:100]}
        if i % 10 == 0:
            cache.write_text(json.dumps(done, indent=0, sort_keys=True))
            if progress:
                progress(i, len(species), sp, done[sp].get("stratum"))
        time.sleep(0.2)
    cache.write_text(json.dumps(done, indent=0, sort_keys=True))
    return done


def ensembl_species() -> list[str]:
    return [
        s["name"]
        for s in _rest("https://rest.ensembl.org/info/species?division=EnsemblVertebrates")["species"]
    ]


def _rest(url: str, tries: int = 4):
    for i in range(tries):
        try:
            req = urllib.request.Request(
                url, headers={"Content-Type": "application/json", "User-Agent": "GenomeOS/0.9 (homology)"}
            )
            with urllib.request.urlopen(req, timeout=90) as r:  # noqa: S310
                return json.load(r)
        except OSError:
            if i == tries - 1:
                raise
            time.sleep(3 * (i + 1))
    return None


ORTHOLOG_TYPES = {"ortholog_one2one", "ortholog_one2many", "ortholog_many2many"}
PARALOG_TYPES = {"within_species_paralog", "other_paralog", "gene_split"}
EVIDENCE = {
    "origin": (
        f"curated: Ensembl Compara release {RELEASE} protein gene trees; the deepest clade with an orthologue"
    ),
    "strata": "curated: Ensembl/NCBI taxonomy classification per species, placed on a fixed clade ladder",
    "paralogues": f"curated: Ensembl Compara release {RELEASE} within-species and older paralogues",
    "caveat": (
        "an origin is a lower bound on age: a gene lost in every sampled outgroup reads younger than it is"
    ),
}


def deepest(strata: list[str | None]) -> str | None:
    """The deepest stratum among those given (the smallest rank); None if none is placed."""
    ranks = [RANK[s] for s in strata if s in RANK]
    return STRATA[min(ranks)] if ranks else None


def load_species_strata(path: Path = KNOWLEDGE / "species_strata.json") -> dict[str, str | None]:
    if not path.exists():
        raise FileNotFoundError(
            f"{path} missing: classify Ensembl's species first (scripts/origin_genome_wide.py)"
        )
    d = json.loads(path.read_text())
    return {sp: v.get("stratum") for sp, v in d.items()}


def load_symbols(path: Path = GENES_TSV) -> dict[str, tuple[str, str]]:
    """Ensembl gene id (no version) → (symbol, gene_type) from the cached GENCODE gene table."""
    out: dict[str, tuple[str, str]] = {}
    if not path.exists():
        return out
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[row["gene_id"].split(".")[0]] = (row["symbol"], row["gene_type"])
    return out


def _open_dump(url: str):
    if url.startswith(("http://", "https://")):
        req = urllib.request.Request(url, headers={"User-Agent": "GenomeOS/0.9 (homology stream)"})
        resp = urllib.request.urlopen(req, timeout=900)  # noqa: S310
        return gzip.open(io.BufferedReader(resp, 1 << 20), "rt")
    return gzip.open(url, "rt") if url.endswith(".gz") else open(url)


def stream(url: str, strata: dict[str, str | None], progress=None) -> tuple[dict[str, dict], dict]:
    """One pass over the dump: per human gene, its deepest orthologue stratum, counts and paralogues."""
    t0 = time.time()
    genes: dict[str, dict] = {}
    unplaced: Counter = Counter()
    rows = 0
    with _open_dump(url) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        col = {name: i for i, name in enumerate(header)}
        gi, ht, hg, hs = (
            col["gene_stable_id"],
            col["homology_type"],
            col["homology_gene_stable_id"],
            col["homology_species"],
        )
        for r in reader:
            rows += 1
            g = genes.setdefault(
                r[gi],
                {
                    "rank": None,
                    "species": set(),
                    "one2one": 0,
                    "paralogues": set(),
                    "paralogue_types": Counter(),
                },
            )
            typ = r[ht]
            if typ in ORTHOLOG_TYPES:
                sp = r[hs]
                g["species"].add(sp)
                if typ == "ortholog_one2one":
                    g["one2one"] += 1
                s = strata.get(sp)
                if s is None:
                    unplaced[sp] += 1
                else:
                    rk = RANK[s]
                    if g["rank"] is None or rk < g["rank"]:
                        g["rank"] = rk
            elif typ in PARALOG_TYPES and r[hs] == "homo_sapiens":
                g["paralogues"].add(r[hg])
                g["paralogue_types"][typ] += 1
            if progress and rows % 500_000 == 0:
                progress(rows, len(genes), time.time() - t0)
    cost = {
        "rows": rows,
        "seconds": round(time.time() - t0, 1),
        "unplaced_species": dict(unplaced.most_common(20)),
    }
    return genes, cost


def distil(
    genes: dict[str, dict], symbols: dict[str, tuple[str, str]], members: dict[str, list[str]] | None
) -> dict:
    """Per gene the origin and paralogues by symbol; per stratum the counts; per library the age."""
    per_gene: dict[str, dict] = {}
    by_stratum: Counter = Counter()
    coding_only = 0
    for gid, g in genes.items():
        sym, gtype = symbols.get(gid, (gid, "unknown"))
        if gtype not in ("protein_coding", "unknown"):
            continue
        coding_only += 1
        origin = STRATA[g["rank"]] if g["rank"] is not None else "Homo"
        by_stratum[origin] += 1
        per_gene[sym] = {
            "origin": origin,
            "ladder": LADDER[origin],
            "species": len(g["species"]),
            "one2one": g["one2one"],
            "paralogues": sorted(symbols.get(p, (p, ""))[0] for p in g["paralogues"]),
            "paralogue_types": dict(g["paralogue_types"]),
        }
    libraries: dict[str, dict] = {}
    if members:
        for lib, syms in members.items():
            placed = [per_gene[s]["origin"] for s in syms if s in per_gene]
            if not placed:
                continue
            dist = Counter(placed)
            ordered = sorted(dist.items(), key=lambda kv: RANK[kv[0]])
            cum, median = 0, None
            for s, n in ordered:
                cum += n
                if median is None and cum * 2 >= len(placed):
                    median = s
            deep = sum(n for s, n in dist.items() if RANK[s] <= RANK["Bilateria"])
            vert = sum(n for s, n in dist.items() if RANK["Chordata"] <= RANK[s] <= RANK["Tetrapoda"])
            mamm = sum(n for s, n in dist.items() if RANK[s] >= RANK["Amniota"])
            libraries[lib] = {
                "members_placed": len(placed),
                "members": len(syms),
                "median_origin": median,
                "ladder": LADDER[median] if median else None,
                "deepest": ordered[0][0],
                "share_animal_or_older": round(deep / len(placed), 3),
                "share_vertebrate_to_tetrapod": round(vert / len(placed), 3),
                "share_amniote_or_younger": round(mamm / len(placed), 3),
                "distribution": {s: n for s, n in ordered},
            }
    return {
        "genes": per_gene,
        "coding_genes": coding_only,
        "by_stratum": {s: by_stratum[s] for s in STRATA if by_stratum[s]},
        "by_ladder": dict(Counter(LADDER[s] for s in by_stratum.elements())),
        "libraries": libraries,
        "strata": list(STRATA),
        "evidence": EVIDENCE,
    }


def library_members() -> dict[str, list[str]] | None:
    try:
        from genomeos.lib.membership import KnowledgeBase

        d = KnowledgeBase.distilled()
        return d.get("members") if d else None
    except Exception:  # pragma: no cover - optional
        return None


def run_and_save(
    url: str | None = None, results_dir: Path = RESULTS_DIR, progress=None, name: str = "origin_genome_wide"
) -> dict:
    strata = load_species_strata()
    genes, cost = stream(url or DUMP_URL.format(release=RELEASE), strata, progress)
    out = distil(genes, load_symbols(), library_members())
    out["cost"] = cost
    out["species_placed"] = sum(1 for v in strata.values() if v)
    out["species_unplaced"] = sorted(sp for sp, v in strata.items() if not v)
    out["source"] = url or DUMP_URL.format(release=RELEASE)
    save_result(name, out, results_dir)
    return out


def origin_of(symbol: str, results_dir: Path = RESULTS_DIR) -> dict[str, Any] | None:
    r = load_result("origin_genome_wide", results_dir) or {}
    return (r.get("genes") or {}).get(symbol)
