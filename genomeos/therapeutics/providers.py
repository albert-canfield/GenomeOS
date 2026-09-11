"""Provider abstractions over the biological databases.

Business logic never talks to an API. It asks a provider, and a provider
either answers with data carrying its provenance or answers that the data is
unavailable and why. That keeps every stage testable offline, makes a source
replaceable without touching the reasoning, and stops a network failure from
turning into a fabricated value.

Wired sources, all public and key-free:

    UniProt / Ensembl / InterPro / PDB / AlphaFold / Reactome / STRING / HPA
        through the existing federated protein compiler (cached on disk)
    Human Protein Atlas search API   per-tissue normal RNA levels (nTPM)
    Open Targets Platform GraphQL    tractability, clinical precedent, safety
    Ensembl comparative genomics     human paralogues for the negative set

Nothing here scrapes a page that has an API, and every response is cached
under data/knowledge/ so an analysis is reproducible without the network.
"""

from __future__ import annotations

import json
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Protocol, runtime_checkable

from .evidence import Evidence, database, today

UA = {"Accept": "application/json", "User-Agent": "GenomeOS/0.1 (therapeutics)"}
CACHE = Path("data/knowledge/therapeutics")

HPA_API = "https://www.proteinatlas.org/api/search_download.php"
OPENTARGETS_API = "https://api.platform.opentargets.org/api/v4/graphql"
ENSEMBL_HOMOLOGY = "https://rest.ensembl.org/homology/id/homo_sapiens"
ENSEMBL_LOOKUP = "https://rest.ensembl.org/lookup/id"
ENSEMBL_SEQUENCE = "https://rest.ensembl.org/sequence/id"

#: HPA consensus-tissue columns, with the weight an on-target hit there carries.
#: Weights are a stated policy of this pipeline, not a measurement: damage to a
#: tissue without a spare and without a transplant option is weighted hardest.
CRITICAL_TISSUES: dict[str, tuple[str, float]] = {
    "t_RNA_heart_muscle": ("heart muscle", 1.0),
    "t_RNA_cerebral_cortex": ("cerebral cortex", 1.0),
    "t_RNA_cerebellum": ("cerebellum", 1.0),
    "t_RNA_hypothalamus": ("hypothalamus", 1.0),
    "t_RNA_bone_marrow": ("bone marrow", 1.0),
    "t_RNA_lung": ("lung", 0.9),
    "t_RNA_liver": ("liver", 0.9),
    "t_RNA_kidney": ("kidney", 0.9),
    "t_RNA_adrenal_gland": ("adrenal gland", 0.7),
    "t_RNA_pancreas": ("pancreas", 0.7),
    "t_RNA_skeletal_muscle": ("skeletal muscle", 0.6),
    "t_RNA_spleen": ("spleen", 0.5),
    "t_RNA_thyroid_gland": ("thyroid gland", 0.5),
    "t_RNA_colon": ("colon", 0.5),
    "t_RNA_stomach_1": ("stomach", 0.5),
    "t_RNA_esophagus": ("esophagus", 0.4),
    "t_RNA_skin_1": ("skin", 0.4),
    "t_RNA_retina": ("retina", 0.6),
    "t_RNA_testis": ("testis", 0.2),
    "t_RNA_placenta": ("placenta", 0.1),
}

HPA_EXTRA_COLUMNS = ("g", "eg", "scl", "scml", "rnats", "rnatd", "rnatsm", "pc", "secretome_location")

HPA_LICENCE = "Human Protein Atlas, CC BY-SA 4.0 (proteinatlas.org)"
OT_LICENCE = "Open Targets Platform, CC0 1.0 (platform.opentargets.org)"


# --- shared plumbing ------------------------------------------------------------------


class UnavailableError(Exception):
    """Raised inside a provider; callers turn it into an explicit absence."""


@dataclass(slots=True)
class Answer:
    """A provider response: data with provenance, or an explained absence."""

    data: Any = None
    evidence: list[Evidence] = field(default_factory=list)
    available: bool = True
    reason: str = ""
    retrieved_at: str = ""

    @classmethod
    def none(cls, reason: str) -> Answer:
        return cls(None, [], False, reason, today())

    def to_dict(self) -> dict[str, Any]:
        return {
            "available": self.available,
            "reason": self.reason,
            "retrieved_at": self.retrieved_at,
            "evidence": [e.to_dict() for e in self.evidence],
        }


@dataclass
class DiskCache:
    """One JSON file per key under data/knowledge/therapeutics/<namespace>/."""

    namespace: str
    root: Path = CACHE

    def path(self, key: str) -> Path:
        safe = urllib.parse.quote(key, safe="")[:120]
        return self.root / self.namespace / f"{safe}.json"

    def get(self, key: str) -> Any | None:
        p = self.path(key)
        if not p.exists():
            return None
        try:
            return json.loads(p.read_text())
        except json.JSONDecodeError:
            return None

    def put(self, key: str, value: Any) -> None:
        p = self.path(key)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(value))


def _fetch_json(url: str, timeout: int = 45, data: bytes | None = None, retries: int = 2) -> Any:
    """GET or POST JSON, backing off on rate limits and server hiccups."""
    delay = 1.5
    last: Exception | None = None
    for attempt in range(retries + 1):
        req = urllib.request.Request(url, data=data, headers={**UA, "Content-Type": "application/json"})
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:  # noqa: S310
                return json.load(r)
        except urllib.error.HTTPError as e:
            last = e
            if attempt == retries or e.code not in (429, 500, 502, 503, 504):
                break
            time.sleep(float(e.headers.get("Retry-After") or delay))
            delay *= 2
        except (urllib.error.URLError, TimeoutError, OSError, ValueError) as e:
            last = e
            if attempt == retries:
                break
            time.sleep(delay)
            delay *= 2
    raise UnavailableError(f"{type(last).__name__}: {str(last)[:120]}") from last


# --- interfaces -----------------------------------------------------------------------


@runtime_checkable
class ProteinAnnotationProvider(Protocol):
    """Sequence, subcellular location, topology, domains, structures, pathways."""

    def annotation(self, gene: str) -> Answer: ...


@runtime_checkable
class ExpressionProvider(Protocol):
    """Where a gene is expressed in healthy tissue, and how strongly."""

    def normal_expression(self, gene: str) -> Answer: ...


@runtime_checkable
class TumourExpressionProvider(Protocol):
    """Expression in this patient's tumour. Population data must not stand in."""

    def tumour_expression(self, gene: str) -> Answer: ...


@runtime_checkable
class CohortExpressionProvider(Protocol):
    """How a gene behaves across a cohort of the same cancer type."""

    def cohort(self, gene: str) -> Answer: ...


@runtime_checkable
class StructureProvider(Protocol):
    """Experimental and predicted structures for a protein."""

    def structures(self, gene: str, annotation: dict[str, Any] | None = None) -> Answer: ...


@runtime_checkable
class TraffickingProvider(Protocol):
    """Endocytosis, endosomal and lysosomal routing, recycling, shedding."""

    def trafficking(self, gene: str, annotation: dict[str, Any] | None = None) -> Answer: ...


@runtime_checkable
class TranscriptSequenceProvider(Protocol):
    """The coding sequence of one transcript, for reconstructing an altered protein."""

    def cds(self, transcript: str) -> Answer: ...


@runtime_checkable
class NeoantigenProvider(Protocol):
    """Peptide/HLA binding and processing prediction."""

    @property
    def name(self) -> str: ...

    def predict(self, peptides: list[str], alleles: list[str]) -> Answer: ...


@runtime_checkable
class CancerEvidenceProvider(Protocol):
    """Therapeutic precedent: approved drugs, trials, tractability, liabilities."""

    def precedent(self, gene: str, ensembl_gene: str | None = None) -> Answer: ...


@runtime_checkable
class HomologyProvider(Protocol):
    """Closest human relatives of a protein, for the negative-target set."""

    def paralogues(self, gene: str, ensembl_gene: str | None = None) -> Answer: ...


# --- implementations ------------------------------------------------------------------


@dataclass
class CompiledProteinProvider:
    """The federated protein compiler as a provider.

    Reads data/knowledge/proteins/<GENE>.json when present; compiles it from
    UniProt, Ensembl, InterPro, PDB, AlphaFold, Reactome, STRING and HPA when
    the network is allowed. The compiler already attaches an evidence grade
    and a confidence to every section, so those are passed through unchanged.
    """

    net: bool = True
    cache_dir: Path = Path("data/knowledge/proteins")

    def annotation(self, gene: str) -> Answer:
        gene = gene.upper()
        cached = self.cache_dir / f"{gene}.json"
        defn: dict[str, Any] | None = None
        if cached.exists():
            try:
                defn = json.loads(cached.read_text())
            except json.JSONDecodeError:
                defn = None
        if defn is None:
            if not self.net:
                return Answer.none(f"no compiled protein definition for {gene} and network disabled")
            try:
                from genomeos.molecules import compile_protein

                defn = compile_protein(gene, cache_dir=self.cache_dir)
            except Exception as e:  # noqa: BLE001 - any source failure is an absence, not a crash
                return Answer.none(f"protein compilation failed for {gene}: {str(e)[:120]}")
        ident = (defn.get("sections", {}).get("identity") or {}).get("items") or {}
        if not ident:
            return Answer.none(f"no reviewed UniProt entry for {gene}")
        ev = [
            database(
                "UniProtKB/Swiss-Prot",
                f"{gene} is {ident.get('name')} ({ident.get('length')} aa), accession "
                f"{ident.get('accession')}",
                0.95,
                "human",
                identifier=ident.get("accession"),
                url=f"https://www.uniprot.org/uniprotkb/{ident.get('accession')}",
            )
        ]
        return Answer(defn, ev, True, "", today())


@dataclass
class HpaExpressionProvider:
    """Normal-tissue RNA levels from the Human Protein Atlas search API.

    HPA publishes consensus nTPM per tissue. Those are population-level
    healthy-tissue measurements, never this patient's tumour, and they are
    labelled as such wherever they are used.
    """

    net: bool = True
    cache: DiskCache = field(default_factory=lambda: DiskCache("hpa"))

    def normal_expression(self, gene: str) -> Answer:
        gene = gene.upper()
        row = self.cache.get(gene)
        if row is None:
            if not self.net:
                return Answer.none(f"no cached HPA record for {gene} and network disabled")
            cols = ",".join([*HPA_EXTRA_COLUMNS, *CRITICAL_TISSUES])
            url = f"{HPA_API}?search={urllib.parse.quote(gene)}&format=json&columns={cols}&compress=no"
            try:
                rows = _fetch_json(url)
            except UnavailableError as e:
                return Answer.none(f"Human Protein Atlas unavailable: {e}")
            match = next((r for r in rows or [] if str(r.get("Gene", "")).upper() == gene), None)
            if match is None:
                return Answer.none(f"no Human Protein Atlas record for {gene}")
            row = match
            self.cache.put(gene, row)
        ev = [
            database(
                "Human Protein Atlas (consensus tissue RNA)",
                f"normal-tissue RNA levels for {gene} across {len(CRITICAL_TISSUES)} queried tissues",
                0.8,
                "human",
                identifier=row.get("Ensembl"),
                url=f"https://www.proteinatlas.org/{row.get('Ensembl', '')}",
            )
        ]
        return Answer(row, ev, True, "", today())


@dataclass
class PatientRnaProvider:
    """Patient tumour RNA-seq, supplied as a two-column gene/value table.

    When present this overrides every population-level assumption about
    tumour expression. When absent the answer is an explicit absence: DNA
    evidence alone never establishes that a protein is made.
    """

    values: dict[str, float] = field(default_factory=dict)
    unit: str = "TPM"
    sample_id: str = ""
    path: str = ""

    @classmethod
    def from_file(cls, path: str | Path, unit: str = "TPM", sample_id: str = "") -> PatientRnaProvider:
        values: dict[str, float] = {}
        p = Path(path)
        for line in p.read_text().splitlines():
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            parts = line.replace(",", "\t").split("\t")
            if len(parts) < 2:
                continue
            try:
                values[parts[0].strip().upper()] = float(parts[1])
            except ValueError:
                continue  # header row
        return cls(values, unit, sample_id or p.stem, str(p))

    @property
    def available(self) -> bool:
        return bool(self.values)

    @property
    def whole_transcriptome(self) -> bool:
        """Enough genes for a within-sample rank to mean anything."""
        return len(self.values) >= 5000

    def rank(self, gene: str) -> float | None:
        """Where this gene sits among all genes measured in this tumour.

        A rank inside one sample needs no unit, so it survives the comparison a
        raw value cannot: a patient's TPM and a cohort's RSEM are different
        scales and must never be compared directly.
        """
        if not self.whole_transcriptome:
            return None
        v = self.values.get(gene.upper())
        if v is None:
            return None
        below = sum(1 for x in self.values.values() if x < v)
        return round(below / len(self.values), 4)

    def tumour_expression(self, gene: str) -> Answer:
        if not self.values:
            return Answer.none("no tumour RNA-seq supplied")
        v = self.values.get(gene.upper())
        if v is None:
            return Answer.none(f"{gene} absent from the supplied tumour RNA-seq table")
        from .evidence import patient

        rank = self.rank(gene)
        claim = f"{gene} measured at {v:g} {self.unit} in this tumour"
        if rank is not None:
            claim += f", at the {rank:.0%} rank of the {len(self.values):,} genes measured in this sample"
        ev = [patient(f"tumour RNA-seq ({self.path or 'supplied table'})", claim, 0.9)]
        return Answer(
            {"value": v, "unit": self.unit, "sample": self.sample_id, "rank": rank},
            ev,
            True,
            "",
            today(),
        )


@dataclass
class CBioPortalCohortProvider:
    """Cohort expression for one cancer type, from cBioPortal.

    A cohort cannot say what this patient's tumour does. It answers a different
    and still useful question: in this cancer type, is the gene raised above
    normal tissue often, rarely, or never? The profile preferred is the one
    scored against the study's own normal samples, which is unit-free and
    therefore comparable; a raw expression profile is used only if no z-score
    profile exists, and is labelled as not comparable across datasets.
    """

    study: str = ""
    net: bool = True
    cache: DiskCache = field(default_factory=lambda: DiskCache("cohort"))
    _warmed: set[str] = field(default_factory=set)

    @property
    def available(self) -> bool:
        return bool(self.study)

    def _key(self, gene: str) -> str:
        return f"{self.study}|{gene.upper()}"

    def warm(self, genes: list[str]) -> None:
        """Fetch a whole gene list in one call instead of one call per gene."""
        if not self.study or not self.net:
            return
        wanted = [g.upper() for g in genes if g and self.cache.get(self._key(g)) is None]
        if not wanted:
            return
        try:
            from genomeos.cancer.cbioportal import CBioPortal

            d = CBioPortal(timeout=180).expression_distribution(self.study, wanted)
        except Exception as e:  # noqa: BLE001 - a cohort is optional context, never a failure
            self._warmed.update(wanted)
            self.cache.put(f"{self.study}|__error__", {"error": str(e)[:160]})
            return
        meta = {k: d[k] for k in ("study", "profile", "kind", "meaning", "samples")}
        for gene, summary in d["genes"].items():
            self.cache.put(self._key(gene), {**meta, "gene": gene, "summary": summary})
        for gene in d["missing"]:
            self.cache.put(self._key(gene), {**meta, "gene": gene, "summary": None})
        self._warmed.update(wanted)

    def cohort(self, gene: str) -> Answer:
        if not self.study:
            return Answer.none("no cohort study selected; pass one to compare this tumour against")
        row = self.cache.get(self._key(gene))
        if row is None and gene.upper() not in self._warmed:
            self.warm([gene])
            row = self.cache.get(self._key(gene))
        if row is None:
            err = self.cache.get(f"{self.study}|__error__") or {}
            return Answer.none(
                f"cohort expression unavailable for {gene}: {err.get('error', 'not retrieved')}"
            )
        if not row.get("summary"):
            return Answer.none(f"{gene} is not measured in the {self.study} expression profile")
        s = row["summary"]
        claim = (
            f"across {s['n']} tumours of {row['study']}, {gene} has median "
            f"{s['median']:+.2f} ({row['meaning']})"
        )
        if "fraction_raised" in s:
            claim += (
                f"; raised above normal tissue (z >= {s['raised_threshold']:g}) in "
                f"{s['fraction_raised']:.1%} of them"
            )
        ev = [
            database(
                f"cBioPortal {row['study']}",
                claim,
                0.8,
                "human",
                identifier=row["profile"],
                url=f"https://www.cbioportal.org/study/summary?id={row['study']}",
            )
        ]
        return Answer(row, ev, True, "", today())


@dataclass
class CompiledStructureProvider:
    """PDB and AlphaFold entries, read from the compiled protein definition."""

    protein: CompiledProteinProvider

    def structures(self, gene: str, annotation: dict[str, Any] | None = None) -> Answer:
        defn = annotation
        if defn is None:
            a = self.protein.annotation(gene)
            if not a.available:
                return Answer.none(a.reason)
            defn = a.data
        sec = defn.get("sections", {})
        exp = (sec.get("structures_experimental") or {}).get("items") or []
        pred = (sec.get("structures_predicted") or {}).get("items") or []
        if not exp and not pred:
            return Answer.none(f"no experimental or predicted structure listed for {gene}")
        ev: list[Evidence] = []
        if exp:
            best = min(
                (x for x in exp if x.get("resolution_A")),
                key=lambda x: x["resolution_A"],
                default=exp[0],
            )
            ev.append(
                database(
                    "PDB (via UniProt cross-references)",
                    f"{len(exp)} experimental structures for {gene}; best {best.get('id')} "
                    f"({best.get('method')}"
                    + (f", {best['resolution_A']} A" if best.get("resolution_A") else "")
                    + ")",
                    0.95,
                    "in_vitro",
                    identifier=best.get("id"),
                    url=f"https://www.ebi.ac.uk/pdbe/entry/pdb/{str(best.get('id', '')).lower()}",
                )
            )
        if pred:
            ev.append(
                Evidence(
                    "AlphaFold DB",
                    "prediction",
                    f"predicted structure available for {gene} ({pred[0].get('id')}); pLDDT per residue",
                    "computational",
                    0.7,
                    identifier=pred[0].get("id"),
                    url=f"https://alphafold.ebi.ac.uk/entry/{pred[0].get('id')}",
                    retrieved_at=today(),
                )
            )
        return Answer({"experimental": exp, "predicted": pred}, ev, True, "", today())


@dataclass
class AnnotationTraffickingProvider:
    """Trafficking read out of curated annotation, never guessed.

    UniProt keywords and subcellular-location terms record endocytosis,
    endosomal and lysosomal routing, recycling and shedding where they have
    been curated. Gene Ontology adds the same facts from a second source.
    Absence of a keyword is absence of evidence, and is reported as unknown.
    """

    protein: CompiledProteinProvider
    knowledge_base: Any = None

    def trafficking(self, gene: str, annotation: dict[str, Any] | None = None) -> Answer:
        defn = annotation
        if defn is None:
            a = self.protein.annotation(gene)
            if not a.available:
                return Answer.none(a.reason)
            defn = a.data
        sec = defn.get("sections", {})
        ident = (sec.get("identity") or {}).get("items") or {}
        fn = (sec.get("function") or {}).get("items") or {}
        payload = {
            "keywords": [k.lower() for k in ident.get("keywords") or []],
            "locations": [x.lower() for x in fn.get("location") or []],
            "processing": (sec.get("processing") or {}).get("items") or [],
            "modifications": (sec.get("modifications") or {}).get("items") or [],
            "go_terms": self._go_terms(gene),
        }
        ev = [
            database(
                "UniProtKB/Swiss-Prot keywords and subcellular locations",
                f"curated trafficking annotation for {gene}",
                0.85,
                "human",
                identifier=ident.get("accession"),
            )
        ]
        return Answer(payload, ev, True, "", today())

    def _go_terms(self, gene: str) -> list[str]:
        kb = self.knowledge_base
        if kb is None:
            return []
        try:
            return sorted(kb.annotations.terms(gene))
        except Exception:  # noqa: BLE001 - the GO files are optional
            return []


@dataclass
class OpenTargetsProvider:
    """Therapeutic precedent and tractability from the Open Targets Platform.

    Gives, per target: which modalities have reached the clinic, the drugs
    and their stages, and curated safety liabilities. An approved antibody
    against a gene is clinical evidence that the gene is reachable by an
    antibody. It is not evidence that the drug suits this patient's tumour.
    """

    net: bool = True
    cache: DiskCache = field(default_factory=lambda: DiskCache("opentargets"))

    QUERY = """query T($id: String!) {
      target(ensemblId: $id) {
        id approvedSymbol approvedName
        tractability { label modality value }
        safetyLiabilities { event datasource }
        subcellularLocations { location source }
        drugAndClinicalCandidates {
          count
          rows { maxClinicalStage drug { id name drugType } diseases { disease { name } } }
        }
      }
    }"""

    def precedent(self, gene: str, ensembl_gene: str | None = None) -> Answer:
        if not ensembl_gene:
            return Answer.none(f"no Ensembl gene id for {gene}; Open Targets is keyed by Ensembl id")
        payload = self.cache.get(ensembl_gene)
        if payload is None:
            if not self.net:
                return Answer.none(f"no cached Open Targets record for {gene} and network disabled")
            body = json.dumps({"query": self.QUERY, "variables": {"id": ensembl_gene}}).encode()
            try:
                res = _fetch_json(OPENTARGETS_API, data=body)
            except UnavailableError as e:
                return Answer.none(f"Open Targets unavailable: {e}")
            payload = (res or {}).get("data", {}).get("target")
            if not payload:
                return Answer.none(f"no Open Targets record for {gene} ({ensembl_gene})")
            self.cache.put(ensembl_gene, payload)
        return Answer(payload, self._evidence(gene, payload, ensembl_gene), True, "", today())

    @staticmethod
    def _evidence(gene: str, payload: dict[str, Any], ensembl_gene: str) -> list[Evidence]:
        from .evidence import clinical

        url = f"https://platform.opentargets.org/target/{ensembl_gene}"
        ev: list[Evidence] = []
        rows = ((payload.get("drugAndClinicalCandidates") or {}).get("rows")) or []
        approved = [r for r in rows if r.get("maxClinicalStage") == "APPROVAL"]
        for kind in ("Antibody", "Antibody drug conjugate", "Small molecule", "Protein", "Cell therapy"):
            hits = [r for r in approved if (r.get("drug") or {}).get("drugType") == kind]
            if hits:
                names = ", ".join(sorted({(r["drug"]["name"] or "").title() for r in hits})[:4])
                ev.append(
                    clinical(
                        "Open Targets Platform / ChEMBL",
                        f"{len(hits)} approved {kind.lower()} therapies target {gene}: {names}",
                        0.95,
                        identifier=ensembl_gene,
                        url=url,
                    )
                )
        in_trials = [r for r in rows if str(r.get("maxClinicalStage", "")).startswith("PHASE")]
        if in_trials:
            ev.append(
                Evidence(
                    "Open Targets Platform / ChEMBL",
                    "clinical_trial",
                    f"{len(in_trials)} agents against {gene} have entered clinical trials",
                    "clinical",
                    0.85,
                    ensembl_gene,
                    url,
                    today(),
                )
            )
        ab = [t for t in payload.get("tractability") or [] if t.get("modality") == "AB" and t.get("value")]
        if ab:
            ev.append(
                database(
                    "Open Targets tractability (antibody modality)",
                    f"{gene} antibody tractability buckets: " + ", ".join(t["label"] for t in ab),
                    0.8,
                    "human",
                    identifier=ensembl_gene,
                    url=url,
                )
            )
        for s in (payload.get("safetyLiabilities") or [])[:6]:
            ev.append(
                Evidence(
                    f"Open Targets safety ({s.get('datasource')})",
                    "publication",
                    f"curated safety liability for {gene}: {s.get('event')}",
                    "human",
                    0.7,
                    ensembl_gene,
                    url,
                    today(),
                )
            )
        return ev


@dataclass
class EnsemblParalogueProvider:
    """Closest human relatives, with sequence identity, from Ensembl Compara.

    These are the proteins a binder is most likely to cross-react with, so
    they belong in the negative-target set with their identity attached.
    """

    net: bool = True
    cache: DiskCache = field(default_factory=lambda: DiskCache("paralogues"))
    symbols: DiskCache = field(default_factory=lambda: DiskCache("ensembl_symbols"))

    def paralogues(self, gene: str, ensembl_gene: str | None = None) -> Answer:
        if not ensembl_gene:
            return Answer.none(f"no Ensembl gene id for {gene}; paralogue lookup is keyed by Ensembl id")
        rows = self.cache.get(ensembl_gene)
        if rows is None:
            if not self.net:
                return Answer.none(f"no cached paralogue record for {gene} and network disabled")
            url = (
                f"{ENSEMBL_HOMOLOGY}/{ensembl_gene}?type=paralogues;target_species=homo_sapiens;"
                "sequence=none;content-type=application/json"
            )
            try:
                res = _fetch_json(url, timeout=60)
            except UnavailableError as e:
                return Answer.none(f"Ensembl comparative genomics unavailable: {e}")
            data = (res or {}).get("data") or []
            homologies = data[0].get("homologies", []) if data else []
            rows = [
                {
                    "ensembl_gene": h["target"]["id"],
                    "identity": h["target"].get("perc_id"),
                    "similarity": h["target"].get("perc_pos"),
                    "type": h.get("type"),
                    "level": h.get("taxonomy_level"),
                }
                for h in homologies
            ]
            rows.sort(key=lambda r: -(r["identity"] or 0))
            self.cache.put(ensembl_gene, rows)
        rows = self._name(rows)
        ev = [
            database(
                "Ensembl Compara (human paralogues)",
                f"{len(rows)} human paralogues of {gene}; closest at {rows[0]['identity']:.0f}% identity"
                if rows
                else f"no human paralogues of {gene}",
                0.9,
                "computational",
                identifier=ensembl_gene,
                url=f"https://www.ensembl.org/Homo_sapiens/Gene/Compara_Paralog?g={ensembl_gene}",
            )
        ]
        return Answer(rows, ev, True, "", today())

    def _name(self, rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
        """Attach gene symbols, looked up in one batch and cached."""
        need = [r["ensembl_gene"] for r in rows if not r.get("symbol")]
        known = {g: self.symbols.get(g) for g in need}
        missing = [g for g, v in known.items() if v is None]
        if missing and self.net:
            body = json.dumps({"ids": missing[:250]}).encode()
            try:
                res = _fetch_json(f"{ENSEMBL_LOOKUP}?content-type=application/json", data=body, timeout=60)
            except UnavailableError:
                res = {}
            for g, rec in (res or {}).items():
                sym = (rec or {}).get("display_name")
                if sym:
                    known[g] = sym
                    self.symbols.put(g, sym)
        for r in rows:
            r["symbol"] = known.get(r["ensembl_gene"]) or r["ensembl_gene"]
        return rows


@dataclass
class EnsemblTranscriptProvider:
    """Coding sequences from Ensembl, keyed by the transcript VEP actually used.

    A frameshift's novel peptide stretch cannot be read off a VCF record: the
    indel has to be applied to the transcript and the result translated. That
    needs the same transcript the consequence was called on, which is why the
    transcript id is carried down from VEP rather than assumed.
    """

    net: bool = True
    cache: DiskCache = field(default_factory=lambda: DiskCache("cds"))

    def cds(self, transcript: str) -> Answer:
        if not transcript:
            return Answer.none("no transcript id on this variant's annotation")
        key = transcript.split(".")[0]
        seq = self.cache.get(key)
        if seq is None:
            if not self.net:
                return Answer.none(f"no cached coding sequence for {transcript} and network disabled")
            try:
                d = _fetch_json(f"{ENSEMBL_SEQUENCE}/{key}?type=cds;content-type=application/json")
            except UnavailableError as e:
                return Answer.none(f"Ensembl sequence unavailable for {transcript}: {e}")
            seq = (d or {}).get("seq")
            if not seq:
                return Answer.none(f"Ensembl returned no coding sequence for {transcript}")
            self.cache.put(key, seq)
        ev = [
            database(
                "Ensembl REST (coding sequence)",
                f"coding sequence of {key}, {len(seq)} nt",
                0.95,
                "human",
                identifier=key,
                url=f"https://www.ensembl.org/Homo_sapiens/Transcript/Summary?t={key}",
            )
        ]
        return Answer(seq, ev, True, "", today())


@dataclass
class NoNeoantigenPredictor:
    """The honest default: no peptide/HLA predictor is wired in.

    IEDB's and NetMHCpan's prediction services have licence or deployment
    terms that this project does not assume. Rather than invent affinities,
    the pipeline enumerates mutation-spanning peptides, records the patient's
    alleles, and reports binding and processing as unavailable. Implement this
    protocol to plug a real predictor in; nothing else has to change.
    """

    @property
    def name(self) -> str:
        return "none configured"

    def predict(self, peptides: list[str], alleles: list[str]) -> Answer:
        return Answer.none(
            "no peptide/HLA binding predictor is configured; GenomeOS does not estimate binding "
            "affinity without one"
        )


# --- the bundle -----------------------------------------------------------------------


@dataclass
class Providers:
    """Everything the pipeline is allowed to ask. Swap any one out freely."""

    protein: CompiledProteinProvider
    expression: HpaExpressionProvider
    structure: CompiledStructureProvider
    trafficking: AnnotationTraffickingProvider
    precedent: OpenTargetsProvider
    homology: EnsemblParalogueProvider
    neoantigen: Any
    patient_rna: PatientRnaProvider
    cohort: CBioPortalCohortProvider
    transcript: EnsemblTranscriptProvider

    @classmethod
    def default(
        cls,
        net: bool = True,
        knowledge_base: Any = None,
        patient_rna: PatientRnaProvider | None = None,
        cohort_study: str = "",
    ) -> Providers:
        protein = CompiledProteinProvider(net=net)
        return cls(
            protein=protein,
            expression=HpaExpressionProvider(net=net),
            structure=CompiledStructureProvider(protein),
            trafficking=AnnotationTraffickingProvider(protein, knowledge_base),
            precedent=OpenTargetsProvider(net=net),
            homology=EnsemblParalogueProvider(net=net),
            neoantigen=NoNeoantigenPredictor(),
            patient_rna=patient_rna or PatientRnaProvider(),
            cohort=CBioPortalCohortProvider(study=cohort_study, net=net),
            transcript=EnsemblTranscriptProvider(net=net),
        )

    def versions(self) -> dict[str, str]:
        """What answered this run, for the dataset's provenance block."""
        return {
            "protein_annotation": "GenomeOS federated protein compiler "
            "(UniProt, Ensembl, InterPro, PDB, AlphaFold, Reactome, STRING, HPA)",
            "normal_expression": HPA_LICENCE,
            "therapeutic_precedent": OT_LICENCE,
            "homology": "Ensembl REST comparative genomics (Apache 2.0 service, EMBL-EBI terms)",
            "transcript_sequence": "Ensembl REST sequence endpoint (Apache 2.0 service, EMBL-EBI terms)",
            "cohort_expression": (
                f"cBioPortal {self.cohort.study} (CC BY-SA 4.0 study data, ODbL portal)"
                if self.cohort.available
                else "no cohort selected"
            ),
            "neoantigen_predictor": getattr(self.neoantigen, "name", "unknown"),
            "patient_rna": self.patient_rna.path or "not supplied",
        }
