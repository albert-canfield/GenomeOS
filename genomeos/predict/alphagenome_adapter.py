"""AlphaGenome adapter (task 3.4): predicted regulatory effects as `predicted` rules.

The `alphagenome` client is an optional extra (`uv sync --extra predict`) and
needs ALPHAGENOME_API_KEY (free, non-commercial). Without a key the adapter
still works with an injected `scorer` callable, which is how the tests run.

Output is always BioIR: each predicted tissue effect becomes a Rule with
evidence kind `predicted`, the model name as source, and no confidence: since
review R4 (R4f, 2026-09-28) an effect carries its log2 fold change in its unit
and a `genomeos.certainty.Certainty` record (model score = |log2 fold change|,
named; probability None with its reason). The size of a predicted effect is not
how sure anyone is, so nothing here converts it into a certainty. Nothing here
is ever labelled experimental.
"""

from __future__ import annotations

import datetime
import hashlib
import importlib.metadata
import importlib.util
import os
from collections.abc import Callable
from dataclasses import dataclass

from genomeos.certainty import Certainty
from genomeos.ir import Action, Entity, Evidence, EvidenceKind, Module, Rule
from genomeos.ir.model import UNSTATED

MODEL_NAME = "AlphaGenome (google-deepmind/alphagenome 0.9)"
KEY_VAR = "ALPHAGENOME_API_KEY"
HOW_TO_ENABLE = (
    "install the client with `uv sync --extra predict` and put ALPHAGENOME_API_KEY=... in a git-ignored "
    ".env file at the project root (or export it); the key is free for non-commercial use from "
    "https://deepmind.google.com/science/alphagenome/account/terms"
)
LICENCE_NOTE = (
    "AlphaGenome is free for non-commercial use under Google DeepMind's terms; its predictions enter "
    "GenomeOS only as `predicted` evidence, each effect with its unit and no probability, and nothing in the"
    " core depends on it"
)
EFFECT_UNIT = (
    "log2 fold change of the gene's predicted RNA-seq expression, alternate allele over reference"
    " (AlphaGenome RNA_SEQ variant scorer, one track)"
)
MODEL_SCORE_NAME = "|log2 fold change| on that track: the magnitude effects are sorted by, not a probability"
NO_PROBABILITY = (
    "no calibration record: no predicted variant effect has been scored against a measured outcome"
    " population by a stated method, so no probability that this variant moves this gene is quoted"
)

# What the key enables. Every feature is always listed; without the key it is loaded but disabled.
FEATURES: tuple[tuple[str, str, str, str], ...] = (
    (
        "a",
        "variant effect",
        "predicted expression change per tissue for any variant (`genomeos predict`, `/api/predict`)",
        "built",
    ),
    (
        "b",
        "regulatory blocks",
        "predicted target gene and tissue for enhancer-like elements by deleting each one in its window "
        "(`genomeos predict --element`, job enhancer_targets_chr21), held against the CTCF-domain inference",
        "built",
    ),
    (
        "c",
        "sequence grammar",
        "splice sites at 1 bp as the donor/acceptor candidates of the segment parser "
        "(`genomeos segments --predicted-sites`), entering as predicted evidence",
        "built",
    ),
    (
        "d",
        "twin",
        "a person's variants inside the elements that reach a gene, each scored for its predicted effect "
        "on the gene and summed per haplotype (`genomeos individual predict`)",
        "built",
    ),
)

Scorer = Callable[[str, int, str, str], list[tuple[str, str, float]]]
"""(chrom, pos_1based, ref, alt) -> [(gene_symbol, tissue, log2_fold_change), ...]"""


@dataclass(frozen=True, slots=True)
class PredictedEffect:
    gene: str
    tissue: str
    log2_fold_change: float

    @property
    def direction(self) -> Action:
        return Action.ACTIVATE if self.log2_fold_change > 0 else Action.INHIBIT

    @property
    def certainty(self) -> Certainty:
        """What the effect rests on (R4f): the effect in its unit, the score named, no probability.
        Only the effect and the score depend on the magnitude; nothing converts either into a
        certainty. It replaces `confidence = min(0.7, 0.2 + 0.25 |log2FC|)`."""
        return Certainty(
            evidence_category=f"predicted: {MODEL_NAME}, one model run",
            effect_estimate=self.log2_fold_change,
            effect_unit=EFFECT_UNIT,
            uncertainty_note="one deterministic model run; no spread computed",
            model_score=abs(self.log2_fold_change),
            model_score_name=MODEL_SCORE_NAME,
            probability_unavailable=NO_PROBABILITY,
        )


API_SERVICE = {
    "service": "google.gdm.gdmscience.alphagenome.v1main.DnaModelService",
    "address": "dns:///gdmscience.googleapis.com:443",
}
UNREQUESTED = (
    "unrequested: dna_client.create was given no model_version, so the server chose; alphagenome 0.9.0 "
    "documents ALL_FOLDS as its default, and no response says which model answered"
)


#: The model version every AlphaGenome client in this project asks for (coordinator decision, ROADMAP
#: section 5 item 11, R9 model-dependency row, 2026-09-28). The name of a member of the installed
#: client's `alphagenome.models.dna_model.ModelVersion` enum (0.9.0: ALL_FOLDS, FOLD_0..FOLD_3), passed
#: as `dna_client.create(api_key, model_version=...)`. ALL_FOLDS is the distilled all-folds model the
#: 0.9.0 client names as its default; asking for it by name makes each request carry it instead of an
#: empty field. Answers made before 2026-09-28 asked for nothing and are labelled "unrequested", since
#: the server's choice then is not guaranteed to be the model that answers a named request.
ALPHAGENOME_MODEL_VERSION = "ALL_FOLDS"


def create_client(api_key: str | None, **kwargs: object) -> object:
    """A live AlphaGenome client that requests ALPHAGENOME_MODEL_VERSION; every place that makes a
    client goes through here. `kwargs` pass on to `dna_client.create` (timeout, address)."""
    from alphagenome.models import dna_client  # type: ignore[import-not-found]

    version = dna_client.ModelVersion[ALPHAGENOME_MODEL_VERSION]
    return dna_client.create(api_key, model_version=version, **kwargs)


def run_metadata(client: object | None = None, date: str | None = None) -> dict:
    """What a live request is made with, kept beside each cached answer (review R9 follow-up): the client
    package and its installed version, the model version the client asks for (None when it asks for
    none, which is how every run before 2026-09-28 asked), the service and the date. Reads the installed
    package and the client object only; makes no request."""
    try:
        version = importlib.metadata.version("alphagenome")
    except importlib.metadata.PackageNotFoundError:
        version = None
    requested = getattr(client, "_model_version", None)  # DnaClient keeps the requested ModelVersion's name
    if client is None:  # no client: the version create_client asks for (a client made with none says None)
        requested = ALPHAGENOME_MODEL_VERSION
    return {
        "client": "alphagenome",
        "client_version": version,
        "model_version": requested,
        "model_version_note": None if requested else UNREQUESTED,
        "api": dict(API_SERVICE),
        "scorer": "variant_scorers.RECOMMENDED_VARIANT_SCORERS['RNA_SEQ']",
        "date": date or datetime.datetime.now(datetime.UTC).date().isoformat(),
    }


def _dotenv_key(name: str, path: str = ".env") -> str | None:
    """Read NAME=value from a local .env file (git-ignored) when the variable is not exported."""
    try:
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if line.startswith(f"{name}="):
                    return line.split("=", 1)[1].strip().strip("'\"") or None
    except OSError:
        return None
    return None


def status(api_key: str | None = None, dotenv: str = ".env") -> dict:
    """Whether the optional AlphaGenome features are enabled, and if not, why and how to enable them."""
    package = importlib.util.find_spec("alphagenome") is not None
    key = bool(api_key or os.environ.get(KEY_VAR) or _dotenv_key(KEY_VAR, dotenv))
    missing = []
    if not package:
        missing.append("the `alphagenome` package is not installed")
    if not key:
        missing.append(f"{KEY_VAR} is not set")
    return {
        "name": "AlphaGenome",
        "enabled": package and key,
        "package": package,
        "key": key,
        "reason": "; ".join(missing),
        "how": HOW_TO_ENABLE,
        "licence": LICENCE_NOTE,
        "model": MODEL_NAME,
        "features": [{"id": i, "name": n, "what": w, "state": st} for i, n, w, st in FEATURES],
    }


class AlphaGenomeAdapter:
    def __init__(self, scorer: Scorer | None = None, api_key: str | None = None) -> None:
        self._scorer = scorer
        self._client = None
        self.last_scan: dict = {}  # genes, tracks and the largest |log2FC| seen by the last live call
        self.api_key = api_key or os.environ.get(KEY_VAR) or _dotenv_key(KEY_VAR)

    @property
    def available(self) -> bool:
        return self._scorer is not None or bool(self.api_key)

    def _live_scorer(self, threshold: float = 0.05) -> Scorer:
        from alphagenome.data import genome  # type: ignore[import-not-found]
        from alphagenome.models import dna_client, variant_scorers  # type: ignore[import-not-found]

        if self._client is None:
            self._client = create_client(self.api_key)
        client = self._client
        model = run_metadata(client)

        def score(chrom: str, pos: int, ref: str, alt: str) -> list[tuple[str, str, float]]:
            variant = genome.Variant(chromosome=chrom, position=pos, reference_bases=ref, alternate_bases=alt)
            interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
            scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS["RNA_SEQ"]
            scores = client.score_variant(interval=interval, variant=variant, variant_scorers=[scorer])
            out: list[tuple[str, str, float]] = []
            for adata in scores:
                genes = list(adata.obs.get("gene_name", []))
                names = list(adata.var.get("biosample_name", adata.var.index))
                gtex = list(adata.var.get("gtex_tissue", [None] * len(names)))
                tissues = []
                for g, n in zip(gtex, names, strict=False):
                    tissues.append(str(g) if g and str(g) not in ("nan", "") else str(n))
                self.last_scan = {"genes": len(genes), "tracks": len(tissues), "max_abs_log2fc": 0.0}
                # the track table this answer was read on, as a checksum: no request, the response's own
                model["tracks"] = len(tissues)
                model["tracks_sha256"] = hashlib.sha256(
                    "\n".join(
                        f"{i}\t{n}\t{t}" for i, n, t in zip(adata.var.index, names, tissues, strict=False)
                    ).encode()
                ).hexdigest()
                for gi, gene in enumerate(genes):
                    for ti, tissue in enumerate(tissues):
                        val = float(adata.X[gi, ti])
                        if abs(val) > self.last_scan["max_abs_log2fc"]:
                            self.last_scan["max_abs_log2fc"] = abs(val)
                        if abs(val) > threshold:
                            out.append((str(gene), str(tissue), val))
            return out

        score.model = model  # type: ignore[attr-defined]  # read by enhancer_target.score_element
        return score

    def predict(
        self, chrom: str, pos_1based: int, ref: str, alt: str, threshold: float = 0.05
    ) -> list[PredictedEffect]:
        """Effects with |log2 fold change| above the threshold; `last_scan` says what was scanned."""
        scorer = self._scorer or self._live_scorer(threshold)
        return [PredictedEffect(g, t, v) for g, t, v in scorer(chrom, pos_1based, ref, alt)]

    def to_module(
        self, chrom: str, pos_1based: int, ref: str, alt: str, effects: list[PredictedEffect]
    ) -> Module:
        vid = f"{chrom}:{pos_1based}{ref}>{alt}"
        m = Module(name=f"predicted.{chrom}_{pos_1based}_{ref}_{alt}")
        m.add(
            Entity(
                id=vid,
                kind="variant",
                attrs={"chrom": chrom, "pos": pos_1based, "ref": ref, "alt": alt},
                evidence=Evidence(EvidenceKind.CURATED, "input"),
                confidence=1.0,
            )
        )
        for e in effects:
            if e.gene not in m.entities:
                m.add(
                    Entity(
                        id=e.gene,
                        kind="gene",
                        evidence=Evidence(EvidenceKind.CURATED, "GENCODE"),
                        confidence=0.9,
                    )
                )
            m.rules.append(
                Rule(
                    id=f"{vid} {e.direction.value} {e.gene} in {e.tissue}",
                    source=vid,
                    action=e.direction,
                    target=e.gene,
                    strength=min(1.0, abs(e.log2_fold_change)),
                    when={"tissue": e.tissue},
                    evidence=Evidence(
                        EvidenceKind.PREDICTED,
                        MODEL_NAME,
                        note=f"log2FC={e.log2_fold_change:+.3f}, probability unavailable",
                    ),
                    confidence=UNSTATED,  # R4f: an effect's size is not a confidence; none is stated
                )
            )
        return m
