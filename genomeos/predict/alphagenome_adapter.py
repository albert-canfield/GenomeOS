"""AlphaGenome adapter (task 3.4): predicted regulatory effects as `predicted` rules.

The `alphagenome` client is an optional extra (`uv sync --extra predict`) and
needs ALPHAGENOME_API_KEY (free, non-commercial). Without a key the adapter
still works with an injected `scorer` callable, which is how the tests run.

Output is always BioIR: each predicted tissue effect becomes a Rule with
evidence kind `predicted`, the model name as source, and confidence derived
from the effect magnitude. Nothing here is ever labelled experimental.
"""

from __future__ import annotations

import importlib.util
import os
from collections.abc import Callable
from dataclasses import dataclass

from genomeos.ir import Action, Entity, Evidence, EvidenceKind, Module, Rule

MODEL_NAME = "AlphaGenome (google-deepmind/alphagenome 0.9)"
KEY_VAR = "ALPHAGENOME_API_KEY"
HOW_TO_ENABLE = (
    "install the client with `uv sync --extra predict` and put ALPHAGENOME_API_KEY=... in a git-ignored "
    ".env file at the project root (or export it); the key is free for non-commercial use from "
    "https://deepmind.google.com/science/alphagenome/account/terms"
)
LICENCE_NOTE = (
    "AlphaGenome is free for non-commercial use under Google DeepMind's terms; its predictions enter "
    "GenomeOS only as `predicted` evidence capped at confidence 0.7, and nothing in the core depends on it"
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
        "predicted target gene and tissue for UNKNOWN regulatory elements, tightening the CTCF-domain "
        "inference",
        "planned",
    ),
    (
        "c",
        "sequence grammar",
        "splice-site and chromatin signals the learned matrices cannot capture, entering as predicted "
        "evidence",
        "planned",
    ),
    (
        "d",
        "twin",
        "predicted expression differences between an individual's haplotypes and the reference",
        "planned",
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
    def confidence(self) -> float:
        """Magnitude-based confidence, capped at 0.7 because it is a prediction."""
        return min(0.7, 0.2 + abs(self.log2_fold_change) * 0.25)


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
            self._client = dna_client.create(self.api_key)
        client = self._client

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
                for gi, gene in enumerate(genes):
                    for ti, tissue in enumerate(tissues):
                        val = float(adata.X[gi, ti])
                        if abs(val) > self.last_scan["max_abs_log2fc"]:
                            self.last_scan["max_abs_log2fc"] = abs(val)
                        if abs(val) > threshold:
                            out.append((str(gene), str(tissue), val))
            return out

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
                        EvidenceKind.PREDICTED, MODEL_NAME, note=f"log2FC={e.log2_fold_change:+.3f}"
                    ),
                    confidence=e.confidence,
                )
            )
        return m
