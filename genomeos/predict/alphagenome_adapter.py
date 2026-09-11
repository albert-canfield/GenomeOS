"""AlphaGenome adapter (task 3.4): predicted regulatory effects as `predicted` rules.

The `alphagenome` client is an optional extra (`uv sync --extra predict`) and
needs ALPHAGENOME_API_KEY (free, non-commercial). Without a key the adapter
still works with an injected `scorer` callable, which is how the tests run.

Output is always BioIR: each predicted tissue effect becomes a Rule with
evidence kind `predicted`, the model name as source, and confidence derived
from the effect magnitude. Nothing here is ever labelled experimental.
"""

from __future__ import annotations

import os
from collections.abc import Callable
from dataclasses import dataclass

from genomeos.ir import Action, Entity, Evidence, EvidenceKind, Module, Rule

MODEL_NAME = "AlphaGenome (google-deepmind/alphagenome 0.9)"

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


class AlphaGenomeAdapter:
    def __init__(self, scorer: Scorer | None = None, api_key: str | None = None) -> None:
        self._scorer = scorer
        self._client = None
        self.api_key = api_key or os.environ.get("ALPHAGENOME_API_KEY") or _dotenv_key("ALPHAGENOME_API_KEY")

    @property
    def available(self) -> bool:
        return self._scorer is not None or bool(self.api_key)

    def _live_scorer(self) -> Scorer:
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
                tissues = list(adata.var.get("biosample_name", adata.var.index))
                for gi, gene in enumerate(genes):
                    for ti, tissue in enumerate(tissues):
                        val = float(adata.X[gi, ti])
                        if abs(val) > 0.05:
                            out.append((str(gene), str(tissue), val))
            return out

        return score

    def predict(self, chrom: str, pos_1based: int, ref: str, alt: str) -> list[PredictedEffect]:
        scorer = self._scorer or self._live_scorer()
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
