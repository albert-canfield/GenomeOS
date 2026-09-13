"""Peptide/HLA binding predictors as optional features.

GenomeOS never invents an affinity. A predictor is wired in only if the user
has installed it, and each one is always *listed* with what it would bring and
how to enable it, so a disabled predictor is visible rather than silent (the
optional-feature rule, docs/DECISIONS.md D34).

Two predictors are supported, deliberately different in licence:

    MHCflurry 2   Apache 2.0, `pip install mhcflurry` plus a model download.
                  Usable by anyone, including commercially. The default.
    NetMHCpan 4.1 academic licence from DTU; a local binary the user installs
                  themselves (https://services.healthtech.dtu.dk). It is the
                  field's reference predictor, so a laboratory that already
                  holds a licence should be able to use it; GenomeOS neither
                  ships nor downloads it, and calls the binary it finds.

Both answer the `NeoantigenProvider` protocol: `predict(peptides, alleles) ->
Answer`. The output is always `computational` evidence with the predictor's
name and version; presentation still needs processing and immunopeptidomics,
which neither predictor provides.
"""

from __future__ import annotations

import importlib.util
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .evidence import prediction
from .providers import Answer

#: A peptide is called a binder below this percentile rank of the allele's
#: score distribution. 2.0 is the convention both tools document for class I;
#: 0.5 is the stricter "strong binder" threshold.
WEAK_RANK = 2.0
STRONG_RANK = 0.5

NETMHCPAN_ENV = "NETMHCPAN"
NETMHCPAN_LICENCE = (
    "NetMHCpan 4.1 is licensed by DTU Health Tech for academic use; other users must obtain a licence "
    "from health-software@dtu.dk. GenomeOS does not ship, download or redistribute it and only calls a "
    "binary the user has installed"
)
MHCFLURRY_LICENCE = "MHCflurry 2 is Apache 2.0 (openvax/mhcflurry); models are downloaded by the user"


@dataclass(frozen=True, slots=True)
class Binding:
    """One peptide against one allele, as the predictor reported it."""

    peptide: str
    allele: str
    percentile_rank: float
    affinity_nm: float | None = None

    @property
    def binder(self) -> bool:
        return self.percentile_rank <= WEAK_RANK

    @property
    def strong(self) -> bool:
        return self.percentile_rank <= STRONG_RANK

    def to_dict(self) -> dict[str, Any]:
        d: dict[str, Any] = {
            "peptide": self.peptide,
            "allele": self.allele,
            "percentile_rank": round(self.percentile_rank, 3),
            "binder": self.binder,
            "strong_binder": self.strong,
        }
        if self.affinity_nm is not None:
            d["affinity_nm"] = round(self.affinity_nm, 1)
        return d


def _summary(
    name: str, bindings: list[Binding], alleles: list[str], processing: bool = False
) -> dict[str, Any]:
    best = min(bindings, key=lambda b: b.percentile_rank) if bindings else None
    return {
        "predictor": name,
        "alleles": alleles,
        "peptides": len({b.peptide for b in bindings}),
        "binders": sorted((b.to_dict() for b in bindings if b.binder), key=lambda d: d["percentile_rank"]),
        "strong_binders": sum(1 for b in bindings if b.strong),
        "best": best.to_dict() if best else None,
        "threshold": {"binder_percentile_rank": WEAK_RANK, "strong_percentile_rank": STRONG_RANK},
        "processing_modelled": processing,
        "caveat": (
            "a predicted binder is not a presented peptide; antigen processing is "
            + ("included in this score but not observed" if processing else "not modelled here")
            + ", and only immunopeptidomics shows what the tumour actually displays"
        ),
    }


class MhcflurryPredictor:
    """MHCflurry 2 (Apache 2.0). Class I, peptide lengths 8-15."""

    licence = MHCFLURRY_LICENCE
    install = (
        "uv sync --extra hla (installs mhcflurry), then "
        "`mhcflurry-downloads fetch models_class1_presentation`"
    )

    def __init__(self) -> None:
        self._predictor: Any = None
        #: MHCflurry 2's presentation score combines binding with its own antigen-processing
        #: model, so processing is predicted here in a way NetMHCpan's rank is not.
        self._presentation = False

    @property
    def name(self) -> str:
        return "MHCflurry 2"

    @staticmethod
    def installed() -> bool:
        return importlib.util.find_spec("mhcflurry") is not None

    @property
    def available(self) -> bool:
        return self.installed()

    def _load(self) -> Any:
        if self._predictor is None:
            from mhcflurry import Class1PresentationPredictor  # type: ignore[import-not-found]

            self._predictor = Class1PresentationPredictor.load()
        return self._predictor

    def predict(self, peptides: list[str], alleles: list[str]) -> Answer:
        if not self.installed():
            return Answer.none(f"MHCflurry is not installed; {self.install}")
        if not peptides or not alleles:
            return Answer.none("no peptides or no HLA alleles supplied")
        try:
            df = self._load().predict(peptides=list(peptides), alleles=list(alleles), verbose=0)
        except Exception as e:  # noqa: BLE001 - a missing model download must not crash an analysis
            hint = self.install
            if "Missing MHCflurry downloadable file" in str(e):
                hint = (
                    "the models are not downloaded; run `genomeos therapeutic --install-hla-models` "
                    "(MHCflurry's own `mhcflurry-downloads fetch` cannot run on Python 3.13+, which "
                    "removed the `pipes` module)"
                )
            return Answer.none(f"MHCflurry could not run ({str(e)[:100]}); {hint}")
        bindings = []
        for row in df.to_dict("records"):
            rank = row.get("presentation_percentile")
            if rank is None:
                rank = row.get("affinity_percentile")
            else:
                self._presentation = True
            if rank is None:
                continue
            bindings.append(
                Binding(
                    peptide=str(row["peptide"]),
                    allele=str(row.get("best_allele", alleles[0])),
                    percentile_rank=float(rank),
                    affinity_nm=float(row["affinity"]) if row.get("affinity") is not None else None,
                )
            )
        data = _summary(self.name, bindings, list(alleles), processing=self._presentation)
        n = len(data["binders"])
        what = "presentation (binding plus antigen processing)" if self._presentation else "binding"
        return Answer(
            data,
            [
                prediction(
                    self.name,
                    f"{n} of {len(peptides)} mutation-spanning peptides predicted for {what} on "
                    f"{', '.join(alleles)} at percentile rank <= {WEAK_RANK}",
                    0.5,
                )
            ],
        )


class NetMHCpanPredictor:
    """NetMHCpan 4.1 through a locally installed binary (academic licence)."""

    licence = NETMHCPAN_LICENCE
    install = (
        "request the academic package at https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/, "
        f"install it, and put the executable on PATH or set {NETMHCPAN_ENV}=/path/to/netMHCpan"
    )

    def __init__(self, binary: str | None = None) -> None:
        self.binary = binary or self.find()

    @property
    def name(self) -> str:
        return "NetMHCpan 4.1"

    @staticmethod
    def find() -> str | None:
        return os.environ.get(NETMHCPAN_ENV) or shutil.which("netMHCpan")

    @property
    def available(self) -> bool:
        return bool(self.binary and Path(self.binary).exists())

    def predict(self, peptides: list[str], alleles: list[str]) -> Answer:
        if not self.available:
            return Answer.none(f"NetMHCpan is not installed; {self.install}")
        if not peptides or not alleles:
            return Answer.none("no peptides or no HLA alleles supplied")
        bindings: list[Binding] = []
        with tempfile.TemporaryDirectory() as tmp:
            pep = Path(tmp) / "peptides.txt"
            pep.write_text("\n".join(peptides) + "\n")
            cmd = [str(self.binary), "-p", str(pep), "-a", ",".join(a.replace("*", "") for a in alleles)]
            try:
                out = subprocess.run(  # noqa: S603 - a binary the user installed, fixed arguments
                    cmd, capture_output=True, text=True, timeout=600, check=False
                )
            except (OSError, subprocess.SubprocessError) as e:
                return Answer.none(f"NetMHCpan could not run ({str(e)[:120]})")
            if out.returncode != 0:
                return Answer.none(f"NetMHCpan exited {out.returncode}: {out.stderr.strip()[:160]}")
            bindings = parse_netmhcpan(out.stdout)
        data = _summary(self.name, bindings, list(alleles))
        n = len(data["binders"])
        return Answer(
            data,
            [
                prediction(
                    self.name,
                    f"{n} of {len(peptides)} mutation-spanning peptides predicted to bind "
                    f"{', '.join(alleles)} at percentile rank <= {WEAK_RANK}",
                    0.55,
                )
            ],
        )


def parse_netmhcpan(text: str) -> list[Binding]:
    """Read NetMHCpan's fixed-column output. Columns differ between versions, so it reads by header."""
    bindings: list[Binding] = []
    cols: list[str] | None = None
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith(("#", "-", "Distance", "Protein")):
            continue
        if line.startswith("Pos") and "Peptide" in line:
            cols = line.split()
            continue
        if cols is None or line.startswith("Pos"):
            continue
        f = line.split()
        if len(f) < len(cols) - 1:
            continue
        row = dict(zip(cols, f, strict=False))
        pep, allele = row.get("Peptide"), row.get("MHC") or row.get("HLA") or row.get("Allele")
        rank = row.get("%Rank_EL") or row.get("%Rank") or row.get("Rank")
        if not (pep and allele and rank):
            continue
        try:
            aff = row.get("Aff(nM)") or row.get("Affinity")
            bindings.append(
                Binding(pep, allele, float(rank), float(aff) if aff not in (None, "NA") else None)
            )
        except ValueError:
            continue
    return bindings


def fetch_mhcflurry_models(name: str = "models_class1_presentation", log: Any = None) -> str:
    """Install MHCflurry's models from the URL its own manifest names.

    MHCflurry ships `mhcflurry-downloads fetch`, but that command imports the
    `pipes` module, which Python 3.13 removed, so it cannot run on a modern
    interpreter. This does the same work: read the release URL out of the
    package's own `downloads.yml`, download it and unpack it where MHCflurry
    looks. Nothing is redistributed; the file comes from the openvax release.
    """
    import re
    import tarfile
    import urllib.request

    import mhcflurry  # type: ignore[import-not-found]
    from mhcflurry.downloads import get_downloads_dir  # type: ignore[import-not-found]

    manifest = (Path(mhcflurry.__file__).parent / "downloads.yml").read_text()
    block = manifest.split(f"name: {name}", 1)
    if len(block) != 2:
        raise RuntimeError(f"{name} is not in MHCflurry's downloads.yml")
    m = re.search(r"url:\s*(\S+)", block[1])
    if not m:
        raise RuntimeError(f"no url for {name} in MHCflurry's downloads.yml")
    url = m.group(1)
    target = Path(get_downloads_dir()) / name
    if target.exists():
        return str(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    if log:
        print(f"mhcflurry: fetching {name} from {url}", file=log, flush=True)
    with tempfile.TemporaryDirectory() as tmp:
        archive = Path(tmp) / "models.tar.bz2"
        with urllib.request.urlopen(url) as r, open(archive, "wb") as fh:  # noqa: S310 - https URL from the package
            shutil.copyfileobj(r, fh)
        with tarfile.open(archive, "r:bz2") as tf:
            tf.extractall(target, filter="data")
    if log:
        print(f"mhcflurry: {name} installed at {target}", file=log, flush=True)
    return str(target)


PREDICTORS: tuple[type[MhcflurryPredictor] | type[NetMHCpanPredictor], ...] = (
    MhcflurryPredictor,
    NetMHCpanPredictor,
)


def status() -> dict[str, Any]:
    """What each predictor would bring and whether it is installed. Always listed, never assumed."""
    rows = []
    for cls in PREDICTORS:
        p = cls()
        rows.append(
            {
                "name": p.name,
                "available": p.available,
                "licence": p.licence,
                "how": p.install,
                "what": (
                    "percentile-rank binding of every mutation-spanning peptide against the patient's "
                    "class I alleles; raises neoantigen strength and presentation confidence"
                ),
            }
        )
    chosen = choose()
    return {
        "name": "peptide/HLA binding",
        "enabled": getattr(chosen, "available", False),
        "using": getattr(chosen, "name", "none configured"),
        "predictors": rows,
        "caveat": (
            "a predicted binder is not a presented peptide; processing and immunopeptidomics stay "
            "unestablished whichever predictor is used"
        ),
    }


def choose(preferred: str = "") -> Any:
    """The predictor to use: the requested one, else the first installed, else the honest default."""
    from .providers import NoNeoantigenPredictor

    want = preferred.strip().lower()
    if want in ("none", "off"):
        return NoNeoantigenPredictor()
    for cls in PREDICTORS:
        p = cls()
        if want and want not in p.name.lower().replace(" ", ""):
            continue
        if p.available:
            return p
    return NoNeoantigenPredictor()
