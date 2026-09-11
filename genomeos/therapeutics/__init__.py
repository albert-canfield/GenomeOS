"""Cancer targeting and therapeutic mechanism reasoning.

From one tumour's molecular alterations: what distinguishes those cells from
healthy ones, which of those differences a therapy could physically reach,
and which therapeutic mechanism the biology supports for each.

    from genomeos.therapeutics import analyse_vcf, text_report, write_outputs

    analysis = analyse_vcf("data/demo/cancer_tumour.vcf", hla=["HLA-A*02:01"])
    print(text_report(analysis))
    write_outputs(analysis, "out/")

Research hypotheses, rankings, evidence and mechanistic reasoning only. Not
clinical recommendations, and no laboratory-ready constructs of any kind.
"""

from .design import (
    SCHEMA,
    SCHEMA_VERSION,
    dataset,
    design_readiness,
    evidence_graph,
    negative_set,
    negative_targets,
    positive_set,
    specification,
)
from .evidence import Evidence, EvidenceLedger, Missing
from .mechanisms import MECHANISMS, MechanismSpec
from .model import (
    Localisation,
    MechanismFit,
    NeoantigenAssessment,
    TargetLogic,
    TherapeuticTargetCandidate,
    Trafficking,
    VariantOrigin,
)
from .pipeline import (
    DISCLAIMER,
    PatientProfile,
    analyse,
    analyse_vcf,
    build_candidate,
    rank_variants,
    write_outputs,
)
from .providers import Providers
from .report import machine_report, target_specification_text, text_report

__all__ = [
    "DISCLAIMER",
    "Evidence",
    "EvidenceLedger",
    "Localisation",
    "MECHANISMS",
    "MechanismFit",
    "MechanismSpec",
    "Missing",
    "NeoantigenAssessment",
    "PatientProfile",
    "Providers",
    "SCHEMA",
    "SCHEMA_VERSION",
    "TargetLogic",
    "TherapeuticTargetCandidate",
    "Trafficking",
    "VariantOrigin",
    "analyse",
    "analyse_vcf",
    "build_candidate",
    "dataset",
    "design_readiness",
    "evidence_graph",
    "machine_report",
    "negative_set",
    "negative_targets",
    "positive_set",
    "rank_variants",
    "specification",
    "target_specification_text",
    "text_report",
    "write_outputs",
]
