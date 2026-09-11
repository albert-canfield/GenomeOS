"""The peptide/HLA route: how an intracellular mutation becomes visible outside.

An intracellular mutation is not a dead end. Proteasomes chop the mutant
protein, TAP moves the fragments, and HLA class I carries some of them to the
surface, where a T-cell receptor or a TCR-mimic antibody can see them. The
target is then the mutant peptide *in complex with a particular HLA molecule*,
not the protein.

Four things are kept strictly apart, because conflating them is how
neoantigen pipelines overpromise:

    peptide exists        derived from the variant and the reference protein
    predicted binding     a predictor's output for one allele
    predicted processing  a different predictor's output
    observed presentation immunopeptidomics from real cells

Not every variant class yields a novel peptide. A stop-gained truncation
leaves every remaining residue identical to wild type, so it creates no
mutation-derived epitope; that is stated rather than papered over.
"""

from __future__ import annotations

import re
from typing import Any

from .evidence import Missing, database, derived
from .model import NeoantigenAssessment, PeptideCandidate

#: Class I peptides are 8-11 residues; 9-mers dominate.
PEPTIDE_LENGTHS = (8, 9, 10, 11)

#: What each variant class does to the protein sequence.
NOVEL_SEQUENCE: dict[str, tuple[bool | None, str]] = {
    "missense_variant": (True, "one residue is replaced, so every peptide spanning it is tumour-specific"),
    "protein_altering_variant": (True, "the protein sequence changes at the altered residue"),
    "inframe_deletion": (True, "residues are removed, creating a novel junction sequence"),
    "inframe_insertion": (True, "residues are inserted, creating a novel junction sequence"),
    "frameshift_variant": (
        True,
        "the reading frame shifts, producing a novel C-terminal stretch that is often the richest "
        "source of neoepitopes",
    ),
    "stop_gained": (
        False,
        "a premature stop truncates the protein but leaves every remaining residue identical to wild "
        "type, so the substitution itself creates no novel peptide",
    ),
    "start_lost": (
        None,
        "translation initiation is lost; any downstream product would need translation "
        "of the alternative start, which is not modelled here",
    ),
    "splice_acceptor_variant": (
        None,
        "a splice change can create a novel junction peptide, but the "
        "retained-intron or exon-skipped sequence is not reconstructed here",
    ),
    "splice_donor_variant": (
        None,
        "a splice change can create a novel junction peptide, but the "
        "retained-intron or exon-skipped sequence is not reconstructed here",
    ),
    "synonymous_variant": (False, "the protein sequence is unchanged"),
    "intron_variant": (False, "the protein sequence is unchanged"),
    "intergenic_variant": (False, "no protein product is involved"),
    "upstream_gene_variant": (False, "the protein sequence is unchanged"),
    "downstream_gene_variant": (False, "the protein sequence is unchanged"),
    "3_prime_UTR_variant": (False, "the protein sequence is unchanged"),
    "5_prime_UTR_variant": (False, "the protein sequence is unchanged"),
}

HLA_PATTERN = re.compile(r"^HLA-[A-DRQP]+[0-9]*\*\d{2}:\d{2}", re.IGNORECASE)
HLA_COMPACT = re.compile(r"^HLA-([A-DRQP]+[0-9]*)\*?(\d{2})(\d{2})$", re.IGNORECASE)
HLA_NO_STAR = re.compile(r"^HLA-([A-DRQP]+[0-9]*)(\d{2}:\d{2}.*)$", re.IGNORECASE)

REQUIRED_FOR_PRIORITISATION = (
    "patient HLA genotype (class I, four-digit)",
    "tumour RNA expression of the mutated gene",
    "a validated peptide/HLA binding and processing predictor",
    "immunopeptidomics from the tumour, to show the peptide is actually presented",
)


def normalise_alleles(alleles: list[str]) -> tuple[list[str], list[str]]:
    """Split supplied HLA strings into recognised alleles and unrecognised text."""
    good, bad = [], []
    for a in alleles:
        s = a.strip().upper().replace("HLA_", "HLA-")
        if not s:
            continue
        if not s.startswith("HLA-"):
            s = "HLA-" + s
        m = HLA_COMPACT.match(s) or HLA_NO_STAR.match(s)
        if m:
            g = m.groups()
            s = f"HLA-{g[0]}*{g[1]}:{g[2]}" if len(g) == 3 else f"HLA-{g[0]}*{g[1]}"
        (good if HLA_PATTERN.match(s) else bad).append(s)
    return good, bad


def peptides_around(
    sequence: str, residue: int, mutant_residue: str, lengths: tuple[int, ...] = PEPTIDE_LENGTHS
) -> list[PeptideCandidate]:
    """Every class-I-length window of the mutant protein that spans the change."""
    i = residue - 1
    mutant = sequence[:i] + mutant_residue + sequence[i + 1 :]
    out: list[PeptideCandidate] = []
    for n in lengths:
        for start in range(max(0, i - n + 1), min(i + 1, len(sequence) - n + 1)):
            end = start + n
            if end > len(mutant):
                continue
            out.append(
                PeptideCandidate(
                    sequence=mutant[start:end],
                    wild_type_sequence=sequence[start:end],
                    length=n,
                    start=start + 1,
                    end=end,
                    mutation_offset=i - start + 1,
                )
            )
    return out


def assess(
    gene: str,
    consequence: str,
    protein_change: str,
    residue: int | None,
    sequence: str | None,
    hla_alleles: list[str],
    predictor: Any,
) -> NeoantigenAssessment:
    """The peptide/HLA route for one variant, layered and honest about gaps."""
    a = NeoantigenAssessment()
    novel, why = NOVEL_SEQUENCE.get(consequence, (None, f"variant class '{consequence}' is not modelled"))
    a.novel_peptide_sequence = novel
    a.reason = why
    a.applicable = novel is not False

    if novel is False:
        a.evidence.append(
            derived(
                "GenomeOS neoantigen reasoning",
                f"{gene} {protein_change or consequence}: {why}",
                0.85,
            )
        )
        return a
    if novel is None:
        a.evidence.append(derived("GenomeOS neoantigen reasoning", f"{gene} {consequence}: {why}", 0.4))
        a.missing.append(
            Missing(
                "mutant protein sequence",
                why,
                (
                    "transcript-level translation of the altered allele (GenomeOS variant engine, needs the "
                    "reference genome and GENCODE for this chromosome)",
                ),
            )
        )

    if novel and consequence in ("missense_variant", "protein_altering_variant"):
        a.peptides, note = _missense_peptides(gene, protein_change, residue, sequence, a)
        if note:
            a.reason = note
    elif novel:
        a.missing.append(
            Missing(
                "mutant peptide sequences",
                f"{consequence} needs the altered transcript translated to recover the novel stretch; "
                "GenomeOS does not reconstruct it from the VCF record alone",
                ("reference genome and GENCODE annotation for this chromosome",),
            )
        )

    if a.peptides:
        a.wild_type_discrimination = (
            "a binder must tell the mutant peptide from the wild-type peptide, which differs at a single "
            "residue inside the same HLA groove; this is the hardest part of the design"
        )
        a.evidence.append(
            derived(
                "GenomeOS neoantigen reasoning",
                f"{len(a.peptides)} mutation-spanning peptides of length "
                f"{min(p.length for p in a.peptides)}-{max(p.length for p in a.peptides)} derived for "
                f"{gene} {protein_change}",
                0.7,
            )
        )

    alleles, unrecognised = normalise_alleles(hla_alleles)
    a.hla_alleles = alleles
    if unrecognised:
        a.missing.append(
            Missing(
                "HLA alleles",
                f"could not parse: {', '.join(unrecognised)}",
                ("four-digit typing, for example HLA-A*02:01",),
            )
        )
    if not alleles:
        a.missing.append(
            Missing(
                "patient HLA genotype",
                "without the patient's class I alleles no peptide can be prioritised; presentation is "
                "allele-specific",
                ("HLA typing from the tumour or normal sample",),
            )
        )
    else:
        a.evidence.append(
            database(
                "patient HLA genotype (supplied)",
                f"class I alleles available for prioritisation: {', '.join(alleles)}",
                0.9,
                "human",
            )
        )

    answer = predictor.predict([p.sequence for p in a.peptides], alleles) if a.peptides else None
    if answer is not None and answer.available:
        a.predicted_binding = answer.data
        a.evidence.extend(answer.evidence)
    else:
        reason = answer.reason if answer is not None else "no peptides to submit"
        a.predicted_binding = None
        a.missing.append(Missing("predicted peptide/HLA binding", reason, ("a wired binding predictor",)))
    a.predicted_processing = None
    a.missing.append(
        Missing(
            "predicted antigen processing",
            "proteasomal cleavage and TAP transport are not modelled; binding affinity alone "
            "overestimates presentation",
            ("a processing predictor",),
        )
    )
    a.observed_immunopeptidomics = None
    a.missing.append(
        Missing(
            "observed presentation",
            "no immunopeptidomics supplied; a predicted binder is not proof that the peptide is "
            "naturally presented",
            ("mass-spectrometry immunopeptidomics of the tumour",),
        )
    )
    a.clinical_evidence = None
    return a


def _missense_peptides(
    gene: str, protein_change: str, residue: int | None, sequence: str | None, a: NeoantigenAssessment
) -> tuple[list[PeptideCandidate], str]:
    if not sequence:
        a.missing.append(
            Missing(
                "reference protein sequence",
                f"no reviewed UniProt sequence available for {gene}",
                ("a UniProt entry for this gene",),
            )
        )
        return [], ""
    if residue is None or len(protein_change) < 3:
        a.missing.append(Missing("protein position", "the variant annotation carries no residue number", ()))
        return [], ""
    if residue > len(sequence):
        note = (
            f"residue {residue} lies beyond the {len(sequence)}-residue canonical sequence: the variant "
            f"was called on a longer {gene} isoform, so peptides are not derived from the canonical one"
        )
        a.missing.append(
            Missing(
                "isoform-matched protein sequence",
                note,
                ("the UniProt isoform matching the transcript the variant was called on",),
            )
        )
        return [], note
    wild, mutant = protein_change[0], protein_change[-1]
    if sequence[residue - 1] != wild:
        note = (
            f"canonical sequence carries {sequence[residue - 1]} at residue {residue}, not the "
            f"{wild} the annotation reports: isoform mismatch, peptides not derived"
        )
        a.missing.append(Missing("isoform-matched protein sequence", note, ()))
        return [], note
    peptides = peptides_around(sequence, residue, mutant)
    a.evidence.append(
        database(
            "UniProtKB/Swiss-Prot canonical sequence",
            f"{gene} residue {residue} is {wild} in the reference, replaced by {mutant} in the tumour",
            0.9,
            "human",
        )
    )
    return peptides, ""


def neoantigen_strength(a: NeoantigenAssessment | None) -> tuple[float | None, str]:
    """How much of the peptide/HLA route this variant actually supports.

    Deliberately capped low without a predictor and an HLA type: the score
    reflects how far the reasoning got, never a predicted immunogenicity.
    """
    if a is None or a.novel_peptide_sequence is False:
        return 0.0, (a.reason if a else "no neoantigen assessment")
    if not a.peptides:
        return None, "no mutant peptide could be derived; " + (a.reason or "sequence unavailable")
    score = 0.3
    parts = ["mutation-spanning peptides derived (0.30)"]
    if a.hla_alleles:
        score += 0.2
        parts.append("patient HLA genotype supplied (+0.20)")
    if a.predicted_binding:
        score += 0.2
        parts.append("binding prediction available (+0.20)")
    if a.observed_immunopeptidomics:
        score += 0.3
        parts.append("observed presentation (+0.30)")
    return round(score, 3), " ".join(parts) + "; ceiling without presentation evidence is 0.70"


def presentation_confidence(a: NeoantigenAssessment | None) -> tuple[float | None, str]:
    """Separate from strength: how sure are we the peptide reaches the surface."""
    if a is None or not a.applicable:
        return None, "peptide/HLA route does not apply to this variant"
    if a.observed_immunopeptidomics:
        return 0.9, "observed in immunopeptidomics"
    if a.predicted_binding and a.predicted_processing:
        return 0.4, "binding and processing predicted, not observed"
    if a.predicted_binding:
        return 0.25, "binding predicted only; processing not modelled"
    return 0.0, "no binding or processing evidence; presentation cannot be established"
