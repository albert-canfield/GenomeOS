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
    hgvsc: str = "",
    transcript: str = "",
    transcripts: Any = None,
) -> NeoantigenAssessment:
    """The peptide/HLA route for one variant, layered and honest about gaps.

    Where the transcript the consequence was called on can be fetched, the
    altered protein is reconstructed from its coding sequence and compared
    with the reference. That replaces a rule about the variant class with a
    measurement of what the protein actually becomes, and it is the only way
    to recover a frameshift's novel C-terminal stretch.
    """
    a = NeoantigenAssessment()
    novel, why = NOVEL_SEQUENCE.get(consequence, (None, f"variant class '{consequence}' is not modelled"))
    a.novel_peptide_sequence = novel
    a.reason = why
    a.applicable = novel is not False

    # the patient's alleles are recorded whatever the variant turns out to be:
    # "not supplied" must mean not supplied, not "not reached"
    alleles, unrecognised = normalise_alleles(hla_alleles)
    a.hla_alleles = alleles

    rebuilt = _reconstruct(gene, consequence, hgvsc, transcript, transcripts, a)
    if rebuilt is not None:
        novel = a.novel_peptide_sequence
        why = a.reason

    if novel is False and not a.peptides:
        a.evidence.append(
            derived(
                "GenomeOS neoantigen reasoning",
                f"{gene} {protein_change or consequence}: {why}",
                0.85,
            )
        )
        return a
    if novel is None and rebuilt is None:
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

    if a.peptides:
        pass  # already recovered from the transcript
    elif novel and consequence in ("missense_variant", "protein_altering_variant"):
        a.peptides, note = _missense_peptides(gene, protein_change, residue, sequence, a)
        if note:
            a.reason = note
    elif novel:
        a.missing.append(
            Missing(
                "mutant peptide sequences",
                f"{consequence} needs the altered transcript translated to recover the novel stretch, "
                "and the transcript's coding sequence was not available",
                ("the Ensembl transcript this consequence was called on",),
            )
        )

    if a.peptides:
        frameshift = any(p.source == "frameshift_derived" for p in a.peptides)
        a.wild_type_discrimination = (
            "the altered reading frame produces residues that exist in no healthy protein, so there is "
            "no wild-type counterpart to discriminate against; this is the easiest discrimination "
            "problem in the whole pipeline, which is why frameshifts are the richest neoantigen source"
            if frameshift and any(not p.wild_type_sequence for p in a.peptides)
            else "a binder must tell the mutant peptide from the wild-type peptide, which differs at a "
            "single residue inside the same HLA groove; this is the hardest part of the design"
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
    if any(p.source == "frameshift_derived" and not p.wild_type_sequence for p in a.peptides):
        score += 0.1
        parts.append(
            "the altered frame produces residues with no wild-type counterpart, so a binder has "
            "nothing self to be confused with (+0.10)"
        )
    if a.hla_alleles:
        score += 0.2
        parts.append("patient HLA genotype supplied (+0.20)")
    if a.predicted_binding:
        score += 0.2
        parts.append("binding prediction available (+0.20)")
    if a.observed_immunopeptidomics:
        score += 0.3
        parts.append("observed presentation (+0.30)")
    return round(score, 3), " ".join(parts) + "; ceiling without presentation evidence is 0.60"


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


# --- reconstructing an altered protein from an indel ----------------------------------

#: HGVS coding-sequence changes this module can apply. Anything else, including
#: changes outside the CDS, is left to the caller as unreconstructed.
HGVS_C = re.compile(
    r"^c\.(?P<start>-?\*?\d+)(?:_(?P<end>-?\*?\d+))?"
    r"(?:(?P<op>del|dup|ins|delins)(?P<seq>[ACGT]*)|(?P<ref>[ACGT])>(?P<alt>[ACGT]))",
    re.IGNORECASE,
)


def apply_hgvs_c(cds: str, hgvsc: str) -> tuple[str, str] | None:
    """Apply one HGVS coding change to a coding sequence.

    Returns the altered sequence and a description, or None when the change is
    outside the coding sequence or in a form this module does not reconstruct.
    Positions with - or * prefixes are UTR coordinates and are declined rather
    than guessed at.
    """
    m = HGVS_C.match((hgvsc or "").strip())
    if not m:
        return None
    raw_start, raw_end = m.group("start"), m.group("end")
    if any(x and ("-" in x or "*" in x or "+" in x) for x in (raw_start, raw_end)):
        return None
    start = int(raw_start)
    end = int(raw_end) if raw_end else start
    if start < 1 or end < start or end > len(cds):
        return None
    i, j = start - 1, end  # 0-based half-open over the CDS
    op = (m.group("op") or "").lower()
    seq = (m.group("seq") or "").upper()
    if m.group("ref"):
        ref, alt = m.group("ref").upper(), m.group("alt").upper()
        if cds[i : i + 1].upper() != ref:
            return None
        return cds[:i] + alt + cds[i + 1 :], f"{ref}>{alt} at c.{start}"
    if op == "del":
        return cds[:i] + cds[j:], f"deletion of {j - i} nt at c.{start}"
    if op == "dup":
        return cds[:j] + cds[i:j] + cds[j:], f"duplication of {j - i} nt at c.{start}"
    if op == "ins":
        if not seq:
            return None
        return cds[:end] + seq + cds[end:], f"insertion of {len(seq)} nt after c.{start}"
    if op == "delins":
        if not seq:
            return None
        return cds[:i] + seq + cds[j:], f"{j - i} nt replaced by {len(seq)} nt at c.{start}"
    return None


def altered_protein(cds: str, hgvsc: str) -> dict[str, Any] | None:
    """Translate the reference and altered coding sequences and compare them.

    The novel stretch is everything from the first residue that differs to the
    new stop codon. For a frameshift that is usually the richest source of
    tumour-specific peptides in a genome, because none of it exists in any
    healthy cell.
    """
    from genomeos.runtime.central_dogma import translate

    applied = apply_hgvs_c(cds, hgvsc)
    if applied is None:
        return None
    mutant_cds, description = applied
    reference = translate(cds)
    mutant = translate(mutant_cds)
    shortest = min(len(reference), len(mutant))
    first = next((i for i in range(shortest) if reference[i] != mutant[i]), shortest)
    # realign from the end, or a single substitution looks like a frameshift
    suffix = 0
    while (
        suffix < shortest - first
        and reference[len(reference) - 1 - suffix] == mutant[len(mutant) - 1 - suffix]
    ):
        suffix += 1
    novel = mutant[first : len(mutant) - suffix]
    return {
        "change": description,
        # the frame, not the protein length, decides whether the mutant still
        # aligns to the reference: a frameshift with no downstream stop can
        # produce a protein of exactly the same length and share no sequence
        "frame_shifted": (len(mutant_cds) - len(cds)) % 3 != 0,
        "reference_length": len(reference),
        "mutant_length": len(mutant),
        "first_altered_residue": first + 1,
        "novel_stretch": novel,
        "novel_residues": len(novel),
        "novel_end_residue": len(mutant) - suffix,
        "reference_protein": reference,
        "mutant_protein": mutant,
        "truncated": len(mutant) < len(reference),
    }


def peptides_over_novel(
    mutant: str,
    first_altered: int,
    novel_end: int | None = None,
    lengths: tuple[int, ...] = PEPTIDE_LENGTHS,
    source: str = "frameshift_derived",
) -> list[PeptideCandidate]:
    """Class-I windows of the altered protein that cover the changed region.

    An in-frame deletion can leave no residue that did not exist before, and
    still be tumour-specific: the two residues either side of it are now
    adjacent, and a peptide spanning that junction exists in no healthy cell.
    So the window covers the changed span, or the junction itself when the span
    is empty.
    """
    out: list[PeptideCandidate] = []
    lo = first_altered - 1
    hi = max(novel_end if novel_end is not None else first_altered, first_altered)
    for n in lengths:
        last_start = len(mutant) - n
        if last_start < 0:
            continue
        for start in range(max(0, lo - n + 1), min(hi, last_start + 1)):
            end = start + n
            out.append(
                PeptideCandidate(
                    sequence=mutant[start:end],
                    wild_type_sequence="",
                    length=n,
                    start=start + 1,
                    end=end,
                    mutation_offset=lo - start + 1,
                    source=source,
                )
            )
    return out


def _reconstruct(
    gene: str,
    consequence: str,
    hgvsc: str,
    transcript: str,
    transcripts: Any,
    a: NeoantigenAssessment,
) -> dict[str, Any] | None:
    """Translate the altered transcript and read the answer off the sequence."""
    if transcripts is None or not transcript or not hgvsc:
        return None
    answer = transcripts.cds(transcript)
    if not answer.available:
        a.missing.append(
            Missing("altered protein sequence", answer.reason, ("the transcript's coding sequence",))
        )
        return None
    rebuilt = altered_protein(answer.data, hgvsc)
    if rebuilt is None:
        a.missing.append(
            Missing(
                "altered protein sequence",
                f"the coding change {hgvsc} is outside the coding sequence or in a form GenomeOS does "
                "not reconstruct",
                ("a transcript-level caller that emits the altered protein directly",),
            )
        )
        return None
    a.evidence.extend(answer.evidence)
    novel = rebuilt["novel_residues"]
    inframe = consequence in ("inframe_deletion", "inframe_insertion")
    a.novel_peptide_sequence = bool(novel) or inframe
    if novel:
        a.reason = (
            f"the altered {transcript} translates to {rebuilt['mutant_length']} residues against "
            f"{rebuilt['reference_length']}, diverging at residue {rebuilt['first_altered_residue']} "
            f"and producing {novel} residue{'s that exist' if novel != 1 else ' that exists'} in no "
            "reference protein"
        )
    elif inframe:
        a.reason = (
            f"the alteration removes or adds residues without shifting the frame, so no residue is new, "
            f"but the sequence either side of residue {rebuilt['first_altered_residue']} is now adjacent "
            "in a way it is not in any healthy cell"
        )
    else:
        a.reason = (
            f"the altered {transcript} translates to {rebuilt['mutant_length']} residues against "
            f"{rebuilt['reference_length']} with no residue that differs from the reference: a "
            "truncation removes protein rather than creating new sequence"
        )
        a.applicable = False
        a.evidence.append(
            derived(
                "GenomeOS neoantigen reasoning",
                f"{gene}: reconstructed from {transcript} and confirmed by sequence, not assumed from "
                "the variant class",
                0.9,
            )
        )
        return rebuilt

    shifted = rebuilt["frame_shifted"]
    substitution = not shifted and rebuilt["mutant_length"] == rebuilt["reference_length"]
    peptides = peptides_over_novel(
        rebuilt["mutant_protein"],
        rebuilt["first_altered_residue"],
        rebuilt["novel_end_residue"],
        source=(
            "frameshift_derived" if shifted else "mutation_derived" if substitution else "junction_derived"
        ),
    )
    reference = rebuilt["reference_protein"]
    if substitution:
        for p in peptides:
            if p.end <= len(reference):
                p.wild_type_sequence = reference[p.start - 1 : p.end]
    a.peptides = peptides
    a.applicable = True
    a.evidence.append(
        derived(
            "GenomeOS neoantigen reasoning",
            f"{len(peptides)} peptides covering the altered region of {gene}, derived by applying "
            f"{hgvsc} to the coding sequence of {transcript} and translating it",
            0.75,
        )
    )
    return rebuilt
