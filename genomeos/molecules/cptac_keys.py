# SPDX-License-Identifier: AGPL-3.0-or-later
"""Gene-keyed CPTAC phosphosites: the gene key and the isoform position, each resolved by a registered rule.

Seven CPTAC phospho profiles on cBioPortal name a site by gene symbol and a position on an unnamed
RefSeq isoform. Taken on the UniProt canonical of the one compiled entry for the symbol, the residue
letter disagreed for 31,424 of 219,035 site entities (14.3%) on 2026-09-28, so the profiles were counted
and not committed (docs/PROTEIN.md, "Relative abundance in tumours (2026-09-28)"). The 14.3% is a residue
mismatch, an isoform problem; separately, 8,085 site entities were dropped because their symbol named no
compiled entry, a gene-key problem. This module holds the rule for both, registered before it was run
(``REGISTRATION`` below; docs/PROTEIN.md, "Gene-keyed profiles (2026-10-10)"):

- gene key: HGNC approved symbol, then HGNC previous symbol, then HGNC alias symbol, each only when it names
  exactly one HGNC record and no earlier tier claims the symbol for another record; the record's UniProt
  accessions among the compiled definitions, else its Ensembl gene among them, must give exactly one
  entry. Anything else is ambiguous (more than one) or unresolved (none). Never a guess.
- position: every candidate isoform of the entry (the UniProt canonical plus GENCODE v50 complete-CDS
  translations of the entry's Ensembl gene) that carries the residue letter at the position is carried
  onto the canonical with ``ptm.transfer_position``. One canonical position from all of them: resolved.
  Disagreement: ambiguous, refused. No candidate carries the letter: unresolved.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from genomeos.molecules import ptm

GENE_KEYED_PROFILES = (
    "brca_cptac_2020_phosphoproteome",
    "ucec_cptac_2020_phosphoproteome",
    "coad_cptac_2019_phosphoprotein_quantification",
    "brca_tcga_phosphoprotein_quantification",
    "ov_tcga_phosphoprotein_quantification",
    "brca_tcga_pan_can_atlas_2018_phosphoprotein_quantification",
    "ov_tcga_pan_can_atlas_2018_phosphoprotein_quantification",
)
REFSEQ_KEYED_PROFILES = (
    "luad_cptac_2020_phosphoproteome",
    "gbm_cptac_2021_phosphoproteome",
    "paad_cptac_2021_phosphoproteome",
    "brain_cptac_2020_phosphoprotein",
)
RESULT = Path("data/results/cptac_genekeyed_result.json")
REGISTRATION_RESULT = Path("data/results/cptac_genekeyed_registration.json")

# Registered 2026-10-10 by lane-cptac, after the 14.3% was reproduced from the cached mapping and before
# the rule below was run on any site, gene-keyed or RefSeq-keyed. Inputs fixed at registration: the
# cBioPortal entity lists and the RefSeq-keyed mapping cached on 2026-09-28 (data/cache/cptac/), the HGNC
# complete set fetched 2026-10-10 (sha256 in the registration result), GENCODE v50 and the compiled
# UniProt definitions (data/knowledge/proteins).
REGISTRATION: dict[str, Any] = {
    "question": (
        "Can the seven gene-keyed CPTAC phospho profiles be committed beside the three RefSeq-keyed ones, "
        "with their gene keys and isoform positions resolved by a fixed rule and never guessed?"
    ),
    "what_the_14_3_percent_was": (
        "a residue-letter mismatch on the UniProt canonical (31,424 of 219,035 single-site entities whose "
        "symbol named one compiled entry), i.e. positions on another isoform; not a gene-key mismatch. The "
        "gene-key drop was separate: 8,085 single-site entities whose symbol named no compiled entry."
    ),
    "gene_key_rule": [
        "1. the symbol is an HGNC approved symbol, and is no other record's previous symbol",
        "2. else the symbol is the previous symbol of exactly one HGNC record",
        "3. else the symbol is an alias of exactly one HGNC record",
        "a symbol that is approved for one record and previous for another, or previous or alias for "
        "several records, is ambiguous and refused; a symbol HGNC does not know is unresolved",
        "the record's UniProt accessions among the compiled definitions must be exactly one entry; with "
        "none, the record's Ensembl gene among the compiled definitions' genes must be exactly one entry; "
        "more than one: ambiguous, refused; none: unresolved (no compiled definition)",
    ],
    "candidate_isoforms": (
        "the UniProt canonical sequence of the compiled entry, plus every GENCODE v50 translation of the "
        "entry's Ensembl gene (its compiled genomic_origin gene, else the HGNC record's) whose CDS is "
        "complete: no cds_start_NF, cds_end_NF, mRNA_start_NF or mRNA_end_NF tag, CDS length a multiple of "
        "three after the phase, translation reaching the final codon (a stop) with no earlier stop"
    ),
    "position_rule": (
        "every candidate that carries the entity's residue letter at its position is carried onto the "
        "canonical with ptm.transfer_position (identical sequence, else the 15-residue window occurring "
        "exactly once), the canonical letter checked; one canonical position from every carrier: "
        "resolved; carriers that disagree, or one carrier with a canonical position and another without: "
        "ambiguous, refused; carriers none of which has a canonical counterpart: carried off the canonical "
        "(not joinable, not counted as unresolved); no carrier: unresolved"
    ),
    "duplicates": "several entities on one canonical site in one profile: all dropped, as on 2026-09-28",
    "calibration": (
        "the same rule, symbol and position only, on the RefSeq-keyed sites whose canonical position is "
        "known from their named RefSeq protein (data/cache/cptac/mapped.json); wrong = resolved to "
        "another accession or position"
    ),
    "decoy": (
        "each gene-keyed single-site entity with a resolved gene key, its position moved by -13 and +13 "
        "residues (same letter), run through the position rule; the estimated wrong resolutions in a "
        "profile are its unresolved sites times (decoys resolved / decoys unresolved), i.e. the chance a "
        "site whose isoform is absent from the candidates is carried by a candidate at that position anyway"
    ),
    "decoy_offsets": (-13, 13),
    "gates": {
        "G0_rule_wrong_rate_on_refseq_truth_below": 0.02,
        "G1_median_of_site_medians_within": (-0.5, 0.5),
        "G2_unresolved_rate_below": 0.02,
        "G3_estimated_wrong_rate_below": 0.02,
    },
    "gates_meaning": {
        "G0": (
            "pooled over the four RefSeq-keyed profiles; if it fails the rule is not fit and all seven "
            "profiles are retired"
        ),
        "G1": (
            "the 2026-09-28 units check, recomputed on the profile's resolved sites on curated phospho "
            "sites; known before this registration: brca_tcga_pan_can_atlas_2018 sat at -0.518 under the "
            "canonical-only mapping"
        ),
        "G2": (
            "unresolved / single-site entities with a resolved gene key; the bound is the 2026-09-28 "
            "gene_keyed_mismatch_bound (0.02), unchanged, and the unresolved count is that registration's "
            "mismatch measured against every candidate isoform instead of the canonical alone"
        ),
        "G3": "estimated wrong resolutions / sites resolved to a canonical position",
    },
    "not_gated": (
        "gene keys that are ambiguous or unresolved cost coverage, not correctness: every one is listed with "
        "its entity count, and its share is reported; ambiguous sites are refused and counted"
    ),
    "committed": (
        "a profile passing G0 to G3 is committed: its resolved sites on curated phospho sites get "
        "relative_abundance_in_tumours from that study, as for the three RefSeq-keyed profiles (tumours "
        "with a value, median log2 ratio to the study's pooled reference; values summarised in memory, "
        "never per tumour in the repository); a profile failing any gate is retired: counted, not committed"
    ),
    "wording_committed": (
        "{profile}: committed. Gene keys by the registered HGNC rule; positions carried from the isoform "
        "that carries the residue onto the UniProt canonical: {resolved} of {sites} single-site entities "
        "resolved, {unresolved} unresolved ({unresolved_rate}), {ambiguous} ambiguous, estimated wrong "
        "resolutions {estimated_wrong_rate}. Values are relative abundance measured in that study's "
        "tumours, not occupancy, not in a person."
    ),
    "wording_retired": (
        "{profile}: retired, not committed: {failed}. Its sites stay counted only; the gene-keyed "
        "positions cannot be placed on an isoform under the registered rule closely enough to commit."
    ),
    "expected": (
        "colon 2019 (canonical mismatch 4.3%) the likeliest to pass G2; breast 2020 and endometrial 2020 "
        "uncertain; the four TCGA quantifications (canonical mismatch 16.2% to 19.2%, RefSeq of about 2013) "
        "likely to fail G2 because GENCODE v50 lacks some of their isoforms; TCGA breast PanCan likely to "
        "fail G1 as before. No coverage expectation is registered."
    ),
}


@dataclass
class Hgnc:
    """HGNC complete set, indexed by approved, previous and alias symbol."""

    approved: dict[str, str] = field(default_factory=dict)  # symbol -> hgnc_id
    previous: dict[str, set[str]] = field(default_factory=lambda: defaultdict(set))
    alias: dict[str, set[str]] = field(default_factory=lambda: defaultdict(set))
    uniprot: dict[str, list[str]] = field(default_factory=dict)  # hgnc_id -> accessions
    ensembl: dict[str, str] = field(default_factory=dict)  # hgnc_id -> ENSG
    symbol: dict[str, str] = field(default_factory=dict)  # hgnc_id -> approved symbol

    @classmethod
    def from_rows(cls, rows: Iterable[dict[str, str]]) -> Hgnc:
        h = cls()

        def split(v: str | None) -> list[str]:
            return [x for x in (v or "").strip('"').split("|") if x]

        for r in rows:
            if (r.get("status") or "Approved") != "Approved":
                continue
            hid = r["hgnc_id"]
            h.approved[r["symbol"]] = hid
            h.symbol[hid] = r["symbol"]
            h.uniprot[hid] = split(r.get("uniprot_ids"))
            if r.get("ensembl_gene_id"):
                h.ensembl[hid] = r["ensembl_gene_id"]
            for s in split(r.get("prev_symbol")):
                h.previous[s].add(hid)
            for s in split(r.get("alias_symbol")):
                h.alias[s].add(hid)
        return h

    @classmethod
    def from_tsv(cls, path: str | Path) -> Hgnc:
        with open(path, newline="") as fh:
            return cls.from_rows(csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE))

    def record(self, symbol: str) -> tuple[str, str | None, str]:
        """(status, hgnc_id, how): status is resolved, ambiguous or unresolved."""
        prev = self.previous.get(symbol, set())
        if symbol in self.approved:
            hid = self.approved[symbol]
            if prev - {hid}:
                return "ambiguous", None, "approved for one record and previous symbol of another"
            return "resolved", hid, "approved symbol"
        if len(prev) == 1:
            return "resolved", next(iter(prev)), "previous symbol"
        if len(prev) > 1:
            return "ambiguous", None, f"previous symbol of {len(prev)} records"
        al = self.alias.get(symbol, set())
        if len(al) == 1:
            return "resolved", next(iter(al)), "alias symbol"
        if len(al) > 1:
            return "ambiguous", None, f"alias of {len(al)} records"
        return "unresolved", None, "symbol unknown to HGNC"


def resolve_gene_key(
    symbol: str,
    hgnc: Hgnc,
    compiled: set[str],
    compiled_by_ensembl: dict[str, list[str]],
) -> tuple[str, str | None, str]:
    """(status, compiled accession or None, how) for one gene symbol, by the registered rule."""
    status, hid, how = hgnc.record(symbol)
    if status != "resolved" or hid is None:
        return status, None, how
    accs = sorted({a for a in hgnc.uniprot.get(hid, []) if a in compiled})
    if len(accs) == 1:
        return "resolved", accs[0], how + ", UniProt accession"
    if len(accs) > 1:
        return "ambiguous", None, f"{how}, {len(accs)} compiled UniProt entries"
    ens = hgnc.ensembl.get(hid)
    by_ens = sorted(set(compiled_by_ensembl.get(ens, []))) if ens else []
    if len(by_ens) == 1:
        return "resolved", by_ens[0], how + ", Ensembl gene"
    if len(by_ens) > 1:
        return "ambiguous", None, f"{how}, {len(by_ens)} compiled entries for the Ensembl gene"
    return "unresolved", None, f"{how}, no compiled definition for the HGNC record"


def resolve_site(
    residue: str, position: int, canonical: str, candidates: Iterable[str], flank: int = 7
) -> tuple[str, int | None, str]:
    """(status, canonical position or None, how) for one gene-keyed site, by the registered rule.

    ``candidates`` are isoform sequences of the entry; the canonical is always added. Status is resolved,
    ambiguous, off_canonical or unresolved.
    """
    seqs = {canonical, *candidates}
    carriers = [s for s in seqs if 0 < position <= len(s) and s[position - 1] == residue]
    if not carriers:
        return "unresolved", None, "no candidate isoform carries the residue at the position"
    targets: set[int | None] = set()
    hows = set()
    for s in carriers:
        pos, how = ptm.transfer_position(s, canonical, position, flank)
        if pos is None or canonical[pos - 1] != residue:
            targets.add(None)
        else:
            targets.add(pos)
            hows.add("canonical" if s == canonical else "isoform window")
    if len(targets) > 1:
        return "ambiguous", None, f"{len(carriers)} carrying isoforms disagree on the canonical position"
    t = next(iter(targets))
    if t is None:
        return "off_canonical", None, "carried only by isoform sequence with no canonical counterpart"
    return "resolved", t, "+".join(sorted(hows))


def estimated_wrong(unresolved: int, decoy_resolved: int, decoy_unresolved: int) -> float | None:
    """Expected chance resolutions: unresolved sites times the decoy odds of resolving by chance."""
    if decoy_unresolved == 0:
        return None
    return unresolved * decoy_resolved / decoy_unresolved
