# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""An individual-level panel for the polygenic scores lane-pgs left UNKNOWN (ROADMAP area D).

lane-pgs (`66a214b` code, `e9f0c95` registered, `0d19277` result) reads every polygenic score as a
position among the person's declared group, modelled from gnomAD v4 group frequencies under
Hardy-Weinberg and independence between variants. On the HPRC panel the independence variance
holds for PGS000004 (PRS313, 1.06 [0.70, 1.77]) and PGS000115 (LDL-C, 0.63 [0.42, 1.06]) and FAILS
for PGS000297 (height, 3,290 SNPs, 1.53 [1.016, 2.56]) and PGS000018 (metaGRS, 1,745,158 variants,
3.19 [2.12, 5.34]), so for those two the person reads UNKNOWN. The remedy is an empirical
distribution: the same score computed person by person on real genotypes of people of the declared
group, which carries the linkage disequilibrium by construction and needs no independence model.

This module fixes, BEFORE any candidate panel's population counts are read and before any score is
computed on any panel, what a panel must meet (GATE) and which candidates are examined
(CANDIDATES). `scripts/pgs_panel_inventory.py` reads each candidate's own sample metadata and
applies `evaluate`. When no candidate meets the gate for the person's declared group, the person
keeps the UNKNOWN lane-pgs gave: the lane stops with the inventory as its result, and no person is
ever re-labelled into a group a panel does cover.

A percentile is a position among people of the stated group on one score; it is not a risk.
"""

from __future__ import annotations

import re
from typing import Any

from genomeos.genome.pgs_reference import DECLARED, GNOMAD, GNOMAD_LABEL
from genomeos.genome.pgs_reference import RULES as REFERENCE_RULES

PERSON = "HG002"
# HG002's first-degree relatives in the GIAB trio: never counted in the person's own panel.
RELATIVES = {"HG002": ("HG003", "HG004", "NA24149", "NA24143", "NA24385")}

LEFT_UNKNOWN = {
    "PGS000297": "height, 3,290 SNPs: R1 1.53 [1.016, 2.56] on the HPRC panel, independence fails",
    "PGS000018": "metaGRS (coronary artery disease), 1,745,158 variants: R1 3.19 [2.12, 5.34], "
    "independence fails, and above the 20,000-variant per-variant lookup cap",
}
CROSS_CHECK = {
    "PGS000004": "PRS313 (breast cancer): R1 1.06 [0.70, 1.77], independence holds",
    "PGS000115": "LDL-C, 223 variants: R1 0.63 [0.42, 1.06], independence holds",
}
REFERENCE_RESULT = "data/results/pgs_reference.json"

# Frozen 2026-10-10, before any candidate's metadata is read and before any score is computed.
GATE: dict[str, Any] = {
    "min_people": 100,
    "min_people_reads": "at least 100 people carrying the declared group's label in the panel, "
    "unrelated to each other and to the person. At 100 people an empirical percentile near the "
    "middle carries a 95% binomial interval of about +/-10 points (1.96 * sqrt(0.25 / 100) = 0.098), "
    "and the chi-square interval of a variance read from them spans about [0.77, 1.35] times the "
    "estimate (df 99), against [0.70, 1.77] for the 39 HPRC people lane-pgs had (df 37), whose "
    "lower bound for height-3290 sat at 1.016",
    "access": "open",
    "access_reads": "the per-person genotypes download without an application, a data-access "
    "committee, a signed letter, a login, or a request to the authors, under terms that permit "
    "research use; anything else fails the gate whatever its size",
    "individual_level": "per-person genotypes (sequence or array calls) are released; aggregate "
    "allele frequencies, which lane-pgs already uses, do not qualify",
    "label": "a panel person belongs to the declared group only through the panel's OWN population "
    "label (the sampling population the panel publishes): for asj, a label naming Ashkenazi Jews "
    "(ASJ_LABEL); for the other gnomAD groups, the frozen 1000 Genomes code table THOUSAND_GROUP. "
    "No person, panel or HG002, is assigned to a group by genetic similarity, and HG002's declared "
    "group (asj, as its donors declared) is never replaced by one a panel covers",
    "unrelated": "people are counted once per family: a panel's own unrelated set where it publishes "
    "one (1000 Genomes phase 3's 2,504), founders only in a trio, and the person and the person's "
    "first-degree relatives (RELATIVES) never",
    "coverage_floor": REFERENCE_RULES["coverage_floor"],
    "coverage": "evaluated only for a candidate that meets access, individual level, label and "
    "size: the score's variants with a directly called genotype in the panel's release (no "
    "imputation) must carry at least coverage_floor of the score's total |weight|, the bar "
    "lane-pgs used",
    "stop": "if no candidate meets the gate for the person's declared group, the lane stops and the "
    "inventory is the result: nothing is registered, no score is computed on any panel, and the "
    "person's UNKNOWN for PGS000297 and PGS000018 stands",
    "other_groups": "for every other gnomAD group the inventory records which candidates meet "
    "access, individual level, label and size (coverage not evaluated): what a later registration "
    "could use for a person who declares that group. Descriptive; no score is computed",
}

ACCESS_LEVELS = ("open", "signed_letter", "controlled", "on_request", "unreachable")

ASJ_LABEL = re.compile(r"ashkenaz", re.IGNORECASE)

# 1000 Genomes phase 3 population codes read against gnomAD's groups (the GATE's other_groups).
THOUSAND_GROUP = {
    **dict.fromkeys(("YRI", "LWK", "GWD", "MSL", "ESN", "ASW", "ACB"), "afr"),
    **dict.fromkeys(("MXL", "PUR", "CLM", "PEL"), "amr"),
    **dict.fromkeys(("CHB", "JPT", "CHS", "CDX", "KHV"), "eas"),
    **dict.fromkeys(("GIH", "PJL", "BEB", "STU", "ITU"), "sas"),
    "FIN": "fin",
    **dict.fromkeys(("CEU", "GBR", "TSI", "IBS"), "nfe"),
}

# The candidates, fixed before their metadata is read. `metadata` is what the inventory script
# fetches for the per-population counts; `count_from` says where the count comes from when there is
# no machine-readable sample sheet.
CANDIDATES: tuple[dict[str, Any], ...] = (
    {
        "key": "1kgp3",
        "name": "1000 Genomes Project phase 3 (IGSR), 2,504 unrelated people; the same people "
        "re-sequenced at 30x by NYGC",
        "data": "whole-genome sequence calls",
        "access": "open",
        "terms": "IGSR: no access restriction (Fort Lauderdale data-release principles)",
        "metadata": "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
        "integrated_call_samples_v3.20130502.ALL.panel",
        "count_from": "metadata",
    },
    {
        "key": "hgdp",
        "name": "HGDP-CEPH whole genomes (Bergstrom et al. 2020, Science 367:eaay5012)",
        "data": "whole-genome sequence calls",
        "access": "open",
        "terms": "public release on the Sanger and ENA servers, no application",
        "metadata": "https://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/metadata/"
        "hgdp_wgs.20190516.metadata.txt",
        "count_from": "metadata",
    },
    {
        "key": "sgdp",
        "name": "Simons Genome Diversity Project (Mallick et al. 2016, Nature 538:201)",
        "data": "whole-genome sequence calls",
        "access": "open",
        "terms": "279 samples public; 21 under a signed letter and 44 (Fan et al.) under a signed "
        "letter with delay, which fail the access gate",
        "metadata": "https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/"
        "SGDP_metadata.279public.21signedLetter.44Fan.samples.txt",
        "count_from": "metadata",
    },
    {
        "key": "hprc2",
        "name": "Human Pangenome Reference Consortium release 2",
        "data": "phased diploid assemblies",
        "access": "open",
        "terms": "HPRC data release: open access on AWS and AnVIL",
        "metadata": "https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/"
        "main/data_tables/sample/hprc_release2_sample_metadata.csv",
        "count_from": "metadata",
    },
    {
        "key": "giab",
        "name": "Genome in a Bottle Ashkenazi trio (HG002 son, HG003 father, HG004 mother)",
        "data": "whole-genome sequence calls",
        "access": "open",
        "terms": "NIST GIAB, public domain; the donors are Personal Genome Project participants",
        "metadata": None,
        "count_from": "the trio's published pedigree: the person and the person's two parents",
    },
    {
        "key": "gnomad_hgdp_1kg",
        "name": "gnomAD v3.1.2 HGDP + 1000 Genomes subset (the same people as 1kgp3 and hgdp, jointly "
        "called), with gnomAD's inferred group per sample",
        "data": "whole-genome sequence calls",
        "access": "open",
        "terms": "gnomAD release on public cloud buckets, no application",
        "metadata": "https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/"
        "gnomad.genomes.v3.1.2.hgdp_1kg_subset_sample_meta.tsv.bgz",
        "count_from": "metadata (the projects' own population labels; gnomAD's inferred group is "
        "recorded beside them as description, never as a label)",
    },
    {
        "key": "aadr_ho",
        "name": "Allen Ancient DNA Resource v66.p1, present-day people on the Affymetrix Human "
        "Origins array (public subset)",
        "data": "array genotypes (about 600,000 SNPs)",
        "access": "open",
        "terms": "Harvard Dataverse doi:10.7910/DVN/FFIDCW, public (.PUB) files",
        "metadata": "https://dataverse.harvard.edu/api/access/datafile/13994526",
        "count_from": "metadata (the .ind file's population column)",
    },
    {
        "key": "ebc_jew2010",
        "name": "Estonian Biocentre free data: Behar et al. 2010, Nature 466:238, Jewish populations",
        "data": "array genotypes (Illumina)",
        "access": "open",
        "terms": "Estonian Biocentre free data repository, download without registration",
        "metadata": "https://evolbio.ut.ee/jew/",
        "count_from": "metadata (the released sample file's population column)",
    },
    {
        "key": "ebc_khazar2013",
        "name": "Estonian Biocentre free data: Behar et al. 2013, Human Biology 85:859, Ashkenazi origins",
        "data": "array genotypes (Illumina)",
        "access": "open",
        "terms": "Estonian Biocentre free data repository, download without registration",
        "metadata": "https://evolbio.ut.ee/khazar/",
        "count_from": "metadata (the released sample file's population column)",
    },
    {
        "key": "tagc",
        "name": "The Ashkenazi Genome Consortium phase I (Carmi et al. 2014, Nat Commun 5:4835), "
        "EGAD00001000781",
        "data": "whole-genome sequence calls (Complete Genomics, 50x)",
        "access": "controlled",
        "terms": "European Genome-phenome Archive, controlled access through data-access committee "
        "EGAC00001000151 (policy EGAP00001000190)",
        "metadata": "https://metadata.ega-archive.org/studies/EGAS00001000664/datasets",
        "count_from": "metadata (EGA's dataset record: num_samples and access_type)",
    },
    {
        "key": "pgp",
        "name": "Harvard Personal Genome Project public profiles",
        "data": "mixed: whole genomes and consumer arrays uploaded per participant",
        "access": "open",
        "terms": "CC0 per participant",
        "metadata": "https://my.pgp-hms.org/google_surveys/1",
        "count_from": "metadata (the participant survey's ancestry answers), when the site answers",
    },
    {
        "key": "controlled_cohorts",
        "name": "application-controlled cohorts holding Ashkenazi Jewish people (dbGaP studies, UK "
        "Biobank, All of Us, TOPMed)",
        "data": "array or sequence",
        "access": "controlled",
        "terms": "each needs an access application; excluded by the brief",
        "metadata": None,
        "count_from": "not counted: excluded by access before any count matters",
    },
)


def group_of(panel: str, label: str) -> str | None:
    """The gnomAD group a panel's own population label names, or None. Never by genetic similarity."""
    if ASJ_LABEL.search(label or ""):
        return "asj"
    if panel in ("1kgp3", "gnomad_hgdp_1kg"):
        return THOUSAND_GROUP.get((label or "").strip())
    return None


def counted(people: list[dict[str, Any]], person: str = PERSON) -> list[dict[str, Any]]:
    """The people a panel contributes: unrelated, never the person or the person's relatives.

    Each entry carries `id`, `label`, and optionally `related` (True for a non-founder or a member
    of a related pair the panel itself flags)."""
    excluded = {person, *RELATIVES.get(person, ())}
    return [p for p in people if p["id"] not in excluded and not p.get("related")]


def group_counts(panel: str, people: list[dict[str, Any]], person: str = PERSON) -> dict[str, int]:
    """People per gnomAD group, through `group_of` on the panel's own label, after `counted`."""
    out = dict.fromkeys(GNOMAD, 0)
    for p in counted(people, person):
        g = group_of(panel, p["label"])
        if g is not None:
            out[g] += 1
    return out


def evaluate(candidate: dict[str, Any], group: str, n_people: int | None) -> dict[str, Any]:
    """Read one candidate against GATE for one group. Coverage is left to the registration."""
    if group not in GNOMAD:
        raise ValueError(f"not a gnomAD group: {group}")
    fails = []
    if candidate["access"] != GATE["access"]:
        fails.append(f"access: {candidate['access']} ({candidate['terms']})")
    if "frequenc" in candidate["data"] and "genotype" not in candidate["data"]:
        fails.append("individual_level: aggregate frequencies only")
    if n_people is None:
        fails.append("size: not counted")
    elif n_people < GATE["min_people"]:
        fails.append(f"size: {n_people} {group} people, fewer than {GATE['min_people']}")
    return {
        "candidate": candidate["key"],
        "group": group,
        "group_label": GNOMAD_LABEL[group],
        "people": n_people,
        "meets_size_access_label": not fails,
        "fails": fails,
        "coverage": "not evaluated: evaluated only at registration, for a candidate that meets "
        "the rest of the gate",
    }


def outcome(evaluations: list[dict[str, Any]], group: str) -> str:
    """PANEL_QUALIFIES when some candidate meets size, access and label for the group; else NONE."""
    ok = [e for e in evaluations if e["group"] == group and e["meets_size_access_label"]]
    return "PANEL_QUALIFIES" if ok else "NO_PANEL_QUALIFIES"


def declared_group(person: str = PERSON) -> str:
    """The person's group as declared (lane-pgs's DECLARED), never inferred."""
    return DECLARED[person]
