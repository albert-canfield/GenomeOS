"""BioLib catalogue v0.1: the reusable "libraries" that exist in the human genome.

This is the reverse-engineering map. Each Library names a conserved toolkit,
its layer, what it does, representative human genes (HGNC symbols) and the
evidence source. Gene lists are representative, not exhaustive; Phase 1 of
the action plan replaces them with GENCODE + GO membership computed by the
compiler.

Layers follow the four questions the project asks of the genome:
    core        what every cell needs to exist and copy itself
    blueprint   how the body plan and each organ is specified
    timer       how construction and maintenance are scheduled and aged
    systems     how tissues and organs communicate and work together
    parts       the inventory of cell types, tissues and organs
"""

from __future__ import annotations

from dataclasses import dataclass

LAYERS = ("core", "blueprint", "timer", "systems", "parts")


@dataclass(frozen=True, slots=True)
class Library:
    id: str  # e.g. "core.replication"
    layer: str
    purpose: str
    genes: tuple[str, ...]  # representative HGNC symbols
    scale: str = ""  # size of the real library in the genome
    source: str = ""  # where the membership comes from
    note: str = ""
    go_terms: tuple[str, ...] = ()  # GO ids whose annotated genes (with descendants) define membership
    reactome: tuple[str, ...] = ()  # substrings of Reactome pathway names that also define membership
    noncoding: tuple[str, ...] = ()  # RNA genes in `genes` that protein-centric GO annotation does not cover


LIBRARIES: dict[str, Library] = {}


def _add(lib: Library) -> None:
    assert lib.layer in LAYERS, lib.layer
    LIBRARIES[lib.id] = lib


# ----------------------------------------------------------------- core ---
_add(
    Library(
        "core.replication",
        "core",
        "copy the genome once per cell cycle",
        ("ORC1", "CDC6", "CDT1", "MCM2", "MCM7", "POLA1", "POLD1", "POLE", "PCNA", "LIG1"),
        "~100 genes",
        "GO:0006260 DNA replication",
        go_terms=("GO:0006260",),
    )
)
_add(
    Library(
        "core.transcription",
        "core",
        "read DNA into RNA",
        ("POLR2A", "TBP", "GTF2B", "MED1", "TAF1", "POLR1A", "POLR3A"),
        "~300 genes incl. Mediator and general TFs",
        "GO:0006351 transcription",
        go_terms=("GO:0006351",),
        reactome=("Generic Transcription Pathway",),
    )
)
_add(
    Library(
        "core.splicing",
        "core",
        "cut introns, join exons",
        ("SF3B1", "U2AF1", "SNRNP70", "PRPF8", "SRSF1", "HNRNPA1"),
        "~300 spliceosome and RNA-binding genes",
        "GO:0000398 mRNA splicing via spliceosome",
        go_terms=("GO:0000398",),
    )
)
_add(
    Library(
        "core.translation",
        "core",
        "read RNA into protein",
        ("RPL3", "RPS6", "EIF4E", "EIF2S1", "EEF1A1", "EEF2"),
        "~80 ribosomal proteins, ~500 tRNA genes, rRNA arrays on chr13/14/15/21/22",
        "GO:0006412 translation",
        go_terms=("GO:0006412",),
        reactome=("Translation",),
    )
)
_add(
    Library(
        "core.folding_proteostasis",
        "core",
        "fold proteins and keep them folded",
        ("HSPA1A", "HSPA8", "HSP90AA1", "HSPD1", "CCT2", "DNAJB1", "HSF1"),
        "~100 chaperones",
        "GO:0006457 protein folding",
        go_terms=("GO:0006457", "GO:0009408"),
    )
)
_add(
    Library(
        "core.degradation",
        "core",
        "remove proteins and organelles: proteasome and autophagy",
        ("UBB", "UBC", "PSMA1", "PSMB5", "ATG5", "ATG7", "BECN1", "SQSTM1"),
        "~600 ubiquitin-system genes",
        "GO:0043161, GO:0006914",
        go_terms=(
            "GO:0006511",
            "GO:0006914",
            "GO:0016567",
        ),
    )
)
_add(
    Library(
        "core.dna_repair",
        "core",
        "maintenance of the genome itself",
        (
            "ATM",
            "ATR",
            "TP53",
            "BRCA1",
            "BRCA2",
            "RAD51",
            "XRCC4",
            "LIG4",
            "MLH1",
            "MSH2",
            "XPA",
            "ERCC2",
            "OGG1",
            "APEX1",
            "PARP1",
        ),
        "~450 genes across BER, NER, MMR, HR, NHEJ",
        "GO:0006281 DNA repair",
        "loss-of-function causes progeroid syndromes: this library IS the maintenance schedule",
        go_terms=("GO:0006281",),
        reactome=("DNA Repair",),
    )
)
_add(
    Library(
        "core.cell_cycle",
        "core",
        "decide whether and when to divide",
        (
            "CCND1",
            "CDK4",
            "CCNE1",
            "CDK2",
            "CCNA2",
            "CCNB1",
            "CDK1",
            "RB1",
            "E2F1",
            "CDKN1A",
            "CDKN1B",
            "CDKN2A",
            "CDC20",
            "PLK1",
            "AURKB",
        ),
        "~600 genes",
        "GO:0007049 cell cycle",
        go_terms=("GO:0007049", "GO:0051726"),
        reactome=("Cell Cycle",),
    )
)
_add(
    Library(
        "core.apoptosis",
        "core",
        "programmed cell death",
        ("CASP3", "CASP8", "CASP9", "BAX", "BAK1", "BCL2", "BCL2L1", "CYCS", "APAF1", "FAS"),
        "~150 genes",
        "GO:0006915 apoptotic process",
        go_terms=("GO:0006915",),
        reactome=("Apoptosis",),
    )
)
_add(
    Library(
        "core.energy",
        "core",
        "glycolysis, TCA cycle, oxidative phosphorylation",
        (
            "HK1",
            "PFKM",
            "GAPDH",
            "PKM",
            "CS",
            "IDH3A",
            "NDUFS1",
            "SDHA",
            "UQCRC1",
            "COX4I1",
            "ATP5F1A",
            "MT-CO1",
            "MT-ND1",
            "MT-ATP6",
        ),
        "~90 OXPHOS subunits (13 encoded by mtDNA), ~1,100 mitochondrial proteins",
        "GO:0006096, GO:0006099, GO:0006119",
        go_terms=("GO:0006096", "GO:0006099", "GO:0006119"),
    )
)
_add(
    Library(
        "core.cytoskeleton",
        "core",
        "shape, movement, internal transport",
        ("ACTB", "TUBA1A", "TUBB", "MYH9", "KRT5", "KRT14", "VIM", "LMNA", "DYNC1H1", "KIF5B"),
        "~500 genes",
        "GO:0005856 cytoskeleton",
        go_terms=("GO:0005856",),
    )
)
_add(
    Library(
        "core.membrane_transport",
        "core",
        "move ions and molecules across membranes",
        ("ATP1A1", "SLC2A1", "SLC6A4", "ABCB1", "CFTR", "AQP1", "KCNQ1", "SCN5A", "CACNA1C"),
        "~400 SLC, ~50 ABC, ~400 ion-channel genes",
        "GO:0055085 transmembrane transport",
        go_terms=("GO:0055085",),
    )
)
_add(
    Library(
        "core.chromatin_epigenome",
        "core",
        "mark and pack DNA: the read/write layer above sequence",
        (
            "DNMT1",
            "DNMT3A",
            "DNMT3B",
            "TET1",
            "TET2",
            "EZH2",
            "KDM6A",
            "HDAC1",
            "KAT2A",
            "H3-3A",
            "SMARCA4",
            "CTCF",
        ),
        "~700 chromatin genes",
        "GO:0006325 chromatin organization",
        "the machinery that epigenetic clocks read",
        go_terms=(
            "GO:0006325",
            "GO:0032259",
        ),
    )
)

# ------------------------------------------------------------ blueprint ---
_add(
    Library(
        "blueprint.pluripotency",
        "blueprint",
        "the uncommitted starting state",
        ("POU5F1", "SOX2", "NANOG", "KLF4", "LIN28A", "DPPA3"),
        "~30 core factors",
        "Takahashi & Yamanaka 2006; Boyer et al. 2005",
        go_terms=(
            "GO:0019827",
            "GO:0043045",
        ),
    )
)
_add(
    Library(
        "blueprint.germ_layers",
        "blueprint",
        "gastrulation: ectoderm, mesoderm, endoderm",
        ("TBXT", "MIXL1", "EOMES", "MESP1", "SOX17", "FOXA2", "GATA6", "SOX1", "PAX6", "NODAL"),
        "~50 genes",
        "Uberon/CL developmental annotations",
        go_terms=(
            "GO:0007369",
            "GO:0001704",
            "GO:0001706",
            "GO:0001707",
            "GO:0001705",
            "GO:0007417",
        ),
    )
)
_add(
    Library(
        "blueprint.axes",
        "blueprint",
        "anterior-posterior, dorsal-ventral, left-right coordinates",
        (
            "HOXA1",
            "HOXA13",
            "HOXB4",
            "HOXC10",
            "HOXD13",
            "CDX2",
            "BMP4",
            "NOG",
            "CHRD",
            "NODAL",
            "LEFTY1",
            "PITX2",
            "SHH",
            "ZIC3",
        ),
        "39 HOX genes in 4 clusters; ~40 axis genes",
        "HOX literature; Uberon",
        "HOX clusters are read in genomic order along the body: literally positional code",
        go_terms=(
            "GO:0009952",
            "GO:0009953",
            "GO:0007368",
            "GO:0048598",
        ),
    )
)
_add(
    Library(
        "blueprint.segmentation",
        "blueprint",
        "somites: the repeated units of the trunk",
        ("HES7", "LFNG", "DLL1", "NOTCH1", "MESP2", "TBX6", "RIPPLY2", "FGF8", "WNT3A"),
        "~30 genes",
        "Pourquié lab; Matsuda et al. 2020 (human ~5 h period)",
        "an oscillator that lays down structure: blueprint and timer in one",
        go_terms=(
            "GO:0001756",
            "GO:0007389",
            "GO:0035282",
        ),
    )
)
_add(
    Library(
        "blueprint.signalling_toolkit",
        "blueprint",
        "the ~12 pathways every organ reuses: the genome's shared API",
        (
            "WNT3A",
            "CTNNB1",
            "SHH",
            "PTCH1",
            "GLI1",
            "BMP4",
            "TGFB1",
            "SMAD4",
            "FGF8",
            "FGFR1",
            "NOTCH1",
            "DLL1",
            "EGFR",
            "KRAS",
            "MAPK1",
            "PIK3CA",
            "MTOR",
            "JAK2",
            "STAT3",
            "YAP1",
            "RARA",
            "ALDH1A2",
            "EPHB2",
            "SEMA3A",
        ),
        "WNT, Hedgehog, BMP/TGF-β, FGF, Notch, RTK/RAS, PI3K/AKT/mTOR, JAK/STAT, Hippo, "
        "retinoic acid, Ephrin, Semaphorin",
        "Reactome signalling pathways",
        "same library, different context (when:) gives different organs",
        go_terms=(
            "GO:0016055",
            "GO:0007224",
            "GO:0030509",
            "GO:0008543",
            "GO:0007219",
            "GO:0007265",
            "GO:0007259",
            "GO:0035329",
            "GO:0048384",
            "GO:0043491",
            "GO:0007169",
            "GO:0007179",
            "GO:0048013",
            "GO:0071526",
            "GO:0000165",
            "GO:0031929",
            "GO:0042573",
            "GO:0045879",
        ),
        reactome=(
            "Signaling by WNT",
            "Signaling by Hedgehog",
            "Signaling by NOTCH",
            "Signaling by BMP",
            "Signaling by FGFR",
            "Signaling by TGF-beta Receptor Complex",
            "MAPK family signaling cascades",
            "PIP3 activates AKT signaling",
            "Signaling by Interleukins",
            "Signaling by Hippo",
            "Signaling by Retinoic Acid",
        ),
    )
)
_add(
    Library(
        "blueprint.organ_heart",
        "blueprint",
        "cardiac specification and morphogenesis",
        ("NKX2-5", "GATA4", "TBX5", "HAND1", "HAND2", "MEF2C", "ISL1", "TBX1", "MYH6", "MYH7"),
        "~100 genes",
        "Reactome; OMIM congenital heart disease",
        go_terms=("GO:0007507",),
    )
)
_add(
    Library(
        "blueprint.organ_nervous",
        "blueprint",
        "neural induction, neurogenesis, gliogenesis",
        ("SOX1", "SOX2", "PAX6", "NEUROG2", "ASCL1", "NEUROD1", "OLIG2", "OTX2", "EN1", "FOXG1"),
        "~500 genes",
        "Uberon/CL nervous system",
        go_terms=("GO:0007399",),
    )
)
_add(
    Library(
        "blueprint.organ_eye",
        "blueprint",
        "eye field and lens",
        ("PAX6", "SIX3", "RAX", "OTX2", "SOX2", "CRYAA", "VSX2", "MITF"),
        "~50 genes",
        "PAX6 is the eye master gene across animals",
        go_terms=(
            "GO:0001654",
            "GO:0030900",
            "GO:0060041",
            "GO:0048066",
        ),
    )
)
_add(
    Library(
        "blueprint.organ_endoderm_gut",
        "blueprint",
        "liver, pancreas, lung, thyroid from the gut tube",
        ("FOXA1", "FOXA2", "HNF1A", "HNF4A", "PDX1", "PTF1A", "NEUROG3", "NKX2-1", "SOX9", "CDX2"),
        "~100 genes",
        "Uberon endoderm derivatives",
        go_terms=(
            "GO:0001889",
            "GO:0031016",
            "GO:0030324",
            "GO:0030878",
            "GO:0048565",
            "GO:0031018",
            "GO:0001714",
        ),
    )
)
_add(
    Library(
        "blueprint.organ_kidney_urogenital",
        "blueprint",
        "kidney and gonads",
        ("PAX2", "PAX8", "WT1", "SIX2", "GDNF", "RET", "SRY", "SOX9", "FOXL2", "NR5A1"),
        "~80 genes",
        "OMIM; SRY/SOX9 vs FOXL2 sex determination switch",
        go_terms=(
            "GO:0001822",
            "GO:0008406",
            "GO:0007530",
        ),
    )
)
_add(
    Library(
        "blueprint.organ_musculoskeletal",
        "blueprint",
        "muscle, cartilage, bone, limbs",
        (
            "PAX3",
            "PAX7",
            "MYOD1",
            "MYF5",
            "MYOG",
            "SOX9",
            "RUNX2",
            "SP7",
            "TBX4",
            "TBX5",
            "SHH",
            "FGF8",
            "HAND2",
            "PITX1",
        ),
        "~150 genes",
        "limb: TBX5 forelimb, TBX4/PITX1 hindlimb, SHH from ZPA, FGF8 from AER",
        go_terms=(
            "GO:0007517",
            "GO:0001501",
            "GO:0060173",
            "GO:0001649",
            "GO:0048641",
            "GO:0035914",
        ),
    )
)
_add(
    Library(
        "blueprint.organ_blood_immune",
        "blueprint",
        "haematopoiesis and lymphocyte generation",
        ("RUNX1", "GATA1", "GATA2", "TAL1", "SPI1", "CEBPA", "PAX5", "RAG1", "RAG2", "IKZF1"),
        "~200 genes",
        "RAG1/2 rearrange V(D)J: the immune system compiles new genes at runtime",
        go_terms=("GO:0030097", "GO:0033151"),
    )
)
_add(
    Library(
        "blueprint.organ_skin",
        "blueprint",
        "epidermis and appendages",
        ("TP63", "KRT5", "KRT14", "KRT1", "KRT10", "IVL", "FLG", "LEF1", "EDAR", "SOX9"),
        "~100 genes",
        "TP63 is the epidermal master regulator",
        go_terms=(
            "GO:0008544",
            "GO:0001942",
        ),
    )
)

# ---------------------------------------------------------------- timer ---
_add(
    Library(
        "timer.telomere",
        "timer",
        "division counter: ends shorten unless telomerase runs",
        ("TERT", "TERC", "DKC1", "TERF1", "TERF2", "POT1", "TINF2", "RTEL1"),
        "~20 genes",
        "Harley 1990; Hayflick 1961; shelterin literature",
        "TERT is off in most somatic cells: the counter is deliberately not reset",
        go_terms=("GO:0000723",),
        reactome=("Telomere Maintenance",),
        noncoding=("TERC",),
    )
)
_add(
    Library(
        "timer.circadian",
        "timer",
        "24-hour clock in nearly every cell",
        ("CLOCK", "BMAL1", "PER1", "PER2", "PER3", "CRY1", "CRY2", "NR1D1", "RORA"),
        "~20 core genes drive ~40% of genes rhythmically",
        "Takahashi 2017 review",
        go_terms=("GO:0007623",),
        reactome=("Circadian clock",),
    )
)
_add(
    Library(
        "timer.segmentation_clock",
        "timer",
        "oscillator that times somite formation",
        ("HES7", "LFNG", "DLL1", "NOTCH1"),
        "shared with blueprint.segmentation",
        "Matsuda et al. 2020, Science",
        go_terms=(
            "GO:0001756",
            "GO:0007389",
            "GO:0035282",
        ),
    )
)
_add(
    Library(
        "timer.developmental_timing",
        "timer",
        "heterochronic control: when stages happen",
        ("LIN28A", "LIN28B", "MIRLET7A1", "IGF2", "H19", "MKRN3"),
        "~20 genes",
        "Ambros/Ruvkun lin-28/let-7; MKRN3 puberty brake",
        "let-7 accumulation is a developmental hourglass conserved from worms to humans",
        go_terms=(
            "GO:0040034",
            "GO:2000631",
            "GO:0048639",
        ),
        noncoding=("MIRLET7A1", "H19"),
    )
)
_add(
    Library(
        "timer.puberty_reproduction",
        "timer",
        "the second construction phase and its end",
        (
            "KISS1",
            "KISS1R",
            "GNRH1",
            "GNRHR",
            "TAC3",
            "MKRN3",
            "FSHR",
            "LHCGR",
            "ESR1",
            "AR",
            "FOXL2",
            "BMP15",
        ),
        "~50 genes",
        "OMIM; growth-plate closure via ESR1; ovarian reserve depletion",
        "menopause is a depletion timer: a fixed oocyte pool, no renewal",
        go_terms=(
            "GO:0032276",
            "GO:0001541",
            "GO:0060009",
            "GO:0022414",
            "GO:0030518",
            "GO:0042698",
        ),
    )
)
_add(
    Library(
        "timer.growth_longevity_axis",
        "timer",
        "growth-hormone / IGF-1 / mTOR / sirtuin axis",
        ("GH1", "GHR", "IGF1", "IGF1R", "INSR", "FOXO3", "SIRT1", "SIRT6", "MTOR", "PRKAA1", "KL"),
        "~40 genes",
        "Kenyon 2010; FOXO3 longevity GWAS; Klotho",
        "the same axis that drives construction accelerates ageing when left on",
        go_terms=(
            "GO:0008286",
            "GO:0048009",
            "GO:0031929",
            "GO:0008340",
            "GO:0031667",
            "GO:0060396",
        ),
    )
)
_add(
    Library(
        "timer.senescence_checkpoint",
        "timer",
        "the brake: stop dividing when damaged or old",
        ("CDKN2A", "CDKN1A", "TP53", "RB1", "SERPINE1", "IL6", "CXCL8"),
        "~30 genes plus the SASP secretome",
        "Campisi; López-Otín hallmarks of aging 2023",
        "p16INK4a expression rises with age in most tissues: a readable age register",
        go_terms=(
            "GO:0090398",
            "GO:2000772",
            "GO:0090399",
        ),
        reactome=(
            "Cellular Senescence",
            "Senescence-Associated Secretory Phenotype (SASP)",
        ),
    )
)
_add(
    Library(
        "timer.progeroid_maintenance",
        "timer",
        "genes whose failure accelerates ageing",
        ("LMNA", "WRN", "BLM", "ERCC6", "ERCC8", "ATM", "TERT", "DKC1", "ZMPSTE24"),
        "~15 genes",
        "OMIM progeroid syndromes",
        "the fastest way to find maintenance modules: see what breaks when each is removed",
        go_terms=(
            "GO:0006281",
            "GO:0000723",
            "GO:0006998",
        ),
    )
)
_add(
    Library(
        "timer.stem_cell_niches",
        "timer",
        "renewal capacity: the maintenance workforce",
        ("LGR5", "KIT", "CD34", "PROM1", "PAX7", "SOX2", "NES", "GFAP", "KRT15", "BMI1"),
        "~50 markers",
        "Cell Ontology stem cell terms",
        go_terms=(
            "GO:0019827",
            "GO:0072089",
            "GO:0048863",
            "GO:0017145",
        ),
    )
)

# -------------------------------------------------------------- systems ---
_add(
    Library(
        "systems.endocrine",
        "systems",
        "hormones: slow, broadcast messages",
        (
            "INS",
            "GCG",
            "GH1",
            "PRL",
            "TSHB",
            "TG",
            "TPO",
            "POMC",
            "CRH",
            "CYP11A1",
            "CYP19A1",
            "LEP",
            "ADIPOQ",
            "EPO",
        ),
        "~200 hormone and receptor genes",
        "Reactome: hormone biosynthesis and signalling",
        go_terms=(
            "GO:0005179",
            "GO:0009755",
            "GO:0042445",
            "GO:0010817",
        ),
        reactome=("Peptide hormone metabolism", "Hormone ligand-binding receptors"),
    )
)
_add(
    Library(
        "systems.nervous",
        "systems",
        "fast, addressed messages",
        ("SCN1A", "KCNA1", "CACNA1A", "GRIN1", "GABRA1", "SLC6A3", "TH", "CHAT", "SYT1", "SNAP25"),
        "~1,000 synaptic and channel genes",
        "SynGO; GO:0007268 chemical synaptic transmission",
        go_terms=("GO:0007268", "GO:0001508"),
    )
)
_add(
    Library(
        "systems.immune",
        "systems",
        "surveillance, defence, and tissue repair signals",
        ("HLA-A", "HLA-DRB1", "B2M", "TLR4", "IL1B", "IL6", "TNF", "IFNG", "CD4", "CD8A", "TRAC", "IGHM"),
        "~1,500 genes; HLA is the most polymorphic region",
        "ImmPort; GO:0006955",
        go_terms=("GO:0006955",),
        reactome=("Immune System",),
    )
)
_add(
    Library(
        "systems.circulation_respiration",
        "systems",
        "transport of oxygen and nutrients",
        ("HBB", "HBA1", "VEGFA", "KDR", "NOS3", "ACE", "AGT", "NPPA", "SFTPC", "TTN"),
        "~300 genes",
        "Reactome; GO:0001525 angiogenesis",
        go_terms=("GO:0001525", "GO:0015671", "GO:0008015", "GO:0007585", "GO:0060047"),
    )
)
_add(
    Library(
        "systems.extracellular_matrix_adhesion",
        "systems",
        "the structural fabric between cells",
        ("COL1A1", "COL2A1", "COL4A1", "FN1", "LAMA1", "ELN", "CDH1", "CDH2", "ITGB1", "ITGA5", "GJA1"),
        "~300 matrisome-core genes, ~1,000 associated",
        "Naba matrisome; GO:0031012",
        go_terms=("GO:0031012", "GO:0007155", "GO:0005921"),
    )
)
_add(
    Library(
        "systems.ligand_receptor_protocol",
        "systems",
        "the message format between any two cells",
        ("WNT5A", "FZD7", "TGFB1", "TGFBR2", "EGF", "EGFR", "CXCL12", "CXCR4", "DLL4", "NOTCH1"),
        "~2,000 curated ligand-receptor pairs",
        "CellPhoneDB; CellChat",
        "this is the dependency graph between tissues, computable from expression",
        go_terms=("GO:0048018", "GO:0038023"),
    )
)
_add(
    Library(
        "systems.metabolic_homeostasis",
        "systems",
        "liver, fat, muscle, gut sharing fuel",
        ("PPARA", "PPARG", "SREBF1", "PCK1", "G6PC1", "LDLR", "APOB", "APOE", "SLC2A2", "SLC2A4"),
        "~500 genes",
        "Reactome metabolism",
        go_terms=(
            "GO:0042593",
            "GO:0055088",
            "GO:0019216",
            "GO:1904659",
        ),
    )
)

# ---------------------------------------------------------------- parts ---
_add(
    Library(
        "parts.cell_types",
        "parts",
        "inventory of cell types",
        (),
        "Cell Ontology ~2,900 classes; Human Cell Atlas finds many more states",
        "http://purl.obolibrary.org/obo/cl.obo",
        "GenomeOS imports these as CellType entities in Phase 1",
    )
)
_add(
    Library(
        "parts.tissues_organs",
        "parts",
        "inventory of anatomical parts",
        (),
        "Uberon ~15,000 terms; ~78 organs; 11 organ systems",
        "http://purl.obolibrary.org/obo/uberon/uberon-ext.obo",
    )
)
_add(
    Library(
        "parts.genome_census",
        "parts",
        "what the 3.1 Gb actually contains",
        (),
        "~19,400 protein-coding genes (1.5% of bases), ~20,000 lncRNA genes, ~14,000 pseudogenes, "
        "~500 tRNA, ~1,600 transcription factors, 518 kinases, ~800 GPCRs (~400 olfactory), "
        "~1 million candidate cis-regulatory elements, ~45% transposable elements, "
        "mtDNA: 37 genes (13 proteins, 22 tRNA, 2 rRNA)",
        "GENCODE 50; ENCODE SCREEN; Lambert et al. 2018 (TFs); Manning 2002 (kinome)",
        "most of the genome is regulation, repeats and history, not protein code",
    )
)


def by_layer(layer: str) -> list[Library]:
    return [l for l in LIBRARIES.values() if l.layer == layer]
