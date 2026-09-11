# Action plan progress

Updated as each task is completed (tests + lint green before moving on).

| Task | Status | Evidence |
|---|---|---|
| 0.0 tooling: ruff lint + format, progress file | done | `uv run ruff check .` clean, 30 tests |
| 1.1 GFF3 → BioIR (GENCODE 50) | done | `genomeos annotate`, `genomeos gene`; chr21+chrM parse in 8 s; 2,703 complete coding transcripts translate with zero internal stops; 2 of 2,705 GENCODE models start off-ATG (NCAM2-201, DSCAM-203) |
| 1.2 random-access FASTA | done | `genomeos index`; samtools-compatible `.fai`; 10,000 random fetches on chr21 in <1 s, identical to in-memory |
| 1.3 VCF → diploid HG002 chr21 | done | `genomeos twin build`; 55,210 PASS variants, 0 reference mismatches; GIAB genotypes are unphased so het ALT goes to hap2 by policy |
| 1.4 alternative start codons | done | AUU/AUA/AUC/GUG initiators for mtDNA; all 13 mitochondrial proteins at published lengths (MT-ND2 347 aa starts AUU) |
| 1.5 GO/Reactome → library membership | done | `genomeos libs --data/--verify/--gene`; Reactome membership propagated up the hierarchy; catalogue agreement 95.5% (20 of ~440 genes are developmental TFs / stem-cell markers that human GO annotates only as "regulation of transcription"); two of my original GO ids were obsolete and were caught by the data |
| 1.6 variant effect on CDS | done | `genomeos variant`; ClinVar chr21: >90% agreement for nonsense, missense, synonymous, frameshift (ClinVar calls truncating frameshifts "nonsense"); APP chr21:25897620 C>T → A673T, the known protective variant |
| 2.2 CellType from Cell Ontology | done (module) | `CellTypes` loads 2,900+ CL terms, lineage, search, compiles to BioIR entities |
| 2.3 SBML engine (BioModels) | done (in-house) | `genomeos/runtime/sbml.py`: MathML evaluator + RK4; BIOMD0000000012 (Elowitz repressilator) parses, assignment rules evaluate to published values, proteins oscillate; libRoadRunner has no Python 3.14 wheel yet, so the adapter is deferred and the engine keeps the same `run()` contract |
| 2.4 Boolean network engine | done (in-house) | `genomeos/runtime/boolean.py`: BoolNet files, sync/async update, attractors, export to BioIR rules; Fauré 2006 mammalian cell cycle gives the G1-arrest fixed point without CycD and a cyclic attractor with it |
| 2.5 process-bigraph spike | go | `genomeos/runtime/compose.py` (optional extra `compose`): NetworkProcess + AgeingProcess run in one Composite on a shared clock; registration is `core.register_link`, updates are additive deltas; decision: adopt as the scheduler for Phase 4 composition |
| 2.6 ligand-receptor protocol | done | `genomeos/knowledge/ligand_receptor.py`: CellPhoneDB → 1,500+ `binds` rules with classification context; CXCL12–CXCR4, EGF–EGFR, DLL4–NOTCH1, TGFB1–TGFBR2 present |
| 2.1 BioLang v0.2 | done | nested `transcript` blocks, `cell_type` (expresses → context silencing), `event` (rate, when, effects), `import` with `bio.std.*` prelude (`cell_types`, `ageing`) and relative paths; v0.1 files unchanged; `genomeos check --context cell_type=X` reports active rules and silenced genes |
| 2.2 CellType CLI | done | `genomeos cells --search/--lineage`; 2,900+ Cell Ontology types compile to BioIR |
| 3.2 telomere length from reads | done (synthetic validation) | `genomeos/genome/telomere.py`: TelSeq-style estimator over FASTQ or BAM (stdlib BGZF decoding); recovers a known 8 kb telomere from synthetic reads; a real 30x HG002 BAM (~100 GB) is not downloaded, so real-data validation is pending |
| 3.3 BioTwin fork / run / diff | done | `genomeos/twin/twin.py`: JSON twins with measured state (age, telomere, epigenetic age) and environment; runs start from the measured state; LoF variants in timer/maintenance library genes modify parameters with `inferred` evidence and lowered confidence; TERT LoF → faster attrition, TP53 LoF → less senescence |
| 3.4 AlphaGenome adapter | done (mock-tested) | `genomeos/predict/alphagenome_adapter.py` (optional extra `predict`, alphagenome 0.9.0 installed); emits `predicted` rules capped at confidence 0.7; live call needs ALPHAGENOME_API_KEY (not set) |
| 3.5 uncertainty per level | done | `genomeos/runtime/uncertainty.py`: molecular / cellular / tissue / organism with UNKNOWN distinct from low; printed by `genomeos run` and `genomeos age`, attached to twin runs |
| 3.1 epigenetic age | done, validated on real data | `genomeos clock`; Horvath 2013 (353 CpGs, intercept 0.6955) and Hannum 2013 (71 CpGs) from published coefficients, no dependencies (biolearn needs torch, unavailable on 3.14); GEO GSE41169 whole blood: Horvath r > 0.8, MAE < 10 y (e.g. 65 → 65.4, 32 → 30.0) |
| 4.1 spatial engine | done (in-house) | `genomeos/runtime/spatial.py`: 2-D diffusion/decay fields, sources, cells with fate rules; Wolpert French flag forms blue/white/red bands in order from one source |
| UI: Twin view + uncertainty | done | web UI gains a Twin tab (create, fork with variants/environment, run A vs B, diff, uncertainty) and uncertainty tables on Program and Ageing |
| 4.2 segmentation clock | done | `data/demo/segmentation_clock.bio` (HES7 loop, FGF wavefront, period 5 h from Matsuda 2020) + `runtime/segmentation.py` clock-and-wavefront; segment count = frozen cells / (speed × period) within one segment; halving the clock rate halves the count |
| 4.3 gastrulation | done | `data/demo/gastrulation.bio` tristable SOX2/TBXT/SOX17 switch read against a clamped NODAL gradient; endoderm → mesoderm → ectoderm in order; proportions within the stated (confidence 0.3) expectations; the module reports itself as low confidence because its mutual-repression rules are inferred |
| 4.4 debugger | done | `runtime/debugger.py`: breakpoints (`TetR > 50`, `divisions >= 10`), step, `explain()` lists the terms driving a species with rule evidence and confidence, ageing events traced to their parameters |
| 5a second organism | done | Ensembl GFF3 dialect supported; C. elegans WBcel235 chromosome III: >2,000 protein-coding genes, >97% of coding transcripts translate cleanly, `lin-12` found |
| 5b minimal organism | done (whole organism) | `genomeos grow data/organisms/celegans/embryo.bio --until 6000 --compare`: one zygote → 2,183 cells → the 959-cell adult hermaphrodite plus Z2/Z3; every cell by name, every terminal fate (961/961) and every programmed death (131/131) as in Sulston 1983 / Sulston & Horvitz 1977; division timing on distilled per-lineage timers median 12 min from the reference; founders decided by mechanism (PAR, SKN-1, PIE-1, PAL-1, POP-1, Wnt from P2, Notch from P2 and MS); the earlier `genomeos organism` view now runs on the same Body runtime |
| 5c BioForge | done | `genomeos forge`: random-restart log-space search over gene/rule/parameter knobs with hard constraints; retunes the repressilator period to a target within 15% while keeping oscillation; all changed values become `predicted: BioForge search` at confidence ≤ 0.3 |
| UI: Space + Debugger | done | web UI gains Space (French flag, segmentation, gastrulation) and Debugger (breakpoints, state, explanation with evidence, event trace) |
| storage: stream, distil, discard | done | `genomeos data status/distil/clean/results`; HG002 telomere estimated from 2 M reads streamed over HTTP in 12 s with zero disk use (2,948 bp TelSeq scale ≈ 6 kb Southern-equivalent); GEO clocks, ClinVar agreement, library membership and HG002 chr21 haplotypes distilled to small JSON summaries under data/results/, committed; raw inputs (~700 MB) deleted; tests read the summaries when raw files are absent |
| embed what was learned | done | clock coefficients and the distilled library membership (114 KB) ship inside the package; GENCODE distilled to a chr21+chrM subset (1.9 MB) and the 153 MB file deleted; data folder now ~130 MB; docs/LESSONS.md gathers the lessons the core carries |
| twin from measured state (track 1) | done | HG002 telomere streamed (TelSeq 2,948 bp ≈ 5.9 kb Southern); clock CpG coordinates distilled from the 450k manifest into the package (418 probes, 13 KB; manifest discarded); HG002 nanopore bedMethyl (hap1+hap2, 1.2 GB, 53 M rows) streamed in 54 s: 409/418 clock CpGs at 36x, Horvath 97.5 y / Hannum 62.7 y (lymphoblastoid cell line signature, recorded as such); HG002 twin now carries measured telomere and epigenetic age; `genomeos methylation`, Twin view "apply measured state" |
| rules from data (track 3) | done (first rule) | `genomeos calibrate`: fits divisions_per_year so the model reproduces a measured telomere at a given age; HG002 at ~45 y → 3.8 divisions/yr, net 91 bp/yr (inferred, conf 0.4); twin runs use saved calibrations |
| sequence grammar (BioLang from the genome's own delimiters) | started | docs/SEQUENCE-GRAMMAR.md; `genomeos signals learn|scan`: PWMs for donor, acceptor and Kozak learned from chr21 (GT 99.1%, AG 99.8%, ATG 220/221, TGA/TAA/TAG 220/221, polyA signal 66%, CpG-island promoters 58%); separability measured: 90% donor recall at ~7 background hits per kb, so segments must be parsed by grammar and context, never asserted |
| UI: friendliness | done | Start tab with what each area does and a five-step walk-through; icons on every tab, card and button; 60 hover explanations (data-tip) on controls and "?" help markers on every panel; evidence-pill legend |
| genome anatomy | done | `genomeos anatomy` and the Anatomy tab: composition (cds / other exonic / intron / intergenic / gap), gene counts by type and strand, exon/intron/gene lengths, single-exon fraction, transcripts per gene, density, spacing, neighbour orientation, overlaps, nesting, CpG islands, homopolymers, microsatellites, telomeric blocks; mtDNA / human chr21 / worm III compared in docs/GENOME-ANATOMY.md with the design budgets a from-zero genome must choose between |
| UI: Progress tab | done | jobs run as background processes with live progress bars and logs (start from the UI; whitelisted only), the whole-genome inventory table fills in as chromosomes finish, the project task log is rendered from this file, and a busy bar shows every request in flight |
| whole-genome anatomy | done | `scripts/anatomy_genome_wide.py` streams each human chromosome (download, count, delete); chr1: 2,066 coding genes, CDS 1.5%, intron 62%, 109 s; all 25 chromosomes counted in ~35 min: 3.088 Gb, 20,107 coding genes, CDS 1.21%, intron 61.3%, intergenic 27.9%, 39,410 CpG islands; peak disk one chromosome; `anatomy_hg38_by_chromosome.json` |
| UI: block viewer / editor | done (v1) | Blocks tab: canvas map of a chromosome's blocks to scale and nested (genes → transcripts → exons/CDS), UNKNOWN spans between genes, CpG islands, repeats, telomeric blocks, gaps; zoom, pan, minimap, search-and-jump, type filters that mute the rest, selection with details; edit mode drags genes/transcripts/exons with server plausibility checks (inside parent, no same-strand overlap, reading frame) and exports the proposal as predicted/0.3 |
| block editor summary + deep links | done | percentage coverage per block type (union of intervals) above the map; `#blocks?locus=chr21:a-b` deep links; ruler spacing |
| cancer section | started | `genomeos cancer distil` (cBioPortal REST, 110-gene driver panel over MSK-IMPACT 2017: 10,945 tumours, 58 types, 30 s, one 100 KB summary kept), `genes`, `gene`, `types`, `compare --normal --tumour` (somatic set → consequences → driver/hotspot grading → suggested cancer types → surface targets → AI-agent packet); docs/CANCER.md; MCP evaluated and not used (needs auth, not reproducible) |
| minimal cell design | done (estimate) | `genomeos design CELLTYPE`: essential (Hart 2017 CEGv2, 684) + identity libraries + master switches; neuron ≈ 2,700 genes, 3.5 Mb compact / 13 Mb dense vs 3.1 Gb; docs/DESIGN-MINIMAL-CELL.md answers the neuron question |
| UNKNOWN investigator | done (v1, chr21: 94% of UNKNOWN bases classified in 55 s; centromere found from CENP-B boxes; Alu density by mismatch-tolerant core search) | `genomeos unknown`: classifies intergenic blocks largest first from sequence alone (gap, telomere, centromere, tandem, low complexity, interspersed repeat via 16-mers learned from the chromosome, Alu core, coding remnant, promoter-like, gene desert); pattern file for named motifs; duplication candidates by shared 20-mers; Blocks tab labels and colours the classes; docs/UNKNOWN.md |
| ENCODE regulatory elements | done | cCRE registry streamed (64 MB) and distilled per chromosome; chr21: 13,356 elements; `regulatory` class in the investigator (96.3% of UNKNOWN classified), ENCODE blocks in the map, CTCF boundary motif in the pattern file; genome-wide UNKNOWN investigation running as a background job |
| UI: class summary + highlight | done | chips with share of bases per block type and UNKNOWN class; click a chip or "Highlight all" in the block panel to mute everything else |
| conceptual model: nodes / reader / writer / executor | recorded | docs/NODES-READER-WRITER.md maps it to TADs and CTCF (nodes), the epigenome (reader), replication/epigenetic/germline writers, runtime (executor); next design step: a `Domain` block above `Gene` |
| molecules: RNA and proteins | started | `genomeos protein SYMBOL` and the Molecules tab: our translation of the canonical transcript checked against UniProt (APP: 770 aa, 100% identical to P05067), UniProt function, location and features (curated), AlphaFold model with per-residue pLDDT (predicted) drawn as a rotating backbone coloured by confidence; structures cached, nothing else stored |

## 2026-09-10 (evening 2) — flow trace, nodes, protein compiler

- **Flow tab and `genomeos/flow`**: `CentralDogmaTrace` maps base ↔ mRNA ↔ codon ↔ residue both ways and traces a single-base change (HGVS c./p.). APP A673T reproduced from the plus-strand coordinate. Tests: `tests/test_flow_trace.py` (plus and minus strand).
- **Nodes above genes**: `genomeos domains` infers domains between CTCF-only ENCODE elements (inferred, 0.4). chr21: 228 domains, median 115 kb, max 5.1 Mb; 137 without coding genes; the RUNX1 domain holds 7 coding genes and 237 enhancers. Domains are a lane in the block map. Result `data/results/domains_chr21.json`.
- **Protein compiler** (`genomeos protein SYMBOL --compile`, Molecules tab): Ensembl + UniProt + InterPro + PDB (with method) + AlphaFold + Reactome + STRING + Human Protein Atlas → one definition keyed by UniProt accession; evidence and confidence per section; ProteinState separate from the definition. TP53: 40 transcripts → 37 protein products → 9 isoforms; 311 experimental structures (266 X-ray, 30 NMR, 15 cryo-EM); 46 pathways; 374 associations, 200 with experimental support. Transient source failures are refetched, never cached for good. Docs: docs/PROTEIN.md.
- **Knowledge to execution**: `genomeos pathway ID --knockout X` runs a Reactome pathway export as a reachability graph (no rate laws in the export, so no integration); a knockout removes every entity containing the protein. Stabilization of p53: 15 reactions, 13 reachable, 6 lost without TP53 (MDM2 binding, CHEK2/ATM phosphorylation, ubiquitination, translocation). Molecules tab: click a pathway chip to run it. Tests on a toy graph and on the committed export `data/models/reactome_R-HSA-69541.sbml`.
- **BioLang importer**: `genomeos protein X --bio` emits a `protein` block with accession, sequence, isoforms, domains, pathways, interactions (physical channel only), structures; the parser accepts these properties.
- **RegulatoryElement in BioIR** with targets: promoters reach the TSS within 1 kb (0.8), enhancers reach the coding genes of their CTCF domain nearest first (0.4 / 0.25), insulators reach nothing; BioLang `element` block; `genomeos regulation --gene G`; regulation layer drawn in the Flow tab and printed by `genomeos flow`. APP: 1 promoter element, 98 reachable enhancers (89 intragenic), 2 bounding insulators, alone in its 175 kb node.
- **The test human on any chromosome**: `data fetch --individual` streams the GIAB HG002 benchmark once and keeps per-chromosome variant files (data/reference, not committed); lookups report whether HG002 carries a variant anywhere, twins can be built for any fetched chromosome (chr22: 18,035 and 50,162 variants applied per haplotype, 0 reference mismatches; 21 chromosome files, 158 MB). Jobs registry: on-disk pid liveness (a restarted server no longer starts twin jobs; five duplicate proteome runs and one UNKNOWN run were found and stopped).
- **Selenocysteine**: the engine reads UGA as Sec on `seleno`-tagged transcripts (SELENOM 47 → 145 aa, three chr22 selenoproteins now exact); compiler prefers the longest primary-name UniProt entry (MIEF1 microprotein case); verification adds gapped similarity to separate boundary differences from real disagreements (chr22: 4 real).
- **Translation verified against UniProt** (`genomeos verify --chrom C`): chrM 13/13 identical; chr21 91.5% canonical identical, 99.1% with some isoform, 2 real disagreements. Caught and fixed a compiler bug (UniProt `gene_exact` matches synonyms: MIF, DMC1, RRP1 and 13 others had another gene's entry; primary-name match now, cache repaired). UNKNOWN job progress now counts chromosomes done with curated repeats (the bar read 24 of 25 during the re-run because the first pass's summary was counted).
- **Proteome chr22 compiled** (432 coding genes from the Ensembl gene list, 21 min): sequence 98.6%, function 86.6%, domains 98.6%, pathways 54.4%, interactions 84.0%, experimental structure 47.9%, expression 99.3%, disease 25.0%. The genome-wide proteome job continues smallest chromosome first.
- **Fetch from the UI**: one `fetch_<chrom>` job per chromosome (hidden until started); Blocks tab has a fetch control; chrY fetched through it in under 90 s.
- **Any chromosome on demand**: `genomeos data fetch --chrom chr22` (sequence, GENCODE rows streamed and kept, ENCODE elements, RepeatMasker); `default_gencode` finds single-chromosome model files; commands follow `--chrom` for the sequence. chr22 fetched in 4 min (17 MB kept): 247 domains, CHEK2 traced end to end, 79,521 curated repeat copies. CI lint fixed (a hand-formatted block in cli.py).
- **`bio` toolchain** (roadmap 1.x, BioLang as an engine): `bio check|compile|run|test|repl`; `# test:` lines inside a module are its own claims (module facts, species final/peaks/min/max, no NaN); a REPL that holds forward references until declared. docs/BIO-TOOLCHAIN.md.
- **Variant lookup through every layer**: `genomeos lookup chr:pos REF>ALT`, `/api/lookup`, Flow tab card. VEP consequence, ClinVar/dbSNP/COSMIC/gnomAD/PubMed, local trace agreement, UniProt features at the residue, structures, pathway loss for truncations. APP A673T resolves to the amyloid-beta start (residue 672) with ClinVar protective; TP53 R175H pathogenic with 65 citations.
- **Knowledge graph**: `genomeos graph [GENE]`, `/api/graph`, neighbourhood card in Molecules (force layout, drag, evidence-coded edges). chr21 + chrM: 4,175 nodes, 32,520 edges, 260 compiled proteins, 133 components. Mitochondrial proteome compiled: 13 of 13 proteins answer every question (`proteome_chrM`).
- **RNA layer**: `genomeos rna GENE` and an RNA card in Molecules: every isoform from the gene models with biotype, exon count, spliced and coding length and canonical/MANE tags (APP: 63 transcripts, 53 coding), plus GTEx median TPM in 54 tissues with a pattern summary (APP: expressed in all tissues, brain frontal cortex highest at 561 TPM). Flow tab links a traced variant to its residue on the structure.
- **Block map shows curated repeats** (RepeatMasker copies with class colours, family and divergence on hover) wherever a chromosome has been distilled with `genomeos repeats`; Progress tab gained the genome-wide UNKNOWN table and the proteome coverage card.
- **Curated repeats in the UNKNOWN classifier**: `genomeos repeats` distils RepeatMasker (UCSC API, 20 MB once, 0.7 MB kept) and the classifier grades blocks by curated coverage first. chr21: long_orf 2.27 → 0.13 Mb, LINE 3.16 Mb curated, classified 96.3% → 97.6%. The genome-wide script fetches repeats per chromosome and discards them.
- **Genome-wide UNKNOWN classification done**: 26,806 blocks, 1.009 Gb, 96.6% classified in 92 min with one chromosome on disk at a time; regulatory 33.6%, long_orf 23.6% (mostly LINE-1 ORF2, an L1 signature is the next pattern), gap 15.8%, unique intergenic 10.5%, centromere 4.8%. Per-chromosome results `unknown_chr*.json`, summary `unknown_genome_wide.json`, table in docs/UNKNOWN.md.
- **Cancer, tumour alone**: `genomeos cancer tumour --vcf T.vcf --packet p.json` and the Cancer tab. Ensembl VEP REST annotates any chromosome (consequence, HGVS, SIFT/PolyPhen, COSMIC, gnomAD; 200 variants per call, cached per variant); gnomAD ≥ 1% set aside as likely germline; cBioPortal grading; mutation burden with its assumption written in; mutant peptide windows from UniProt; truncated drivers run through their Reactome pathways (demo RUNX1 Y480*: 64 reactions lost over 10 pathways). Positions are 1-based at the VEP and reporting boundary (the internal Variant is 0-based; caught by the demo disagreeing with the local classifier).
- **Proteome chr21 compiled** (216 coding genes, three passes because transient Ensembl failures had to be refetched; the compiler now treats never-attempted sections as failed): origin 100%, sequence 98.6%, domains 96.3%, function 89.8%, interactions 82.4%, pathways 62.0%, experimental structure 45.4%, predicted structure 98.1%, expression 100%, disease 25.5%. Result `proteome_chr21`.

## 2026-09-11 — therapeutic targeting and the design dataset

`genomeos/therapeutics/`: from a tumour's alterations to what physically
distinguishes those cells, what a therapy could reach, and which mechanism the
biology supports. Research hypotheses with evidence per claim; docs/THERAPEUTICS.md.

- **Pipeline** (`pipeline.py`): consequence → localisation → accessibility →
  expression and selectivity → trafficking → neoantigen/HLA → pathway-induced →
  normal-tissue exclusion → structure and epitope → mechanism selection →
  ranking. Each stage is a function over explicit inputs and a provider, so any
  one is replaceable; the whole pipeline runs offline against the caches.
- **Accessibility is not a membrane word.** The correction that mattered most:
  UniProt annotates SRC as `Cell membrane`, and SRC sits on the inner leaflet on
  a myristoyl anchor. External accessibility is now claimed only with positive
  evidence of an outward-facing part (extracellular topological domain,
  transmembrane segment with a plasma-membrane location, GPI anchor, cell-surface
  annotation, or signal peptide). A lipid anchor with no transmembrane segment is
  reported as *cytoplasmic face*; a bare membrane annotation as *side not
  established*. Before the fix the demo called SRC, YAP1 and RARA surface targets;
  after it, only APP and CSF1R, which is correct.
- **One level deeper**: a surface protein whose *mutation* sits in the
  cytoplasmic tail has no mutation-specific extracellular epitope. The demo's
  APP N770K is exactly that (UniProt topology: extracellular 18–701,
  transmembrane 702–722, cytoplasmic 723–770), and the report says so.
- **Truncations do not make neoepitopes.** A premature stop leaves every
  remaining residue identical to wild type, so RUNX1 Y480* is reported as
  yielding no mutation-derived peptide, with nonsense-mediated decay flagged as a
  further reason the truncated protein may not exist. The missense APP N770K does
  yield peptides (4 mutation-spanning windows; the change is the last residue).
- **15 mechanisms** (`mechanisms.py`) with declarative gates and weighted
  factors: antibody blocking and agonism, ADCC, ADCP, complement, T-cell and
  NK-cell engagers, ADC, targeted radionuclide, immune-marker delivery, RNA
  delivery, tumour-suppressor restoration, genome editing, TCR and TCR-mimic. A
  failed gate scores zero and prints why; an unknown factor is dropped and listed
  as a blocking unknown, and compatibility is multiplied by input coverage so a
  mechanism cannot rank highly on ignorance. Target, binder, mechanism, cargo and
  effector are separate objects, so `payload = none` is a complete mechanism.
- **Scoring** (`scoring.py`): 13 component scores, each with the sentence that
  produced it, and one overall whose formula and adjustments are printed.
  Unknown dimensions are excluded, never defaulted; missing normal-tissue data
  caps the overall at 0.6; a target no mechanism can address caps design
  readiness at 0.15.
- **Evidence** (`evidence.py`): source, source type, claim, level
  (clinical / human / preclinical / in vitro / computational / inferred) and
  confidence, grouped by level in the report, never mixed. `Missing` records what
  could not be established and what would resolve it.
- **Design dataset** (`design.py`): `TherapeuticDesignDataset` v1.0 with
  provenance (GenomeOS version, pipeline version, provider versions, input
  hashes), per-target recognition specification, positive set, negative set,
  tumour-versus-normal differential, structure identifiers, binder requirements,
  desired action, encoding-strategy classification and design readiness. Five
  layers written side by side plus the text report. The negative set is what
  makes selectivity tractable: wild-type protein, wild-type and similar
  peptide/HLA complexes, every healthy tissue above threshold with its level, the
  closest human paralogues with identity (RUNX1 → RUNX3 65%, RUNX2 60%;
  APP → APLP2 50%, APLP1 37%; CSF1R → KIT 39%, FLT3 30%, PDGFRA/B 29%), and
  polymorphism flagged as not screened. No sequence, construct, vector, cassette,
  formulation or protocol is produced anywhere.
- **New sources**: Human Protein Atlas search API for per-tissue consensus nTPM
  across 20 tissues (CC BY-SA 4.0); Open Targets Platform GraphQL for antibody
  tractability, approved drugs and trials by modality, and curated safety
  liabilities (CC0); Ensembl Compara for human paralogues with sequence identity.
  All cached under `data/knowledge/therapeutics/`. No peptide/HLA predictor is
  wired in: `NoNeoantigenPredictor` reports binding unavailable rather than
  inventing affinities.
- **Patient data levels 1–9** are reported, so a conclusion is never read above
  its data. The demo runs at level 1.
- **Integration**: `genomeos therapeutic --tumour ... [--normal --rna --hla
  --purity --out --report --spec GENE --offline]`; `genomeos cancer tumour
  --therapeutic`; `POST /api/therapeutics`; a **Targets** tab in `genomeos serve`
  with candidate cards, score bars, a mechanism panel and expandable evidence.
- **Demo** (`data/demo/cancer_tumour.vcf`, 5 variants, level 1): APP
  direct surface (score 0.50, accessibility 1.00, normal-tissue safety 0.00 —
  669 nTPM in cerebral cortex, 465 in kidney, 356 in heart — so the safety cap
  binds); CSF1R reached indirectly from disrupted RUNX1 (0.50, best mechanism
  blocking antibody 61%); SOD1, RUNX1 and NRIP1 intracellular with every surface
  mechanism gated off and the reason printed. `APP AND CSF1R` lowers the worst
  weighted healthy-tissue load 12.5-fold, with the bulk-RNA caveat attached.
  Outputs committed under `data/demo/therapeutics/`.
- **Tests**: `tests/test_therapeutics.py`, 41 offline tests with stubbed
  providers: nuclear protein keeps the peptide route and loses the surface one,
  lipid-anchored protein is not a surface target, GPI protein is, mutation in the
  cytoplasmic tail is not a surface epitope, high essential-organ expression is
  penalised, missing expression stays unknown, internalisation raises ADC and
  lowers it when absent while ADCC survives, a mechanism cannot score well on
  missing inputs, no HLA limits the peptide route, combination logic computes its
  gain, the dataset keeps the DNA origin and contains no nucleotide run, readiness
  is gated, and the evidence graph links a mechanism back to the variant.
  Whole suite: 170 passed, 4 skipped.

## 2026-09-11 — from one cell to an organism (BioLang v0.3, Body runtime)

- **Research and plan**: docs/ORGANISM-FROM-ONE-CELL.md (state of the art, the data that exists, design, plan). Nobody has an end-to-end one-cell-to-organism engine; the parts exist; C. elegans is the only organism with a complete cell-by-cell ground truth.
- **Reference lineage distilled**: `genomeos data distil --only celegans_lineage` streams the complete timed lineage (WormWeb, CC BY 1.0; Sulston 1983 + Sulston & Horvitz 1977) into `data/results/celegans_lineage_cells.json` (2,183 cells, 131 deaths, 961 alive in the adult = 959 somatic + Z2/Z3) and generates three BioLang modules: `timers.bio` (75 cycle timers per founder lineage and generation, mean and sd), `lineage_embryo.bio` and `lineage_larva.bio` (every division with daughter names, every death, every terminal fate, all cited). Packer 2019 (GEO GSE126954, 86,024 cells) distilled to a lineage → cell-type table and checked against the lineage: 173 of 232 comparable lineage ids agree (74.6%); the disagreements are pharynx and glia internal nodes where the two sources label at different depths (`data/results/celegans_packer2019.json`).
- **BioLang v0.3** (docs/BIOLANG-v0.3.md): `organism` (species, genome, root, bootstrap cell type, maternal factors, environment, tempo, observe, assert, reference), `stage`, `timer` (duration, sd, lengthening, when), `signal` (contact / gradient / systemic; from, to, sets), `decision` (divide, differentiate, migrate, quiesce, die; when, daughters, lineages, asymmetric, timer, after, fraction) and `domain`; `when` clauses accept alternatives and comparisons; diamond imports merge once. BioIR 0.3 carries them and round-trips.
- **Body runtime** (`genomeos/runtime/body.py`): discrete-event growth from one cell; per cell read → decide → wait → write; first matching decision per action wins (mechanism written before lookup); signals re-decide a pending division; asymmetric factor inheritance; populations as cells with a count; UNKNOWN stops reported, never filled in; organism-level uncertainty populated for the first time.
- **Result**: the whole worm grows in 0.4 s: 2,183/2,183 cells by name, 0 parent mismatches, 961/961 fates, 131/131 deaths; timing median 12 min (embryo 8 min, p90 25 min) on distilled timers; the alive-cell curve is under the reference between 350 and 500 min because per-generation mean timers smooth the real spread (194 vs 274 at 350 min), which the `assert` in the program flags. `genomeos grow FILE --compare` prints the diff; `tests/test_body.py` holds the numbers.
- **What is mechanism and what is lookup**: the founder divisions, the germline (PIE-1), MS/E (Wnt from P2, POP-1 asymmetry, END-3) and ABa/ABp (Notch from P2, second Notch from MS) are decided by cited mechanistic rules in `founders.bio`; every later fate is the observed lineage program (experimental, Sulston). The uncertainty report says so per level.

## 2026-09-11 (later) — cohort expression, altered-protein reconstruction, copy number

The three developments named as most valuable at the end of the therapeutics
work, implemented.

- **Cohort expression reference.** `genomeos cancer expression --study X`
  distils, per gene, how a cancer type behaves: the distribution of
  tumour-versus-normal z-scores across the study's tumours and the fraction in
  which the gene is raised above normal tissue. cBioPortal's reference-normal
  z-score profile is preferred because it is scored against the study's own
  normal samples and is therefore unit-free. Breast cancer separates MKI67
  (raised in 65.8% of 1,082 tumours) from CEACAM5 (38.8%), ERBB2 (11.0%), APP
  (2.2%) and RUNX1 (0.8%); EGFR sits at median -4.20, lower in tumour than in
  normal breast. `genomeos therapeutic --cohort X` uses it as the selectivity
  prior in place of the tissue-specificity class, capped at 0.6 because a cohort
  describes a cancer type and not a patient.
- **Three expression questions kept apart**: healthy tissue (Human Protein
  Atlas), this cancer type (cohort), this tumour (patient RNA). A patient TPM is
  never compared against a cohort's z-score or RSEM, because they are different
  scales; where a whole transcriptome is supplied the gene's rank within the
  patient's own sample is reported instead, which needs no unit.
- **Altered protein reconstruction.** The transcript id is now carried down from
  VEP, the coding sequence fetched from Ensembl, the HGVS coding change applied
  and both proteins translated and compared. Frameshifts yield their novel
  C-terminal stretch, in-frame indels their junction, substitutions their
  wild-type counterpart, and a premature stop is confirmed by sequence to create
  no new residue. Alignment is decided by the reading frame, not by protein
  length: a frameshift with no downstream stop can produce a protein of the same
  length that shares no sequence with the reference, and an earlier version
  mislabelled exactly that case. A frameshift now scores above a substitution on
  neoantigen strength because its residues have no wild-type counterpart at all.
  This also settles isoform disagreements: RUNX1 Y480* is called on a
  481-residue transcript against a 453-residue canonical UniProt entry, and
  rebuilding from ENST00000675419 confirms the truncation rather than reporting
  a mismatch.
- **VEP cache schema check**: a cached record written before a field existed is
  refetched instead of read as a missing value, so adding the transcript id does
  not silently lose it on an existing cache.
- **Copy number** (`--cnv gene<TAB>copies`): recorded as patient-derived,
  reported in the design dataset, and scored as at most 0.5 on tumour expression
  with the basis line stating that copy number bounds what a cell could display
  and never shows that it does. Patient RNA outranks it whenever both exist.
- **Data level** now reports both the level at which the chain of inputs first
  breaks and the higher-level inputs that were supplied anyway, so a run with
  HLA and copy number but no matched normal is not described as level 1 with no
  further explanation.
- **Tests**: 58 in tests/test_therapeutics.py, including the reconstruction of
  every variant class against a synthetic transcript whose translation can be
  read by eye, the frame-versus-length distinction, the cohort prior against a
  tissue-specificity class, unit-free ranking, and copy number bounding rather
  than standing in for expression. Suite excluding another session's in-progress
  files: 201 passed, 4 skipped.
- **Human body as counted populations** (step 6 of docs/ORGANISM-FROM-ONE-CELL.md): `genomeos data distil --only human_cell_turnover` streams Sender & Milo 2021 (26 cell types: counts, lifespans, turnover; 2.89e13 cells, 3.43e11 replaced per day) and generates `data/organisms/human/tissues.bio` (16 tissue populations with Cell Ontology ids, germ-layer shares, adult-count caps, daily loss and replacement, childhood growth rates from each tissue's share at birth). `data/organisms/human/body.bio` adds the Carnegie-stage timers, cleavage, trophectoderm / hypoblast / epiblast, the three germ layers and `resolution: populations`. `genomeos grow data/organisms/human/body.bio --until "20 yr"`: 2-cell at day 1, 8-cell at day 3, 1.8e12 cells at birth, 2.87e13 at 18 years, 2.84e13 at 20 years with 3.19e11 cells replaced per day (published 3.3e11; the totals are Sender & Milo's own counts, reproduced by growth and turnover rather than asserted). Runtime additions: population nodes (`count`, `fraction`, in-place growth, recurring fractional loss with the chain restarted after quiescence, stage events that make populations decide again, a recorded total for past-time counts and asserts). The report labels the whole thing low confidence: curated counts and inferred shares, no mechanism. Tests: `tests/test_human_body.py`.
- **Organism tab** in the web UI (`/api/grow`): choose the program, grow, see the alive-cell curve against the reference, fates, asserts, the diff against Sulston's lineage, the mechanism that fired, the tree and the uncertainty report.

## 2026-09-11 — project review and reorganised roadmap

- **Review against the code**: 68 commits in two days, 42 commands plus the `bio` toolchain, 16 web views, 220 tests passing and 5 skipped, lint clean on the whole working tree including the other sessions' uncommitted work. Per-chromosome coverage measured from `data/results`: fetched 18 of 25, UNKNOWN classified 25 of 25 (curated repeats on 14 of 24), CTCF domains 2, proteome 4, translation verified 3, knowledge graph 1, HG002 twin 1.
- **docs/ROADMAP.md rewritten** as the one organising document: scope and non-goals, the architecture with the BioLang/BioVM engine kept separable from the GenomeOS application (two boundary leaks named: `runtime/central_dogma.py` and `runtime/variant_effect.py` import the genome layer), eight areas each with goal, code, design, data, requirements, what is missing and ordered next steps, the table of data jobs still to run, milestones 0.9 (whole genome) to 2.0 (BioLang standalone), an improvement pool, and the working rules for several sessions in one checkout with file ownership per session.
- **Found and handed over**: five copies of the proteome job and two of the UNKNOWN job were running at once because the jobs registry de-duplicated only in memory across server restarts; genomeos-fe added the on-disk pid liveness check and killed the duplicates. The shared `.git/index` holds a stale tree (core modules marked deleted), so every session commits through a private index. The CI failure on PR #10 (embryo assert at 350 min, 194 cells against 240..310, on mean timers) was fixed by running timers with their measured spread; later runs are green.
- **Gaps recorded for the plan**: `experiment` block (in progress), Evidence explorer and Cell views never built, no haematopoietic mechanism module, human program is counts only, sequence grammar has no segment parser, the reader (open nodes per cell type) has no code, no retrospective therapeutic benchmark, no release tags (pyproject still 0.1.0), README family table stale, `data/jobs` tracked in git.
