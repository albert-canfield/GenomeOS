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
- **The proteome as a library**: `genomeos/lib/data/proteome.json.gz` (19,283 proteins, 2.3 MB, packaged) distilled from the 290 MB cache; `genomeos protein X --lib [--bio]` and `genomeos libs --proteome` work offline; a BioLang protein block per protein on demand.
- **Libraries tab shows the proteome library**: a card at the top with the packaged totals (proteins, function, pathway, partners, experimental structure) and a protein lookup that prints the record and its BioLang block; `/api/proteome_lib?gene=`.
- **Enhancer to gene, predicted (AlphaGenome feature b)**: `genomeos predict --element chr:START-END` deletes an element in its 1 Mb window and reads which gene moves in which tissue; job `enhancer_targets_chr21` scored 200 distal enhancers (63.5% move a gene by ≥ 0.1 log2; the strongest coding gene is the nearest TSS in the node for 67.8% of them and inside the node for 86.7%; 43 behave as silencers); predictions cached per element, shown as `predicted` next to `inferred` in `genomeos regulation` and the Flow tab; Progress card per chromosome; docs/NODES-READER-WRITER.md records it as the first test of the node model.
- **Your own genome as a local individual**: `genomeos individual import FILE.vcf[.gz] --name ME` splits any GRCh38 VCF into per-chromosome PASS files under the git-ignored `data/individuals/`; `lookup` carriers, the gene `report`, `twin build` and the new `individual genes` walk (every coding gene with the consequence of each coding SNV) see every imported person; Twin tab card with import, list and gene table; nothing leaves the machine.
- **Splice sites from AlphaGenome in the segment parser (feature c)**: `genomeos segments --chrom chr21 --predicted-sites` swaps the donor/acceptor matrices for the model's 1 bp splice-site tracks (calibrated: last exonic base / first exonic base, both strands) and nothing else; exact CDS segments 4.7%/5.2% → 40.1%/46.2%, exact splice sites 7.5%/6.9% → 51.6%/49.2% on chr21; gene precision unchanged, so the delimiters were the limit and the starts and coding model are the next; enhancer targets also scored on chr22 and chr1 (200 each, 92% and 87% inside the node).
- **A person's regulatory variants on one gene, predicted (AlphaGenome feature d)**: `genomeos individual predict --name ME --gene APP --chrom chr21` collects the person's variants inside the elements that reach the gene, scores each for its effect on the gene per tissue and sums per haplotype where phased; demo × APP: 39 variants in 99 elements, 4 move APP, the strongest in the enhancer feature b had already tied to APP; per-variant cache under data/knowledge, per-person results never under data/results; Twin tab button; all four AlphaGenome features built.
- **Reader lane in the block map**: a `reader` selector on the Blocks tab lists the cell types read on the chromosome; nodes fill with their open fraction (red edge when silent) and every coding gene carries a read/silent dot for that cell; the attributes travel with `/api/blocks` (`K562_open_fraction`, `HepG2_read`, …).
- **Gene starts anchored at ENCODE promoters**: `genomeos segments --promoter-anchored [REACH]` confines start codons to 5 kb downstream of a promoter-like element; with predicted splice sites, chr21 candidates fall 1,003 → 158 and precision doubles or triples (exons 61.9%, sites 66.0%, genes 53.8%) while sensitivity drops to the genes the registry marks (59.3%); both results committed side by side with the grammar baseline.
- **Reader on eleven cell types, every chromosome**: `scripts/reader_genome_wide.py` resumes per cell type and takes any ENCODE biosample; hepatocyte, cardiac muscle cell, H1, SK-N-SH, IMR-90, astrocyte, CD14-positive monocyte, GM12878 and keratinocyte join K562 and HepG2 (42% to 77% of coding genes read per cell); all eleven are reader lanes on the Blocks tab and columns on the Progress tab.
- **Reference check for an imported genome**: `genomeos individual check --name ME` applies the person's variants to GRCh38 haplotype by haplotype (statistics only), counts reference mismatches and gives an assembly verdict (≤ 0.5% matches; most calls disagreeing means GRCh37/hg19); stored under the person's directory, shown in the Twin tab list and button.
- **Individuals on the Flow tab**: one row of ticks per local person under the DNA lane (grey = variant inside the gene, red = protein-changing coding SNV), the changes listed with genotype and traceable with a click; `/api/flow` carries `individuals`.
- **Segment parser scored against canonical transcripts too**: `canonical_exon_sensitivity` in every segments result; with predicted splice sites 74.8% of chr21's 1,907 canonical CDS segments are found exactly (8.7% with the matrices, 56.3% anchored), so most of the exons still missing are isoform-specific.
- **Predicted enhancer targets as one genome-wide job**: `enhancer_targets_genome_wide` runs the chromosome jobs one after another, smallest first, retries on quota and never gives up; chr21, chr22, chr1–chr5 done (200 each).
- **Predicted enhancer targets on all 24 chromosomes**: 4,800 distal enhancers deleted one by one; 2,291 move a coding gene, 90.2% of those inside the element's CTCF node (79.8%–91.8% per chromosome), 71.2% exactly the nearest TSS, 1,157 silencer-like; `enhancer_targets_genome_wide.json` and a genome-wide block on the Progress card; the node model holds at genome scale as the model reads it.
- **Gene dossier extended**: `genomeos report` now says how many of the gene's enhancers AlphaGenome scored and which one moves it most, and which cell types read the gene (promoter open) from every reader result on the chromosome; per-chromosome enhancer jobs fold into one line on the Progress tab like the fetch jobs.
- **Kinetic pathways from BioModels**: `genomeos pathway R-HSA-… --kinetic [--model BIOMD… --knockout SYMBOL]` finds a curated ODE model by the pathway's name, runs it on the in-house SBML engine (now with functionDefinitions, rateRules, annotations and amount/concentration semantics) and compares a knockout run with the baseline on final levels and peaks; RAF1 out of the Hornberg2005 ERK cascade erases the MEK-P and ERK-P transients, MEK out of Kholodenko2000 stops the ERK-PP oscillation; Schoeberl2002 recorded as not reproducing.
- **Kinetic pathways on the Molecules tab**: under the pathway card, a kinetic model by name or BIOMD id, a knockout, the changed species table and a time-course chart with knockout curves dashed; `/api/pathway_kinetic`.
- **The second mammal**: `genomeos mouse --chrom chr19` fetches mouse mm10 chr19 (UCSC sequence, GENCODE vM25 models, ENCODE mouse cCREs) and infers its nodes with the human code; of 91 mouse nodes with ≥ 2 symbol-matched genes, 54% land in one human node and 97% in the same human neighbourhood (one node or adjacent nodes), 3% scattered; the node is a unit in both genomes and the CTCF-only boundary is the resolution limit (docs/NODES-READER-WRITER.md).
- **Second-mammal card on the Progress tab**: mouse chromosome totals and the node comparison (one node, same neighbourhood, scattered) with the evidence of each input.
- **Translation disagreements triaged by mechanism** (96 of 19,310 genes): 38 isoform/frame choice, 20 exon-boundary, 18 hg38 frameshift alleles, 13 hg38 nonsense alleles, 7 different product; the frameshift signature flags 29 canonical transcripts with the six mitochondrial mid-codon genes as the sanity check. Compiler: UniProt stop-word symbols (WAS) fall back to a plain gene query; 168 definitions marked symbol_match false. chrY verified. Jobs: a complete job reads done.
- **Genome-wide graph served from a cached build** (2.5 s to load 19k definitions, then 70 ms per neighbourhood); proteome coverage card shows one row per chromosome with a genome-wide total.
- **Milestone: the whole human proteome compiled** (19,478 coding genes, 25 chromosomes; sequence 99.1%, function 86.3%, pathways 58.2%, experimental structure 45.2%); translation verified on 19,249 genes (97.8% exact for some isoform, 99 real disagreements); knowledge graph 41,982 nodes / 327,024 edges. Results `proteome_genome_wide`, `graph_genome`, per-chromosome tables. chr4 repaired after the laptop sleep (92.6% → 99.9% sequence).
- **Proteome compile no longer calls Ensembl**: the genomic origin comes from the local GENCODE models (`origins_from_annotation`), which removes the one source that timed out at 60 s per gene; measured per source on a chr19 gene: UniProt 0.2 s, STRING 0.2 s, HPA 0.2 s, Ensembl lookup timed out.
- **Licence boundary closed on my side**: the engine imports nothing from the application. `runtime/central_dogma.py` takes Locus and Strand from `genomeos.coords` and has its own reverse complement; `bio` no longer imports the CLI: `genomeos/lang/tools.py` (Apache-2.0) carries load, check, compile, run, SBML and Boolean runs, and `bio.py` uses them.
- **Jobs are self-healing** (docs/JOBS.md): heartbeat per step, stall detection at 20 min, a supervisor thread in the server that restarts a dead or stalled auto-heal job from its saved results, a Heal button, and scripts that loop until every unit has its result; the proteome compile runs eight genes at a time.
- **Jobs show progress inside the current step** (proteome: 4.23 of 25 · chr19: 320/1,397) and logs are append-only.
- **Progress tab**: cards for the test human on every chromosome and the reader on every chromosome (`/api/individual`); finished fetch jobs fold into one line.
- **Reader on every chromosome** (K562 vs HepG2): K562 reads 14,828 of 20,094 coding genes, HepG2 13,015, 11,987 in both; K562 reads 0 of 61 on chrY (female line), HepG2 9. Summary `reader_genome_wide`; per-chromosome files local.
- **HG002 twin on every autosome** (`scripts/twin_genome_wide.py`, statistics only): 4,048,342 PASS variants applied, 0 reference mismatches, 840 overlaps skipped; result `hg002_twin_by_chromosome`.
- **Genome-wide UNKNOWN classification re-done with curated repeats**: 26,806 blocks, 1.009 Gb, 98.5% classified (from 96.6%); regulatory 33.6%, LINE 28.5%, gap 15.8%, unique intergenic 7.4%; long_orf gone. Table in docs/UNKNOWN.md.
- **Segment parser v1** (`genomeos segments --chrom chr21`): genes from learned signals alone, Viterbi over candidate starts/donors/acceptors/stops with codon log-odds; chr21 in ~2 min: gene sensitivity 89% precision 22%, splice sites 7.5%/6.9%, exact CDS segments 4.7%/5.2%. Recorded as the honest baseline in docs/SEQUENCE-GRAMMAR.md; predictions stay `predicted` and out of the block map.
- **Block map shows the reader**: open-chromatin blocks per cell type (one colour per cell type, chips to highlight) wherever `genomeos reader` has run.
- **Reader v1** (`genomeos reader --cell-type K562 --versus HepG2 --chrom chr21`): ENCODE DNase peaks per cell type over nodes, promoters and enhancers; K562 reads 135/221 coding genes on chr21, HepG2 118/221, 106 in both; RUNX1/ITGB2 only in K562, TFF1/TFF3/FTCD only in HepG2. docs/NODES-READER-WRITER.md.
- **Gene dossier**: `genomeos report GENE --chrom C` (Markdown or JSON), `/api/report`, Dossier button in Molecules: origin, DNA→RNA→protein with UniProt identity, regulation node and reachable enhancers, GTEx expression, protein facts, pathway loss, and the test human's variants inside the gene with coding ones traced.
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
- **Experiments (BioLang v0.3 `experiment` block)**: knockouts of maternal factors and signals run against the wild type with the same program and horizon; `when` clauses gain `absent`; `founders.bio` now states the founder fates as the genetics established them (ABa/ABp by Notch, EMS by SKN-1 without PIE-1, MS by POP-1, E without POP-1, C and D by PAL-1 without SKN-1, germline by PIE-1). `genomeos grow data/organisms/celegans/mutants.bio --experiments`: pop-1 → MS takes the E programme (192 descendants), skn-1 → MS, E and C take the C programme, pie-1 → P2 takes the EMS programme and no germline, apx-1 and glp-1 → ABp becomes ABa-like (528 descendants), pal-1 → C and D lose their fates; each with the published phenotype in `expect:` and asserts on the mutant. Terminal fates below the founders still come from the observed lineage program, which the report says. `genomeos/organism/experiment.py`, `tests/test_experiment.py`.

## 2026-09-11 — project review and reorganised roadmap

- **Review against the code**: 68 commits in two days, 42 commands plus the `bio` toolchain, 16 web views, 220 tests passing and 5 skipped, lint clean on the whole working tree including the other sessions' uncommitted work. Per-chromosome coverage measured from `data/results`: fetched 18 of 25, UNKNOWN classified 25 of 25 (curated repeats on 14 of 24), CTCF domains 2, proteome 4, translation verified 3, knowledge graph 1, HG002 twin 1.
- **docs/ROADMAP.md rewritten** as the one organising document: scope and non-goals, the architecture with the BioLang/BioVM engine kept separable from the GenomeOS application (two boundary leaks named: `runtime/central_dogma.py` and `runtime/variant_effect.py` import the genome layer), eight areas each with goal, code, design, data, requirements, what is missing and ordered next steps, the table of data jobs still to run, milestones 0.9 (whole genome) to 2.0 (BioLang standalone), an improvement pool, and the working rules for several sessions in one checkout with file ownership per session.
- **Found and handed over**: five copies of the proteome job and two of the UNKNOWN job were running at once because the jobs registry de-duplicated only in memory across server restarts; genomeos-fe added the on-disk pid liveness check and killed the duplicates. The shared `.git/index` holds a stale tree (core modules marked deleted), so every session commits through a private index. The CI failure on PR #10 (embryo assert at 350 min, 194 cells against 240..310, on mean timers) was fixed by running timers with their measured spread; later runs are green.
- **Gaps recorded for the plan**: `experiment` block (in progress), Evidence explorer and Cell views never built, no haematopoietic mechanism module, human program is counts only, sequence grammar has no segment parser, the reader (open nodes per cell type) has no code, no retrospective therapeutic benchmark, no release tags (pyproject still 0.1.0), README family table stale, `data/jobs` tracked in git.
- **Haematopoiesis, the first human mechanism module** (`genomeos/organism/haematopoiesis.py` → `data/organisms/human/haematopoiesis.bio`): HSC → MPP → CMP/CLP → MEP/GMP → erythroblast, myeloblast, monoblast, pro-B, thymocyte → erythrocytes, neutrophils, monocytes, B and T cells. Each compartment requires the transcription factors the genetics showed to be necessary (TAL1, RUNX1, GATA2, MYB, IKZF1, TCF3, GATA1, SPI1, CEBPA, KLF1, GFI1, IRF8, EBF1, PAX5, NOTCH1, TCF7, GATA3; mouse knockouts cited per compartment). Flows are solved from the measured daily outputs and lifespans upward (Sender & Milo 2021): outflow = output, residence sets the pool, amplification (inferred divisions per transit) sets what is needed from upstream; the stem-cell pool of 1e5 (Lee-Six 2018) then implies one HSC division per ~300 days (Catlin 2011: ~280). After three years every pool sits within 15% of its design value, 2.8e11 cells are replaced per day against 2.83e11 measured, and the program's asserts hold. Nine knockouts (`haematopoiesis_mutants.bio`, `genomeos grow ... --experiments`) reproduce the published phenotypes: GATA1 and KLF1 → no erythrocytes; SPI1 → no neutrophils, monocytes or B cells; CEBPA → no neutrophils; IRF8 → no monocytes; PAX5 → no B; NOTCH1 → no T; IKZF1 → no lymphoid; TAL1 → no blood. Runtime additions for it: recurring differentiation flows between populations (`differentiate` with `fraction` and `after`), pools that merge, growth steps above doubling, quarter-day steps because the steady state of a high-amplification compartment is a small difference between large flows. Tests: `tests/test_haematopoiesis.py`.

## 2026-09-11 — AlphaGenome as an optional feature, loaded disabled without the key

- **Optional features pattern** (decision D34): an external key or package never becomes a dependency; the feature is always listed, disabled without the key with the reason and the way to enable it, enabled the moment the key exists. `genomeos predict --status`, `/api/features`, an "Optional features" card in the Progress tab; the key lives in a git-ignored `.env` (read after the environment) and never reaches a result, a log or a job file. docs/DATA.md documents the four features the key enables: a variant effect (built), b regulatory blocks, c sequence grammar signals, d twin haplotypes (planned).
- **`genomeos predict chr:pos REF>ALT`** (feature a): the recommended RNA-seq variant scorer over a 1 Mb window, every gene by tissue track, effects above a log2 threshold as `predicted` rules capped at confidence 0.7, and a line saying what was scanned when nothing passes. Live checks with Albert's key: APP A673T (a coding change) has no expression effect (largest |log2FC| 0.015 across 24 genes and 371 tracks, as expected); the SORT1 enhancer variant rs12740374 raises PSRC1 (+0.81 log2 in liver), CELSR2 (+0.51) and SORT1 in liver tracks, the published liver eQTL of Musunuru 2010. Summary committed as `data/results/alphagenome_rs12740374.json`; tests in tests/test_alphagenome.py cover the status, the `.env` reader and the mocked scorer.
- **The Body as a process-bigraph process** (`genomeos/runtime/compose.py`, `BodyProcess`, `build_body_composite`): an organism program advances per composite interval beside the regulatory-network and ageing processes on one clock; counts, deaths, births, turnover and the census by cell type are reported as additive deltas. A network gates an organism through named factors (`gates={"TetR_high": {"species": "TetR", "threshold": 1}}`): while the species is at or above the threshold the factor is present in the organism's context, and `Body.set_factor` makes every population and resting cell decide again. The C. elegans embryo runs unchanged inside the composite (610 cells and 110 deaths at 14 h). Tests in `tests/test_compose.py` (skipped without the `compose` extra).

## 2026-09-11 — peptide/HLA binding, and everything credited

- **Two predictors, neither assumed** (`genomeos/therapeutics/binding.py`, decision D35). **MHCflurry 2** (Apache 2.0, `uv sync --extra hla`) is the default when installed; **NetMHCpan 4.1** (DTU academic licence) is called as a binary the user installed, found on `PATH` or at `$NETMHCPAN`, never shipped or downloaded. Both answer the existing `NeoantigenProvider` protocol, so nothing else in the pipeline changed. `--hla-predictor mhcflurry|netmhcpan|none`; `/api/features` and the Progress tab list both with their licence and state.
- **Percentile rank, with the field's thresholds**: binder at <= 2.0, strong binder at <= 0.5. Validated on known HLA-A*02:01 epitopes: SLYNTVATL scores 0.13 and KLVFFAEDV 0.23, while a poly-aspartate control is not a binder. The demo's APP N770K moves from "presentation unestablished" to "presentation predicted" with four mutation-spanning peptides.
- **Presentation is not binding, and the two predictors differ**: MHCflurry 2's presentation percentile includes its own antigen-processing model, so `predicted_processing` is recorded for it and presentation confidence reaches 0.4; NetMHCpan's rank is binding alone and stays at 0.25 with processing listed as missing. Immunopeptidomics remains the only evidence of real presentation, whichever predictor is used.
- **MHCflurry's own downloader cannot run on Python 3.13+** (it imports the removed `pipes` module), so `genomeos therapeutic --install-hla-models` reads the release URL from the package's own manifest and installs the models where MHCflurry looks. Nothing is redistributed.
- **ACKNOWLEDGEMENTS.md**: every library, optional model and data source credited with its licence, following AlphaGenome's own practice; plus the use terms, which are two separate things: the code is MIT, while AlphaGenome, NetMHCpan, COSMIC and some cBioPortal studies are non-commercial or academic, so a commercial user must check every source they enable. Whether the project's own licence should become non-commercial is recorded as an open question (D36).
- **Tests**: `tests/test_binding.py`, 9 tests that pass with or without either predictor installed (listing and licences, the honest default, an absent predictor answering unavailable with its reason and never an affinity, NetMHCpan output parsed by header rather than column position, the thresholds, the model-download hint, a working predictor's binders and computational evidence, and a binding-only score refusing to claim processing). Whole suite 244 passed, 4 skipped.
- **Housekeeping**: an orphaned 44-line `genomeos/therapeutics/hla.py` draft left untracked by the departed therapeutics session duplicated this design and stopped mid-sentence; it was removed in favour of the finished module.

## 2026-09-11 — what AlphaGenome is worth to GenomeOS

- **docs/ALPHAGENOME.md**: the model researched from its documentation, client source, the Nature paper (Avsec et al. 2026) and the model card, and mapped onto this project. What it provides: up to 1 Mb of sequence, 11 output modalities at 1 bp for expression, accessibility and splicing (128 bp for ChIP, 2,048 bp for contact maps), 5,930 human tracks each annotated with assay, biosample, ontology CURIE and GTEx tissue, 19 variant scorers returning genes by tracks with an effect size and a quantile rank against ~300,000 common variants, and in-silico mutagenesis.
- **Seven uses, ordered by value.** Variant effect per tissue (built); **enhancer to gene**, which is the biggest gain because 33.6% of the UNKNOWN space is regulatory and no element currently names its gene with evidence, only "nearest coding gene in the CTCF domain, inferred 0.4"; splicing at the 1 bp resolution our learned matrices cannot reach; domains from a predicted contact map as second evidence for a CTCF boundary; the reader extended to tissues ENCODE never assayed; in-silico mutagenesis to ask which bases of an unclassified block matter; and the twin, with the unphased-single-sequence caveat attached.
- **Shaped by the quota**: one variant per request and a service documented for thousands, not millions, of predictions, so the unit of work is a chromosome's regulatory elements as a resumable job, never a genome-wide sweep. Responses cache under data/knowledge/alphagenome/; only distilled summaries are committed. Open weights since January 2026 would lift that ceiling later without changing the evidence rules or the non-commercial terms.
- **DeepMind's science-skills collection reviewed and not adopted** (D38): GenomeOS already implements most of those databases natively with evidence and confidence, so wrapping them as skills would duplicate the work without the evidence model. Two things were taken from it: a source checklist of what we do not use yet (JASPAR, which the sequence-grammar plan already wants for promoter motifs, plus UniBind, Foldseek and the EBI Ontology Lookup Service), and the reverse idea, now in the roadmap pool: publish a skill that drives `genomeos` and `bio` so an agent reaches every layer through one tool.

## 2026-09-11 — relicensed: open, and protected (D39)

MIT was replaced by two licences split along the line the architecture already
draws, so the engine can spread and the application cannot be taken private.

- **BioLang engine** (`genomeos/lang`, `ir`, `runtime`, `std`, `coords.py`, `bio.py`, and the `.bio` programs): **Apache 2.0**. A language has to be embeddable; anyone may put BioLang inside their own tool, open or closed, as they would use Python or Node. Apache also grants patent rights, which MIT does not.
- **GenomeOS application** (everything else): **AGPL-3.0-or-later**. Use it, change it, publish research with it; run a modified version as a service and its users get the source. This is the protection: improvements come back instead of disappearing into a private product. The web interface now carries the source link AGPL section 13 expects.
- **Documentation** CC BY 4.0; `data/results/` summaries keep their upstream terms; the author reserves commercial dual licensing, and contributions are inbound-equals-outbound with a relicensing grant (CONTRIBUTING.md).
- **Non-commercial licences were considered and rejected** (PolyForm Noncommercial, CC BY-NC): they are not open source, they would bar exactly the internal research use by a company that adoption depends on, and they protect nothing the AGPL does not already protect.
- **The data question is separate from the code question**, and saying so plainly is the honest part: AlphaGenome, NetMHCpan, COSMIC and some cBioPortal studies are non-commercial or academic whatever GenomeOS's licence says, so a commercial user must check every source they enable. LICENSING.md holds the path map and the reasoning, NOTICE the summary, ACKNOWLEDGEMENTS.md every library, model and source with its terms.
- Both licence texts are byte-identical to the FSF and Apache originals; SPDX headers mark the package entry points; `pyproject` carries the classifiers and ships both texts with the wheel. Commits published before today remain MIT, which the documents state rather than hide.

## 2026-09-11 — the licence split has a dependency direction (D40)

Raised by the session that owns the genome and molecules layers, and correct:
the engine/application boundary stopped being only an architectural preference
the moment the two halves took different licences.

- **The rule**: the application may import the engine; the engine may never import the application. Apache 2.0 code sits happily inside an AGPL work, so `genomeos/genome` importing `genomeos/runtime` is fine. The reverse is not: an Apache-marked file that imports AGPL code is a work based on it and could not honestly be distributed as Apache. Adding such an import relicenses that file by accident. Stated in LICENSING.md and ROADMAP.md §2.
- **`runtime/variant_effect.py` fixed**: `Locus` and `Strand` now come from `genomeos.coords` (they always lived there; the genome layer only re-exports them), a variant is a structural `Protocol` that the application's VCF-backed class satisfies without the engine importing it, and reverse-complement is a four-line local helper instead of the application's `Sequence` class. No behaviour change; variant-effect, flow-trace, index-variant and twin tests unchanged.
- **Still crossing the line, in files owned by another session**: `runtime/central_dogma.py` imports `Sequence` from the genome layer, and `bio.py` imports three command functions from `genomeos.cli`, so the toolchain that is meant to stand alone currently depends on the application's command line. Both are recorded in the roadmap as the remaining engine-boundary work.
- **Space in the Body** (`genomeos/runtime/body.py`, using the spatial engine's fields): an organism may declare a grid and `field` blocks; cells hold sites, fields are stepped between events within the stable step, daughters take the parent's site and the nearest free site (along a stated direction first), a division with no free site waits another cycle (contact inhibition, counted), gradient `signal`s read a field at the cell's site and set factors that cells re-read every `sense` interval, and `migrate` moves cells along a direction or up a gradient into free sites, once or recurring. `data/demo/flag_organism.bio` grows Wolpert's French flag from one founder: the 30-site strip fills by 60 h, the gradient forms, and the bands come out Blue 0-7, White 7-12, Red 12-30, in order, with the program's asserts holding; knocking out the high-threshold signal removes the blue band. SPDX headers added to the files created this week (Apache for the engine, AGPL for genomeos/organism). Tests: `tests/test_space.py`.

## 2026-09-11 — the engine boundary holds, and a test keeps it that way

The last import from the Apache-licensed engine into the AGPL application is
gone, which closes the structural blocker for packaging BioLang separately.

- **Verified against the pushed commits, not a working tree**: reading `genomeos/lang`, `ir`, `runtime`, `coords.py` and `bio.py` straight out of dev, the runtime package imports nothing from the application. One crossing survived in `bio.py`: `from genomeos import __version__` reads `genomeos/__init__.py`, which is the AGPL package root, so a file marked Apache-2.0 was not honestly Apache.
- **`genomeos/version.py`** (D41) now holds the version string on the engine side and the package root re-exports it, so every existing `from genomeos import __version__` keeps working while `bio` reads the engine. It is also one fewer tie to cut at milestone 2.0, because the package root does not travel with the engine.
- **`tests/test_engine_boundary.py`** makes the rule permanent in three forms: every engine file is parsed for imports of the application, including imports inside functions; reading the package root is forbidden by name; and `import genomeos.bio` runs in a subprocess so the test asserts what is actually loaded rather than what is written. Nine modules load, none from the application.
- **Why it is a test and not a note**: since the relicensing the direction is a licensing requirement (D40), and an accidental import would relicense the file that made it. A grep in a review would eventually miss one; the run-time check cannot.

## 2026-09-11 — tree review: what was loose, and what it actually was

Asked to review and commit what the sessions had left behind. The shared git
index made it look like hundreds of files; against a clean index the truth was
six modified and 66 untracked, and most of that was not meant to be committed
at all.

- **`data/jobs` is now ignored** (63 untracked files, 304 KB, plus 11 that had been tracked since the Progress tab was built). Job metadata is per-machine runtime state: process ids, timings, heartbeats and console output that mean nothing in another checkout, and that change on every run. The registry that defines the jobs lives in `genomeos/jobs.py` and everything a job is worth keeping goes to `data/results`. Checked before untracking that a checkout without the directory still lists all 33 jobs without error, since the runtime creates it on demand.
- **`.env.example`** committed: the key goes in a git-ignored `.env`, and a template with an empty value is the difference between a documented feature and a guessing game.
- **`uv.lock`** updated for the `hla` extra added with the binding predictors.
- **`data/results/ccres_chrM.bed.gz`**, so the mitochondrial cCRE file is tracked like the 25 others.
- **One line in `cmd_verify`**, left uncommitted by the session that owns verification: it prints the same-protein-different-boundaries count that `molecules/verify.py` already writes. Confirmed with that session before committing.
- **A stray file named `T`**, created by a shell redirect of `C>T` in a variant argument, deleted. The lesson is in the command, not the repository: quote a variant.
- **BioForge takes experiment blocks as input** (`design` block, `genomeos/organism/forge.py`, `genomeos grow --design`): a design names the perturbations BioForge may try (factors to knock out or add, up to `at_most` combined, continuous knobs on timers and decisions), the targets to reach (assert grammar, loss = normalised distance from holding) and the constraints to keep; every candidate is an experiment run against the same program, seed and horizon; the answer is emitted as an `experiment` block with predicted evidence. `data/organisms/celegans/designs.bio`: "two intestinal founders" → POP-1, "two EMS founders and no germline" → PIE-1 (7 runs each). `data/organisms/human/haematopoiesis_designs.bio`: "lose neutrophils, keep the rest" → GFI1 (C/EBPa also removes monocytes because the GMP needs it); "no lymphocytes at all" → IKZF1 or TCF3 alone (both required by the common lymphoid progenitor), ranked above any pair. Runtime fix found on the way: a population changed only by an outgoing flow never decided again, so a pool that overshot its cap stopped amplifying for good; flows and pours now re-decide, and a decision swap keeps the pending step's time so the growth/outflow phase the balance depends on is not shifted. A toy design tunes a timer to a cell count. Tests: `tests/test_design.py`.

## 2026-09-11 — roadmap updated for the whole-genome protein layer

The protein milestone landed (3a64298, the session that owns molecules). The
numbers below were read back out of the committed results rather than taken
from the report, and the roadmap now carries them.

- **Whole human proteome compiled**: 19,478 coding genes over 25 chromosomes, from seven public databases. Sequence 99.1%, domains 99.0%, function 86.3%, interactions 81.5%, pathways 58.2%, predicted structure 98.2%, experimental structure 45.2%, disease 25.5%. The coverage numbers are measured, not assumed: they say what the databases actually answer for a human gene, and the honest reading of 45.2% is that fewer than half of human proteins have an experimental structure at all.
- **Translation verified against UniProt on 19,249 genes**: 90.9% identical to the canonical entry, 97.8% exact for some isoform, 99 real disagreements. That is the engine reading genes the way the curators do, at genome scale.
- **Knowledge graph genome-wide**: 41,982 nodes, 327,024 edges, 19,283 compiled proteins, largest component 15,165, 3,924 components.
- **Two gaps recorded rather than rounded away**: chrY is the one chromosome with no verification result (61 coding genes), and the 99 disagreements are counted but not yet explained one at a time. They concentrate mildly on the gene-dense chromosomes (chr19 12, chr11 11, chr1 10, chr6 10), which is what a per-gene isoform-choice problem would look like rather than a systematic bug.
- **Milestone 0.9 is now partial rather than pending**: fetch, UNKNOWN, domains, proteome and graph are at 25 of 25; the HG002 twin, chrY verification and the version bump are what remain.

## 2026-09-11 — the proteome becomes a library, and the shared index is made safe

- **The compiled proteome is now a packaged library** (de046e7, 451d01f): 19,283 human proteins distilled into `genomeos/lib/data/proteome.json.gz`, 2.3 MB, shipped inside the package. Function on 87.0%, pathways 58.7%, partners 63.6%, experimental structure 45.6%, disease 25.7%, and the 168 symbol mismatches marked rather than hidden. `genomeos protein TP53 --lib` answers in full with no network: name, function, location, domains, pathways, physical partners, 311 experimental structures and its diseases. That is the difference between a protein layer that needs seven live APIs and one that works on a laptop, in CI and on a plane.
- **Translation verification is at 25 of 25**, and the disagreements are triaged by mechanism rather than counted, with the frameshift signature and its mitochondrial sanity check built into `verify.py`.
- **The shared git index was made safe.** Three sessions share one `.git/index`, and it had drifted 77 entries from HEAD, carrying deletions of 52 files including `graph_genome.json` and every `proteome_chr*.json`: a plain `git commit` from any session would have removed the whole genome-wide protein result from the repository, 248,899 lines of it. Verified first that the working tree matched HEAD exactly and that every threatened file was present and committed, then reset the index. Nothing was lost, because a mixed reset touches staging and not content. The lesson is the one the project already knew and had not fully acted on: in a shared checkout the index is not a scratch pad, it is shared mutable state, and the only safe commit is through a private `GIT_INDEX_FILE`.

## 2026-09-11 — enhancer to gene, and the whole-genome data milestone

- **Feature b landed** (d19ec00): deleting a regulatory element in its 1 Mb window and reading which gene moves turns the target from "nearest coding gene inside the CTCF node, inferred at 0.4" into a named gene with a tissue and a magnitude. On 200 sampled chr21 distal enhancers, 63.5% move some gene by at least 0.1 log2; where a coding gene is named it is the nearest TSS in the node 67.8% of the time and inside the node 87.4%.
- **Two findings worth more than the headline.** A third of the time the nearest-gene heuristic names the wrong gene, which is the first quantified error rate we have for an inference the project has been making since domains existed. And 43 of the 127 elements that the ENCODE registry calls enhancer-like behave as silencers: expression rises when they are deleted. The registry names a class of element; it does not say which direction the element pushes, and now we can ask.
- **Read honestly**, the 87.4% agreement is a prediction agreeing with an inference, not a measurement, and the model's own stated weakness beyond 100 kb shows up in the data: elements within 20 kb of their inferred target name a gene 52% of the time, those beyond 100 kb only 32%. The document says so in the same table as the result.
- **The whole-genome data milestone is complete.** Every row of the roadmap's job table is now at full coverage: fetch and analysis, UNKNOWN classification with curated repeats, CTCF domains, the chromatin reader, the compiled and verified proteome, the knowledge graph, and the HG002 twin across all 22 autosomes with 4.05 M variants applied and zero reference mismatches. chrX and chrY are recorded as out of scope rather than pending, because the GIAB benchmark does not phase them. What is left of 0.9 is the version bump, the tag and a README refresh, which is release work rather than science.

## 2026-09-11 — better splice sites, and the thing they did not fix

- **Feature c landed** (a9d49a5): the segment parser can take AlphaGenome's 1 bp splice-site tracks in place of the learned donor and acceptor matrices, calibrated to the last and first exonic base on both strands, with nothing else changed. On chr21, scored against every coding transcript: exact splice sites rise from 7.5% to 51.6% and exact CDS segments from 4.7% to 40.1%, an eight-fold gain on both.
- **Gene-level precision does not move**: 21.9% before, 20.0% after, with sensitivity slightly lower (89.1% to 84.2%). That negative result is worth more than the positive one. Splice sites were never the binding constraint on finding a gene; the start codon and the coding model are, and the next work goes there instead of into better sites. Holding both baselines side by side is what makes the claim checkable, so `segments_chr21.json` stays committed next to `segments_chr21_predicted_sites.json`.
- **A person's own genome is now a first-class input** (8654354): `genomeos individual import FILE.vcf` splits any GRCh38 VCF into per-chromosome PASS files under a git-ignored directory, and carrier lookup, the gene report, twin build and the gene-by-gene walk all see it. Nothing leaves the machine. That was the last item standing between HG002 as a test subject and GenomeOS being usable on the genome a geneticist actually has.
- **Predicted enhancer targets extended** to chr22 and chr1 (200 elements each): the strongest coding gene sits inside the CTCF node 91.8% and 86.7% of the time, against 87.4% on chr21, so the node model holds across three chromosomes rather than one. The job is registered for every chromosome.

## 2026-09-11 — all four AlphaGenome features, and privacy enforced by the repository

- **Feature d landed** (9a8ca6d): `genomeos individual predict --name ME --gene APP` collects the person's own variants inside the elements that actually reach the gene, promoters first, then the enhancers feature b tied to that gene, then near ones, then the rest of the node, scores each for its effect on the gene per tissue, and sums per haplotype where the genotype is phased. On the demo person, 39 variants across 99 elements, 4 of which move APP; the strongest is a 6 bp insertion in the same enhancer the deletion job had independently tied to APP, which is two different questions agreeing.
- **All four features are now built**: variant effect per tissue, enhancer to gene, splice sites for the segment parser, and a person's own regulatory variants.
- **A privacy gap closed by measurement rather than by trust.** Checked, rather than accepted, that nothing person-shaped is tracked: per-person predictions are cached under a git-ignored directory and never written to `data/results`, and the imported genomes are ignored too. But `data/twins` was not ignored, and a twin built from an imported genome holds that person's chronological age, telomere length, epigenetic age, variants and modifiers. It is now ignored except for the two open-consent public test-subject files, so a personal twin cannot be committed by an absent-minded `git add`. The rule the project keeps relearning: a privacy property that depends on everyone remembering is not a property, it is a hope; put it in `.gitignore` where the tool enforces it.

## 2026-09-11 — two signals, one conclusion: the coding model is the limit

- **Gene starts anchored at curated promoters** (f0cb643): confining a gene start to 5 kb downstream of an ENCODE promoter-like element cuts chr21 candidates from 1,003 to 158 and raises precision 2.5-fold at gene level (21.9% → 53.8%), while sensitivity falls to 59.3%, because it can only find the genes the registry already marks. A trade, recorded as one, with all three parser runs kept side by side.
- **What the three runs say together.** Better splice sites multiplied exon and site accuracy eight-fold and left gene finding untouched. A better start signal finally moved gene precision but paid for it in sensitivity. Neither signal found a gene the grammar was missing. The codon model, which decides between a real reading frame and a merely plausible one, is the limit; that is now a measured claim rather than a suspicion, and it is where area B goes next.
- **The reader on eleven cell types, genome-wide**: 42.0% of coding genes read in keratinocyte up to 77.1% in hepatocyte, resumable per cell type and accepting any ENCODE biosample. The block map gained a reader lane: nodes fill with their open fraction, silent ones take a red edge, and every coding gene carries a read or silent dot for the chosen cell. The same chromosome, seen through a different cell, is a different map, which is the point the node model has been making in prose since it was written.
- **The AlphaGenome daily quota ran out, and nothing broke.** The waiting chromosome job retried rather than failing, which is exactly the behaviour the quota-shaped job design was built for and the first time it was tested by the real limit rather than by a unit test.

## 2026-09-11 — correction: roadmap edits that were logged but never landed

The previous two entries described roadmap changes that were not in the file.
An edit script asserted its way out part-way through, writing nothing, while
the commit that followed carried only the progress entry, so the log described
a state the roadmap did not have. Everything is now applied and each edit was
verified present after the fact rather than assumed:

- area B: the three-run segment-parser table, the reader lane, and a next list
  headed by the coding model;
- the §4 reader row (eleven cell types, 42.0% to 77.1% of coding genes read);
- the §4 enhancer row (three chromosomes, and the quota exhaustion the job
  handled by waiting);
- a §4 row for the segment parser itself, which had none.

The lesson is the same one the shared index taught this morning: a step that
reports success in one place and silently does nothing in another is worse
than a step that fails loudly. A script that edits a document should verify
its own result, and the person running it should read the file, not the exit
line. Both mistakes today were caught by reading the artefact instead of the
report, which is the habit to keep.

## 2026-09-11 — the node model tested on the whole genome, and a scoring artefact removed

- **Enhancer to gene, all 24 chromosomes** (84c3318): 4,800 distal enhancers deleted one at a time. 2,994 move a gene, 2,291 name a coding gene, and **90.2% of those sit inside the element's own CTCF node** (79.8% to 91.8% per chromosome), 71.2% exactly the nearest transcription start. The node model has been argued for since docs/NODES-READER-WRITER.md was written; it is now tested at genome scale, and the boundary predicts where an element acts nine times out of ten. 1,157 elements that the ENCODE registry calls enhancer-like rise when deleted, so they act as silencers: the registry names a class of element, not a direction of effect.
- **A scoring artefact removed** (4c034f6): the segment parser was being scored against every coding transcript, so an exon that exists only in a minor isoform counted as a miss. Scored against each gene's canonical transcript, the same run finds 74.8% of CDS segments exactly, against the 40.1% the earlier figure reported. Most of what looked like failure was the parser not reproducing isoform-specific exons, which is a much smaller problem than not finding the gene's exons at all. The roadmap table now carries both columns, because the pessimistic number is still the right one for "reproduce GENCODE" and the optimistic one is the right one for "find the gene".
- Neither number changes the conclusion from this morning: gene-level precision is still the weak axis and the coding model is still the limit.
- **Not recorded here**: kinetic pathways through the SBML engine were announced as about to land and are not on dev yet, so the roadmap says nothing about them. After this morning's slip, a claim goes in the documents when it is in the repository and has been read back out of it, not when it is reported.

## 2026-09-11 — pathways that run, and a failure kept next to them

- **Kinetics where a curated model exists** (e7b438a, read back from the committed results). `genomeos pathway ID --kinetic` finds a curated ODE model by the Reactome pathway's name, runs it on the in-house SBML engine, and compares a knockout against the baseline on final levels and peaks; a gene symbol reaches species named `x1`/`x2` through the model's MIRIAM annotations, which is what makes a curated model addressable by biology rather than by variable name. Taking RAF1 out of Hornberg2005 erases the downstream transients completely (peaks 0.201 and 0.562 to zero). Taking MEK out of Kholodenko2000 stops the ERK oscillation: Erk2-PP peak 298.8 to 10.0, final 37.0 to 0. The engine gained function definitions, rate rules, species names and annotations, and amount-versus-concentration semantics.
- **This closes the gap the pathway layer has carried since it was built.** Reactome's export has no rate laws, so pathways ran as reachability graphs and a knockout could only say which reactions became unreachable. Where BioModels has a curated model of the same biology, a knockout can now say how much of which species, when, and whether the oscillation survives.
- **The failure is committed next to the successes.** Schoeberl2002, 100 species and 125 reactions, does not reproduce on the in-house engine: the cascade species stay flat. Keeping that result file is what turns "adopt libRoadRunner eventually" from an aspiration into an item with a reason and a size attached, and it marks where our engine stops.
- Recorded only after reading all four result files out of the repository, following this morning's slip. The announcement came before the push; the roadmap waited for the push.

## 2026-09-11 — the second mammal, and what it does and does not show

- **`genomeos mouse --chrom chr19`** (64f2f98) brings mouse mm10 chromosome 19 through the same three doors as any human chromosome, UCSC sequence, GENCODE vM25 models and the ENCODE mouse element registry, and then infers its nodes with the *unchanged* human code. That reuse is the first result: the node inference was written against human CTCF elements and needed no species-specific change to produce 371 mouse nodes, 256 of them with coding genes, median 109 kb, against a human median of the same order.
- **The comparison, read carefully.** Of the 91 mouse nodes holding two or more genes whose symbol matches a human gene, 49 land inside one human node and 39 more split across *adjacent* human nodes, with 3 scattered. Synteny between mouse and human is textbook, so "the genes stay together" is not news and would not have been worth running. The number that carries information is the shape of the failures: when a mouse node does not map to a single human node, its pieces are next to each other 39 times out of 42. A unit that fragments into neighbours is a unit whose boundaries are drawn slightly differently, not a unit that does not exist. The node behaves as a real object in both genomes and the CTCF-only boundary is the resolution limit rather than the biology.
- **Stated as a lower bound, on purpose.** Orthology here is gene-symbol identity, which misses every renamed and every one-to-many orthologue, and the result file says so in its evidence block. A stronger run through Ensembl Compara would move the numbers up, not down, which is why the conclusion survives the weak method.
- Also landed: the kinetic pathway view on the Molecules tab (a443cec), with the knockout curves drawn dashed against the baseline.

## 2026-09-11 — 0.9

Version 0.9.0. The whole genome is decoded end to end, and the parts of the
model that could be tested have been tested rather than asserted.

- **Every chromosome**: 25 fetched and analysed, UNKNOWN space classified at 98.5% with curated repeats, CTCF nodes on 24, open chromatin read for eleven cell types (42.0% to 77.1% of coding genes read, depending on the cell).
- **The whole proteome**: 19,478 coding genes compiled from seven public databases, 19,249 verified against UniProt (90.9% identical to the canonical entry, 97.8% exact for some isoform), the disagreements triaged by mechanism, and the result packaged as a 2.3 MB library that answers with no network at all.
- **The structural model tested twice, from different directions**: 4,800 regulatory elements deleted one at a time, and the gene that moves sits inside the element's own node 90.2% of the time; in mouse, a node that splits across human nodes splits into adjacent ones 39 times out of 42. The node stopped being a design argument and became a measured object.
- **A person's own genome** as a first-class input that never leaves the machine, with predicted effects of their own regulatory variants on a named gene.
- **Pathways that run**: Reactome as reachability, and as kinetics where a curated model exists, with the model size at which our engine fails recorded next to the two that work.
- **The engine is separable**: BioLang imports nothing from the application, enforced by a test rather than by care, and licensed apart from it (Apache 2.0 against AGPL-3.0-or-later).
- Suite: 274 passed, 4 skipped; 22 BioLang programs testing themselves, 62 checks; lint and format clean.

What 0.9 does not contain, stated plainly: no benchmark showing the therapeutic
pipeline recovers known targets, no published perturbation reproduced, and the
human body above haematopoiesis is still counts rather than mechanism. Those
are 1.0 to 1.2, and they are all validation rather than construction.
