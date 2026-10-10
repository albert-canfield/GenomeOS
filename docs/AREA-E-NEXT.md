# Area E, the next mechanism: the choice, written before anything was built

Area E (From one cell to an organism) has had no owner and no model result since 2026-09-22.
Two candidates were on the table. This note chooses between them and says why, before a rule
was written or a number fitted. Written 2026-10-10 by lane-organism; the authoritative,
machine-readable form of the list below is `REGISTRATION` in `scripts/intestine_registration.py`,
committed before the module's numbers were fitted so that the closed list provably cannot grow.

## The two candidates

**(a) The worm's next real number.** Area E's standing count is *522 with the lookup as fallback,
381 honest, 136 free, **0 decided without the lookup having written the fate first*** (`d5215bf`),
and the area's own record says the next real number "is not in the Ma 2021 atlas but in a
measurement about the decision to divide". The terminality gate was asked and answered in
September: no (`de7a135` registration, `9d266a3` result). What the worm needs next is therefore
either contact and geometry at division — new runtime machinery, which this lane is forbidden to
write and which is not a one-lane job — or a measured cell-cycle regulator, which needs a data
source nobody in the project has yet found. Both are blocked on something that is not a model.

**(b) A second human mechanism module after haematopoiesis.** Outside haematopoiesis the human body
program is counts: `body.bio` and the generated `tissues.bio` are Carnegie stages, germ-layer shares
and sixteen adult populations from Sender & Milo 2021, with no factor deciding anything. The
haematopoiesis module is the one place where a lineage choice is gated by a named factor, and the
pattern it set — a requirement per compartment taken from published loss-of-function genetics, flows
solved from measured outputs, published mutants as `experiment` blocks run by `bio test` — is
reusable wherever a human lineage has that quality of genetics.

## Chosen: (b), and specifically the small-intestinal crypt-villus lineage

The coordinator's recommendation was (b) and it is taken, with the lineage chosen as follows.

| Why this lineage | Evidence it rests on |
| --- | --- |
| The gut epithelium is already the largest-turnover population in the counts-only half of the body program (1.97e11 cells, 4.03e10 replaced per day in `tissues.bio`), so mechanism lands where the program is weakest | `data/organisms/human/tissues.bio`, from Sender & Milo 2021 |
| Its allocation genetics is of the same quality as blood's: one stem compartment, one binary absorptive/secretory branch, and single-gene loss-of-function phenotypes published for each branch and for three of the four secretory types | the eight entries of the closed list |
| Its quantities come from the same source the rest of the human program uses, so a fitted number can be stated against a published one | Sender & Milo 2021: small-intestine epithelia 1.2e11 cells, lifespan 5.3 d, 2.2e10 cells/day |
| Erythroid-adjacent systems are already covered by haematopoiesis, and this lineage shares no factor with it | `data/organisms/human/haematopoiesis.bio` |

Not chosen within (b): the colon (its one well-sourced factor, KLF4, is colon-specific and the
composition differs — no Paneth cells), the stomach, and the tuft lineage. The reasons are in the
exclusions below.

## The closed, sourced list, written before any rule was tuned

Eight entries are used as gates. Each names a factor, the step it gates, and a published
loss-of-function phenotype; every citation was fetched and read in this lane rather than recalled.

| # | Factor | Step it gates | Published phenotype of its loss | Source |
| --- | --- | --- | --- | --- |
| S1 | TCF7L2 (Tcf4) | the stem compartment exists | no proliferative compartment is maintained in the prospective crypts; the neonatal epithelium is entirely differentiated villus cells | Korinek et al. 1998, Nat Genet 19:379 |
| S2 | ASCL2 | the stem compartment exists | induced deletion in adult small intestine: Lgr5 stem cells disappear within days | van der Flier et al. 2009, Cell 136:903 |
| S3 | RBPJ (Notch effector) | the transit-amplifying compartment exists | conditional removal converts proliferative crypt cells into post-mitotic goblet cells | van Es et al. 2005, Nature 435:959 |
| S4 | ATOH1 (Math1) | the secretory progenitor | goblet, enteroendocrine and Paneth cells are depleted; enterocytes are unaffected | Yang et al. 2001, Science 294:2155 |
| S5 | NEUROG3 (Ngn3) | enteroendocrine output | null mice generate no intestinal endocrine cells and lack endocrine progenitors; the other epithelial types differentiate properly | Jenny et al. 2002, EMBO J 21:6338 |
| S6 | GFI1 | goblet and Paneth output | Gfi1-/- lack Paneth cells, have fewer goblet cells and supernumerary enteroendocrine cells | Shroyer et al. 2005, Genes Dev 19:2412 |
| S7 | SPDEF | goblet and Paneth output | loss severely impairs maturation of goblet and Paneth cells and accumulates immature secretory progenitors | Gregorieff et al. 2009, Gastroenterology 137:1333 |
| S8 | SOX9 | Paneth output | in its absence Paneth cells are not formed, other types unaffected (Bastide: Paneth disappear, goblet decreased) | Mori-Akiyama et al. 2007, Gastroenterology 133:539; Bastide et al. 2007, J Cell Biol 178:635 |

Excluded, with a source and a reason — the list is closed and these do not enter the module:

- **POU2F3** — Pou2f3-/- mice lack intestinal tuft cells (Gerbe et al. 2016, Nature 529:226). The
  factor is sourced; the *quantities* are not. Neither the tuft share of daily output nor a tuft
  residence time has a measurement this lane could cite, so the lineage would be a population whose
  every number was invented. Excluded for want of a sourced quantity, not for want of genetics.
- **KLF4** — Klf4-/- mice show a 90% decrease in colonic goblet cells (Katz et al. 2002,
  Development 129:2619). Colon-specific; this module is the small intestine.
- **HES1** — Hes1-/- shows precocious and excessive endocrine differentiation in the developing
  stomach and gut (Jensen et al. 2000, Nat Genet 24:36). Excluded **as a gate**: the phenotype is an
  *increase*, and a population module made of flows can delete a lineage but cannot re-specify one,
  so gating the absorptive branch on HES1 would predict a loss no paper reports.
- **CDX2** and **LGR5** were considered and are excluded with no citation claimed, because nothing
  was fetched for them in this lane: CDX2's published phenotype is a change of organ identity rather
  than the loss of a lineage inside the epithelium, and LGR5 marks the stem cell without a
  loss-of-function phenotype that removes the compartment.

## What the formalism cannot express, registered in advance

Two of the eight phenotypes are re-allocations, not deletions: GFI1 loss *increases* enteroendocrine
cells and RBPJ loss *converts* progenitors into goblet cells. A flow module has no fate
re-specification, so these two can only be registered as partial matches, and they are registered
that way before the run — not explained afterwards.

## Which numbers are fitted, and to what

Stated before the fit, so that the result cannot be read as a prediction:

- **Fitted (one knob):** the villus residence shared by enterocytes and enteroendocrine cells, fitted
  so that the module's total epithelial pool equals the published human count (1.2e11 cells). The
  mouse-measured absorptive turnover is about 3 d; whatever the fit returns is a fitted number and
  the module says so.
- **Inferred, from sourced bounds, not fitted:** the lineage shares of daily output (enterocytes
  >95% of villus cells and enteroendocrine <1% of epithelial cells, Cheng & Leblond 1974,
  Am J Anat 141:461; goblet ~4% in the duodenum, Karam 1999, Front Biosci 4:D286 — a review; Paneth
  as the remainder) and the transit-amplifying transit time and division count.
- **Measured, used as given:** daily output and the epithelial count (Sender & Milo 2021, Nat Med
  27:45); goblet turnover ~3 d (Cheng 1974, Am J Anat 141:481); Paneth turnover ~15 d (Cheng 1974,
  Am J Anat 141:521); stem-cell division about daily, mouse (Snippert et al. 2010, Cell 143:134).
- **Not claimed:** the human stem-cell pool. It is a consequence of the flows, reported with the
  crypt count it would imply, and checked against nothing.

The knockouts are the result, because the gates come from the genetics above and nothing about a
mutant was fitted. What the module reproduces is therefore the *pattern of which lineages disappear*,
and the total-count assert only checks that the one fitted knob closed.

## The `reader` construct

`docs/BIOLANG-v0.5-READER.md` (`5682e6c`) holds a rule through a closed element and reports the
element-absent state. Nothing in this module is gated on a *genomic element* whose closure is
measured — the gates are factor-presence requirements, exactly as in haematopoiesis — so the
construct is not used here. Using it would be decoration.
