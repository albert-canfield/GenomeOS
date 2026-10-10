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

---

# The result (2026-10-10, lane-organism, `460c46f` registration, `7646bd7` module)

`data/organisms/human/intestine.bio`, generated by `genomeos/organism/intestine.py`, with its
knockouts in `intestine_mutants.bio` and its tests in `tests/test_intestine.py`.

## The wild type: the fit closing, not a prediction

One year of the module, against the numbers it was built from:

| | module | published | kind |
| --- | --- | --- | --- |
| epithelial cells | 1.2008e11 | 1.2e11 (Sender & Milo 2021) | the fit closing |
| cells replaced per day | 2.2e10 | 2.2e10 (Sender & Milo 2021) | the fit closing |
| villus residence | 4.70 d, **fitted** | about 3 d in the mouse; 5.45 d is what the two numbers above imply | the one knob |
| stem pool | 6.88e8 | nothing: it is a consequence of the flows | claimed of nothing |
| crypts implied at 14 Lgr5 cells each | 4.91e7 | nothing was fetched to check it against | reported only |

Every compartment lands within 0.1% of its design (stem 6.88e8, transit-amplifying 1.588e10,
secretory progenitor 1.101e9, enterocytes 9.828e10, goblet 2.642e9, enteroendocrine 8.276e8, Paneth
6.604e8), and all five `assert` lines pass. That is the closure of one fitted number and it is worth
exactly what the brief says it is worth: nothing, as a result.

## The knockouts: 8 registered, 8 reproducing the registered pattern

No mutant number was fitted, so this is the module's prediction. 32 lineage checks, no exception.

| mutant | registered | run | registered as |
| --- | --- | --- | --- |
| TCF7L2 | no lineage at all | all four at 0 | full |
| ASCL2 | no lineage at all | all four at 0 | full |
| RBPJ | the proliferative compartment lost | all four at 0 | **partial in advance**: the published conversion into goblet cells is not expressible |
| ATOH1 | goblet, Paneth, enteroendocrine gone; enterocytes kept | exactly that | full |
| NEUROG3 | enteroendocrine alone gone | exactly that | full |
| GFI1 | goblet and Paneth gone; enterocytes and enteroendocrine kept | exactly that | **partial in advance**: "fewer goblet cells" is overstated by a module that can only delete |
| SPDEF | goblet and Paneth gone, and the secretory progenitors accumulate | exactly that: 6.88e9 against the wild type's 1.10e9 | full, with the registered extra |
| SOX9 | Paneth alone gone | exactly that | full |

## Two things the registration got wrong about the formalism

Both are pinned in `tests/test_intestine.py` so they cannot be read later as claims about the gut.

1. **The registration said a flow module cannot express the supernumerary enteroendocrine cells of
   Gfi1-/-. It can, and it does.** The progenitors that cannot leave by the blocked routes leave by
   the open one: enteroendocrine cells stand at 6.25 times the wild-type pool without GFI1 (and the
   same without SPDEF). The direction is the published one and the mechanism is mass conservation,
   not fate re-specification; the 6.25-fold is a property of the module's shares and no measured fold
   change was fetched to compare it with.
2. **The same conservation inflates what is left when a branch is deleted.** Without ATOH1 the
   enterocytes stand at 1.64 times the wild type, while Yang et al. 2001 report them unaffected. The
   registered check (the lineage is kept) passes; the magnitude does not match the paper. It happens
   because each compartment's growth term was solved for the wild-type outflows, so removing an
   outflow leaves the growth term too large until the backs-up cap stops it.

## What stays open

- The colon and the stomach are not modelled, so the body program's `GutEpithelium` population
  (1.97e11) is still counts; this module covers the small-intestinal part of it (1.2e11).
- The tuft lineage waits on a sourced share of output or a residence time (POU2F3 has the genetics).
- Two published phenotypes in the closed list are only half-expressible, and the reason is the
  formalism: a `differentiate` flow cannot re-specify a fate, so a conversion reads as a deletion.
- A graded gate. Every gate here is presence-or-absence, so a hypomorph, a heterozygote or a partial
  knockdown cannot be stated at all — the same limit the mesoderm lane found in the rule numbers
  (`b1e8839`).
- Composition with `body.bio` was not wired, deliberately and for the reason haematopoiesis is not
  wired either: `body.bio` counts the same cells as populations of `tissues.bio`, so importing a
  mechanism module beside it would count the small-intestinal epithelium twice. The module stands
  beside the body program exactly as `haematopoiesis.bio` does.

---

# The colon, written before anything was built (2026-10-10, lane-colon)

The line left open above — "the colon and the stomach are not modelled, so the body program's
`GutEpithelium` population (1.97e11) is still counts" — is what this section closes for the colon.
The authoritative, machine-readable form of everything below is `REGISTRATION` in
`scripts/colon_registration.py`, committed before a single number of the module was fitted, so that
the closed list provably cannot grow. Every citation here was fetched in this lane (PubMed
E-utilities and PMC full text) and the sentence each gate rests on is quoted in that file.

## How the colon differs from the small intestine

| | small intestine (the module above) | colon (this module) |
| --- | --- | --- |
| geometry | crypt and villus; cells are born in the crypt and shed from the villus tip | crypt and a flat surface epithelium, no villus: "The colon has crypts but no villi" (van Es et al. 2012, Mol Cell Biol 32:1918) |
| mature lineages | four: enterocyte, goblet, enteroendocrine, Paneth | three: colonocyte, goblet, enteroendocrine (Katz et al. 2002: "differentiating to form one of three cell types: colonocytes, goblet cells and enteroendocrine cells") |
| the secretory-cell niche at the crypt base | Paneth cells, residence ~15 d | no Paneth cells; Reg4+ **deep crypt secretory cells** are intermingled with the Lgr5+ stem cells and serve as the Paneth equivalent (Sasaki et al. 2016, PNAS 113:E5399) |
| goblet share | ~4% of duodenal epithelial cells (Karam 1999, a review) | **20% of colonic epithelial cells are mature goblet cells** in the wild-type mouse colon (Katz et al. 2002, measured) |
| goblet terminal differentiation | GFI1 and SPDEF | **KLF4**, excluded from the small-intestinal module for being colon-specific, is a gate here |
| stem cells per crypt | ~14 Lgr5+ (Snippert et al. 2010) | ~20 Lgr5+ (Sasaki et al. 2016), plus 12-14 Reg4+ deep crypt secretory cells in the caecum and distal colon |
| counts (the same source the body program uses) | 1.2e11 cells, 2.2e10/day, lifespan 5.3 d | **6.0e10 cells, 1.7e10/day, lifespan 3.4 d** (Sender & Milo 2021, Nat Med 27:45) |

Two differences decide the shape of the module. First, with no villus and no Paneth cells there is no
long-lived resident lineage: every mature colonic cell rides the same crypt-to-surface conveyor and is
shed at the surface, which is why one shared residence is a defensible inference here and was not in
the small intestine. Second, the goblet lineage is 20% of the epithelium rather than 4%, so every
consequence of deleting it is five times larger.

## The closed, sourced list: seven gates, each with a fetched COLON phenotype

A factor well sourced in the small intestine is not a source for the colon. Each gate below rests on
a sentence about the colon, quoted in the registration.

| # | Factor | Step | Published colon phenotype of its loss | Source | Registered as |
| --- | --- | --- | --- | --- | --- |
| C1 | TCF7L2 (Tcf4) | the stem compartment | an adult inducible deletion abolishes the proliferating cells of the colon within 7 days | van Es et al. 2012, Mol Cell Biol 32:1918 | **contested**: Angus-Hill et al. 2011, PNAS 108:4914 report adult colonic crypts lacking Tcf4 *enlarged*, with proliferative and differentiated cells present |
| C2 | RBPJ (Notch effector) | the transit-amplifying compartment | the proliferative crypt compartment is lost completely and converts into goblet cells, in the colon as in the proximal intestine | Pellegrinet et al. 2011, Gastroenterology 140:1230 | partial: a conversion is not expressible |
| C3 | ATOH1 | the secretory progenitor | the secretory lineages are nearly absent in the colon | Lo et al. 2016, Cell Mol Gastroenterol Hepatol 3:51 | full |
| C4 | NEUROG3 | enteroendocrine output | no enteroendocrine cells anywhere along the intestinal tract, the colon included, goblet cells remaining | Mellitzer et al. 2010, J Clin Invest 120:1708; Jenny et al. 2002, EMBO J 21:6338 | full |
| C5 | GFI1 | goblet output | about a third the normal number of goblet cells, the reduction also seen in the colon; colonic PYY cells sixfold up | Shroyer et al. 2005, Genes Dev 19:2412 | partial: a reduction read as a deletion |
| C6 | SOX9 | goblet output | a strong reduction of the colonic goblet population, enteroendocrine and absorptive cells comparable | Bastide et al. 2007, J Cell Biol 178:635 | partial: a reduction read as a deletion |
| C7 | KLF4 | goblet output | Klf4-/- newborns have 90% fewer colonic goblet cells, colonocytes and enteroendocrine cells normal | Katz et al. 2002, Development 129:2619 | **contested**: Pellegrinet et al. 2011 found an adult inducible deletion changed no goblet number in small intestine or colon |

Excluded, each with the fetched reason — the list is closed and none of these enters the module:

- **ASCL2** (van der Flier et al. 2009, Cell 136:903) — the record states the *small intestine*;
  carrying the gate across would be an assumption, not a source.
- **SPDEF** (Gregorieff et al. 2009, Gastroenterology 137:1333) — a gate of the small-intestinal
  module; its abstract says "intestinal epithelium" without naming the colon and the publisher's full
  text refused the fetch (HTTP 403), so there is no colon sentence to quote.
- **REG4** (Sasaki et al. 2016) — the colon's Paneth equivalent and the one lineage this module would
  add to the small-intestinal pattern, but the experiment ablates the *cells* with diphtheria toxin
  rather than removing a factor, and no share of output or residence time exists for them. The
  measured 12-14 per crypt is recorded in the registration and used for nothing.
- **POU2F3** (Gerbe et al. 2016, Nature 529:226) — tuft cells: sourced factor, unsourced quantities.
- **HES1** (Jensen et al. 2000, Nat Genet 24:36) — the phenotype is an increase.
- **NOX1** (Coant et al. 2010, Mol Cell Biol 30:2636) — colon-specific and well sourced, but the
  phenotype is a conversion towards goblet cells at the cost of colonocytes; gating the colonocyte
  branch on it would predict a loss with no compensating gain.
- **DLL1 with DLL4** (Pellegrinet et al. 2011) — single Dll1 loss is an increase (1.83-fold goblet),
  and the double loss is the RBPJ phenotype already gated at C2.

## Which numbers are fitted, and to what

- **Fitted (one knob):** `surface_residence_days`, shared by colonocytes, goblet and enteroendocrine
  cells, fitted so that every compartment's pool sums to the published 6.0e10. The published mean
  residence is 3.4 d and the two measured numbers imply 3.53 d on their own.
- **Measured, used as given:** 6.0e10 cells, 1.7e10 replaced per day, lifespan 3.4 d (Sender & Milo
  2021); 20% goblet share of colonic epithelial cells (Katz et al. 2002); ~20 Lgr5+ cells per colonic
  crypt (Sasaki et al. 2016).
- **Inferred, not fitted:** that the measured 20% *standing* share can be used as the goblet share of
  daily *output* — which holds because all three mature residences are the one shared number; the
  enteroendocrine output share 0.004, inside the <1% bound; the colonocyte share as the remainder;
  the stem division interval (1 d), the transit-amplifying transit (2.5 d) and divisions (5), and the
  secretory progenitor residence (1 d), none with a primary colon citation fetched here.
- **Not claimed:** the colonic stem pool and the crypt count it implies; and the 20% goblet share,
  which comes back **by construction** and is not a check.

## Registered in advance about the formalism

R1 the two contested gates cannot be settled by a run, so the result must count agreement with a
named paper and never "seven of seven". R2 **a number, registered before the run**: blocking the
goblet route (gfi1, sox9, klf4) will leave the enteroendocrine pool at about **51 times** its
wild-type value, because the secretory progenitor backs up until its one open outflow carries the
whole inflow: 1 / (0.004 / 0.204) = 51.0. R3 without ATOH1 the colonocytes will stand *above* the
wild type (direction only; the transit-amplifying compartment carries a growth term and a cap).
R4 RBPJ's published conversion into goblet cells will read as a deletion of every lineage. R5 three
of the seven phenotypes are reductions and are stated as deletions, which is worse here than in the
small intestine because the lineage is 20% of the epithelium and not 4%. R6 three mature lineages,
not four, is the whole structural difference between the two programs once the villus is gone.

## Wiring into `body.bio` / `tissues.bio`: not wired, and why

Not wired, for exactly the reason lane-organism did not wire the small intestine. `tissues.bio` says
of its one gut population: `cell_type GutEpithelium { ... evidence: curated "Cell Ontology; grouped
from Small intestine epithelia, Colon epithelia, Stomach epithelia" }`, 1.97e11 cells and 4.03e10
replaced per day. The colon's 6.0e10 cells are already inside that number, so importing this module
beside it would count the colonic epithelium twice. It stands beside the body program as
`haematopoiesis.bio` and `intestine.bio` do.

---

# The colon's result (2026-10-10, lane-colon, `cf9c09a` registration, `3750411` module)

`data/organisms/human/colon.bio`, generated by `genomeos/organism/colon.py`, with its knockouts in
`colon_mutants.bio` and its tests in `tests/test_colon.py`.

## The wild type: one fit closing, not a prediction

| | module | published | kind |
| --- | --- | --- | --- |
| colonic epithelial cells | 5.9983e10 | 6.0e10 (Sender & Milo 2021) | the fit closing |
| cells replaced per day | 1.7e10 | 1.7e10 (Sender & Milo 2021) | the fit closing |
| surface residence | 2.57 d, **fitted** | 3.4 d published mean; 3.53 d is what the two numbers above imply | the one knob |
| goblet share of the mature cells | 0.200 | 0.20 of colonic epithelial cells (Katz et al. 2002) | **by construction**: the share was inferred as the output share and the three residences are one number. Not a check |
| goblet share of every epithelial cell | 0.146 | nothing: Katz counted total epithelial cells in sections | reported only |
| stem pool | 5.311e8 | nothing: it is a consequence of the flows | claimed of nothing |
| crypts implied at 20 Lgr5+ cells each | 2.66e7 | nothing was fetched to check it against | reported only |

Every compartment lands within 0.03% of its design (stem 5.311e8, transit-amplifying 1.226e10,
secretory progenitor 3.467e9, colonocytes 3.481e10, goblet 8.745e9, enteroendocrine 1.749e8), all
four `assert` lines pass and nothing is unknown.

## The knockouts: 7 registered, 7 holding the registered pattern, 5 counted as agreement

25 lineage checks, no exception, and no mutant number was fitted. The count is deliberately not
"7 of 7 reproduced": two gates are contested in the literature and the registration said so first.

| mutant | registered | run | counted as |
| --- | --- | --- | --- |
| TCF7L2 | no lineage at all | all three at 0 | **agreement with van Es et al. 2012**, disagreement with Angus-Hill et al. 2011 (contested in advance) |
| RBPJ | the proliferative compartment lost | all three at 0 | partial in advance: the published goblet conversion is not expressible |
| ATOH1 | goblet and enteroendocrine gone, colonocytes kept | exactly that | full |
| NEUROG3 | enteroendocrine alone gone | exactly that | full |
| GFI1 | goblet gone, colonocytes and enteroendocrine kept | exactly that | partial in advance: the published phenotype is a reduction to about a third, not a deletion |
| SOX9 | goblet gone, the other two kept | exactly that | partial in advance: "a strong reduction", not a deletion |
| KLF4 | goblet gone, the other two kept | exactly that | **agreement with Katz et al. 2002** (90% fewer, newborn), disagreement with Pellegrinet et al. 2011 (adult, unchanged) |

## The registered readings: all six held, and R2 held as a number

- **R2 held.** Blocking the goblet route leaves the enteroendocrine pool at **50.92** times the wild
  type against the **51.0** registered before the run, and the secretory progenitor backs up by the
  same factor (1.77e11 against 3.47e9). The mirror case confirms the mechanism: without NEUROG3 the
  goblet pool stands at 1.02 times the wild type, which is 1 / 0.9804.
- **R3 held in the direction registered.** Without ATOH1 the colonocytes stand at 1.77 times the wild
  type while Lo et al. 2016 report the absorptive lineage remaining, not expanding. Only the
  direction was registered, and that was right: the magnitude is not derivable from the shares.
- **R1, R4, R5, R6 held** as written, and are the reason the table above counts five agreements and
  not seven reproductions.

Worth saying plainly, because it is the difference from the small-intestinal result: **no registered
reading turned out wrong here, and that is not because this lane predicted better.** R2 and R3 are
precisely the two things lane-organism's registration got wrong and then pinned in
`tests/test_intestine.py` (`ad8c662`). This lane registered their consequence in advance because that
pin existed, and R2's 51.0 was arithmetic on the registered shares, not insight. A reading that can
be computed from the shares before the run is a weak prediction by construction.

## What stays open

- The stomach is still counts: with the colon modelled, `tissues.bio`'s `GutEpithelium` (1.97e11) now
  has mechanism for 1.2e11 + 6.0e10 = 1.8e11 of it, and the remaining 1.7e10 is the stomach.
- The Reg4+ deep crypt secretory cells — the colon's Paneth equivalent — wait on a share of output or
  a residence time, exactly as the tuft lineage does in the small intestine. The genetics exists
  (Sasaki et al. 2016) and the quantity does not.
- SPDEF waits on a colon sentence. Gregorieff et al. 2009 may well contain one; the publisher refused
  the fetch (HTTP 403) and a gate is not written on an abstract that says "intestinal epithelium".
- The two contested gates would need what the module has not got: an age axis (KLF4 is a newborn
  phenotype and an adult non-phenotype) and a graded gate (TCF7L2's two papers differ in how much of
  the factor is left). The second is the same limit the small-intestinal module and the mesoderm lane
  reported.
- A reduction still reads as a deletion, and in the colon that mis-states 20% of the epithelium.
- Not wired into `body.bio`, for the reason the small intestine is not: `tissues.bio` counts the
  colonic epithelium inside `GutEpithelium`, so importing the module beside it would count those
  6.0e10 cells twice.

`bio test data/demo genomeos/std data/organisms`: 54 files and 354/354 checks before, **56 files and
383/383 after** — 4 new checks in the module and 25 in the mutants.

---

# The stomach, written before anything was built (2026-10-10, lane-stomach)

The line the colonic section left open — "the stomach is still counts: with the colon modelled,
`tissues.bio`'s `GutEpithelium` (1.97e11) now has mechanism for 1.8e11 of it, and the remaining
1.7e10 is the stomach" — is what this section closes. The authoritative, machine-readable form of
everything below is `REGISTRATION` in `scripts/stomach_registration.py`, committed before a single
number of the module was fitted, so that the closed list provably cannot grow. Every citation here was
fetched in this lane (PubMed E-utilities, PMC full text, EBI OLS4 for the Cell Ontology ids) and the
sentence each gate rests on is quoted in that file.

## The gastric gland's lineages, and how the stomach differs from both gut modules

The oxyntic (corpus) gland unit has four regions — pit, isthmus, neck, base — and the isthmus, not the
base, is where the division happens. Karam & Leblond's five-part series on the mouse corpus counted
and timed every type:

| | small intestine | colon | stomach (this module) |
| --- | --- | --- | --- |
| geometry | crypt and villus | crypt and flat surface | pit, isthmus, neck, base; the proliferative zone is in the MIDDLE and cells migrate in both directions from it |
| where the stem cell sits | crypt base | crypt base | the isthmus: "granule-free cells, which are devoid of secretory granules and have the highest labeling index, 32.4%" (Karam & Leblond 1993, Anat Rec 236:259) |
| mature lineages | 4: enterocyte, goblet, enteroendocrine, Paneth | 3: colonocyte, goblet, enteroendocrine | 4 in this module: pit (foveolar), parietal, chief (zymogenic), enteroendocrine, with the mucous neck cell an intermediate on the chief route |
| residences | one fitted villus residence plus two measured | ONE shared fitted residence | four residences a factor of 60 apart: pit 3.1 d and parietal 54 d both measured, endocrine 86.6 d from a fetched half-life, chief 190 d inferred |
| the amplifying step | transit-amplifying, 5 inferred divisions | the same | NONE: the division is in the isthmus itself and no division count per transit was fetched, so the root divides |
| counts (the body program's source) | 1.2e11 cells, 2.2e10/day, 5.3 d | 6.0e10, 1.7e10/day, 3.4 d | **1.7e10 cells, 1.3e9/day, 13.0 d** (Sender & Milo 2021, Nat Med 27:45) |
| the gates | lineage losses | lineage losses | mostly NOT lineage losses: see below |

Two differences decide the shape of this module. First, the stomach's lineages do not share a
residence, so its standing composition and its daily-output composition are not the same thing, and
they diverge by an order of magnitude. Second — and this is the finding that shaped the closed list —
**most well-sourced gastric transcription-factor phenotypes are maturation defects with the cell still
present, not lineage losses.** In the intestine and the colon, deleting the gate deleted the lineage.
In the stomach, deleting XBP1, HNF4A, SPDEF, FOXQ1 or KLF4 leaves the cell standing there, smaller or
emptier or more numerous. A module whose gates are presence-or-absence has nothing to say about that,
so those five factors are excluded, with the fetched sentence as the reason in each case.

## The closed, sourced list: five gates, each with a fetched STOMACH phenotype

A factor well sourced in the intestine or the colon is not a source for the stomach. Each gate below
rests on a sentence about the stomach, quoted in the registration.

| # | Factor | Step | Published stomach phenotype of its loss | Source | Registered as |
| --- | --- | --- | --- | --- | --- |
| S1 | RBPJ | the isthmus compartment | pan-Notch inhibition leaves "no Ki67+ cells in most glands" of corpus and antrum | Kim & Shivdasani 2011, J Exp Med 208:677 | **contested**: the same paper's genetic arm (Osr1-Cre; RBP-Jk Fl/Fl) leaves the units alive and merely "reduced by 47%" by competition |
| S2 | ESRRG | parietal output | "whole swaths of gastric units completely lacking PCs ... the absence of PCs from birth onward" | Adkins-Threats et al. 2024, Cell Stem Cell 31:886 | **full**; the same mice's steady-state increase in proliferation is named as not expressible |
| S3 | MIST1 (Bhlha15) | chief output | the mature chief cells are "uniformly structurally abnormal" and the immature transitional cells accumulate markedly; the neck cells differentiate normally | Ramsey et al. 2007, Development 134:211 | partial: the mature cell is present and spoiled, not absent |
| S4 | NEUROG3 | endocrine output | the glandular stomach's A, D and G cells are absent and the enterochromaffin cells "decreased dramatically", with intestinal metaplasia | Lee et al. 2002, Genes Dev 16:1488 | partial, and with NO kept lineage registered: the metaplasia means the paper does not report the others normal |
| S5 | ASCL1 | endocrine output | newborn nulls lose the neurogenin-3-dependent endocrine markers "with normal expression of enterochromaffin-like cells, mucous cells, chief cells, and parietal cells" | Hibdon et al. 2023, Am J Physiol Gastrointest Liver Physiol 325:G458 | partial: the ECL subtype is normal and this module has no subtypes; a newborn readout |

Excluded, each with the fetched reason — the list is closed and none of these enters the module:

- **XBP1** (Huh et al. 2010, Gastroenterology 139:2038) — "did not affect survival of ZCs but
  prevented formation of their structure". The paper denies outright the loss a gate would state.
- **HNF4A** (Moore et al. 2016, Am J Physiol Gastrointest Liver Physiol 311:G267) — reduced chief-cell
  size and ER, "phenocopying ... Xbp1 deletion", plus *increased* proliferation in the isthmal zone.
- **SPDEF** (Horst et al. 2010, J Biol Chem 285:35047) — a candidate the task named; "impaired terminal
  maturation of antral mucous gland cells ... reduced numbers of secretory granules", plus hyperplasia.
  In the colon SPDEF was excluded for want of a sentence; here the sentence exists and says the wrong
  kind of thing.
- **FOXQ1** (Verzi et al. 2008, Gastroenterology 135:591) — the best pit-cell candidate there was:
  "virtual absence of ... MUC5AC ... a paucity of foveolar cell secretory vesicles". The cell is
  present and empty.
- **KLF4** (Katz et al. 2005, Gastroenterology 128:935) — "increased proliferation and altered
  differentiation", no lineage missing. An increase, which the formalism cannot express.
- **ATOH1** — the secretory gate of BOTH gut modules, and no stomach loss-of-function sentence was
  found for it in this lane. Carrying the gut gate across would be an assumption, not a source: the
  ground on which the colonic registration excluded ASCL2.
- **SOX2** (Que et al. 2007, Development 134:2521) — the candidate the task named first. The fetched
  phenotype is the squamous forestomach boundary in hypomorphic embryos, not an adult glandular lineage.
- **GATA4** (DeLaForest et al. 2021, Cell Mol Gastroenterol Hepatol 12:1391) — developmental
  morphogenesis of the hindstomach's columnar epithelium, before the gland's lineages exist.
- **LGR5 with AXIN2** (Barker et al. 2010, Cell Stem Cell 6:25) — lineage tracing and organoids: a
  marker of the antral stem cell, not a loss of function.
- **POU2F3 / the caveolated (tuft) lineage** (Karam & Leblond 1993, Anat Rec 236:333) — "Little
  information is available on the dynamics of caveolated cells". The stomach says the quiet part out
  loud: no share of output, no residence, no lineage.

## Which numbers are fitted, and to what

- **Fitted (one knob):** `isthmus_residence_days`, the isthmal granule-free cell's turnover time,
  fitted so that every compartment's pool sums to the published 1.7e10. It is held against a bound
  fetched **before** the fit: "less than 2.6 days for granule-free and pre-pit cells" (Karam & Leblond
  1993, part I).
- **Measured, used as given:** 1.7e10 cells and 1.3e9 replaced per day (Sender & Milo 2021); per
  gastric unit in the mouse corpus, 37 pit cells, 12.6 neck cells, 67.5 zymogenic cells and 26.0
  parietal cells (Karam & Leblond 1993, parts II–IV); the pit residence 3.1 d and the parietal
  residence 54 d (parts II and IV); the gastric endocrine half-life of about 2 months (Hibdon et al.
  2023).
- **Inferred, not fitted:** that a mouse corpus composition can stand for a human stomach that also has
  an antrum; the neck residence 10.5 d, the midpoint of the fetched 1-to-2-week range; the endocrine
  residence 86.6 d = 60 / ln 2, the mean the fetched half-life implies under the exponential loss law
  the module applies; the chief residence 190 d, for which **no number was fetched in this lane**; the
  endocrine standing count of 1.0 cell per unit, inside the carried <1% bound; and that no
  transit-amplifying compartment is needed.
- **Not claimed:** the standing composition, which comes back **by construction** because it is what
  the module was built from; the published total, which is the fit target; the implied number of
  gastric units; and the published 13.0-day lifespan, which is 1.7e10 / 1.3e9 and not independent.

## Registered in advance about the formalism

R1 one gate full, three partial, one contested, so no run may report "five of five". R2 **a number**:
the one fitted number will come out at about **2.50 d** — arithmetic on the registered inputs, not
insight — and the content is whether it lands inside the 2.6-day bound fetched before the fit. R3 **a
number**: the pit cells will be about 22% of the standing epithelium and about 93% of the daily output,
the chief cells about 40% standing and 2.8% of output; in the colon the two compositions were identical
by construction. R4 blocking the chief route is not the colon's backing-up case: the neck compartment
has exactly one outflow, so nothing leaves and it accumulates for as long as the run lasts — only the
direction is registered, and the magnitude is a property of the run length, not of the stomach.
R5 **a number**: moving the unfetched chief residence from 50 d to 1000 d will move the fitted number
by less than 1% (2.504 → 2.497 d), so nothing in this module rests on it. R6 a knockout's effect on the
standing count is not its effect on turnover — parietal loss is 16% of the cells and 4% of the output,
endocrine loss 0.6% and 0.09%. R7 **the lineage that carries 93% of the stomach's turnover has no
gate**, because every fetched pit-cell factor spoils the cell instead of removing it. R8 there is no
corpus/antrum axis and no gland geometry: the G cell sits inside the one endocrine lineage.

## Wiring into `body.bio` / `tissues.bio`: not wired, and why

Not wired, for exactly the reason the small intestine and the colon are not. `tissues.bio` says of its
one gut population: `cell_type GutEpithelium { ... evidence: curated "Cell Ontology; grouped from Small
intestine epithelia, Colon epithelia, Stomach epithelia" }`, 1.97e11 cells and 4.03e10 replaced per
day. **`tissues.bio` counts these 1.7e10 stomach cells inside that number**, so importing this module
beside it would count them twice. It stands beside the body program as `haematopoiesis.bio`,
`intestine.bio` and `colon.bio` do.

---

# The stomach's result (2026-10-10, lane-stomach, `f345ffe` registration, `c0e0c85` module)

`data/organisms/human/stomach.bio`, generated by `genomeos/organism/stomach.py`, with its knockouts in
`stomach_mutants.bio` and its tests in `tests/test_stomach.py`.

## The wild type: one fit closing and one construction, neither a prediction

| | module | published | kind |
| --- | --- | --- | --- |
| stomach epithelial cells | 1.7013e10 | 1.7e10 (Sender & Milo 2021) | **the fit closing**: it is the target the one knob was solved against |
| cells replaced per day | 1.301e9 | 1.3e9 (Sender & Milo 2021) | **by construction**: the output was used as given and distributed by the derived shares. Not a check |
| isthmus residence | **2.4987 d, fitted** | less than 2.6 d for the isthmal granule-free cells (Karam & Leblond 1993, part I) | the one knob — and it lands **inside a bound fetched before the fit** |
| standing composition per gastric unit | 37 : 12.6 : 67.5 : 26.0 exactly | the measured mouse corpus unit (Karam & Leblond 1993) | **by construction**: it is what the module was built from. Not a check |
| pit share of the epithelium / of the output | 0.221 / 0.934 | nothing measured for the human stomach | registered in advance as R3, and held |
| chief share of the epithelium / of the output | 0.404 / 0.028 | nothing | registered in advance as R3, and held |
| implied gastric units | 1.017e8 | nothing was fetched to check it against | reported only |

Every compartment lands within 0.09% of its design (isthmus 3.251e9, neck 3.797e8, pit 3.766e9,
parietal 2.646e9, chief 6.869e9, endocrine 1.018e8), all five `assert` lines pass and nothing is
unknown. What is worth stating plainly: **two of the three headline numbers are not results.** The
total is the fit's own target and the daily output was used as given. The one number that could have
come out wrong is the fitted residence, and the reason it is worth anything is that the 2.6-day bound
was fetched and committed before the fit was run.

## The knockouts: 5 registered, 5 holding the registered pattern, 1 counted as a full agreement

22 lineage checks, no exception, and no mutant number was fitted. The count is deliberately not
"5 of 5 reproduced": one gate is registered full, three partial and one contested, and the registration
said so first.

| mutant | registered | run | counted as |
| --- | --- | --- | --- |
| RBPJ | no lineage at all | all four at 0 | **agreement with the pan-Notch (DBZ) arm of Kim & Shivdasani 2011**, disagreement with the same paper's genetic arm (contested in advance) |
| ESRRG | parietal gone, the other three kept | exactly that | **full agreement with Adkins-Threats et al. 2024**; their steady-state rise in proliferation is not expressible |
| MIST1 | chief gone, the other three kept | exactly that | partial in advance: Ramsey et al. 2007 report the mature chief cell present and spoiled, not absent |
| NEUROG3 | endocrine gone, nothing registered as kept | exactly that | partial in advance: the enterochromaffin subtype is a decrease, and the intestinal metaplasia is not expressible |
| ASCL1 | endocrine gone, the other three kept | exactly that | partial in advance: the enterochromaffin-like cells are normal and this module has no subtypes; a newborn readout |

## The registered readings: all eight held, and the three numbers held as numbers

- **R2 held.** The one fitted number came out at **2.4987 d** against the **2.50** registered before the
  run, and below the **2.6 d** bound fetched before the fit. The arithmetic was never in doubt; the
  comparison was.
- **R3 held as four numbers.** Pit 0.221 standing and 0.934 of the output, chief 0.404 and 0.028, each
  within 0.01 of the registered value. The stomach's standing share of a lineage says almost nothing
  about its share of the turnover, which the colonic module could not show because one shared residence
  made the two compositions identical by construction.
- **R4 held in the direction registered, and only that.** With the chief route gated off the neck
  compartment stands at **32.7 times** its wild-type pool, and no multiple was predicted: it has one
  outflow, so nothing leaves and the pool is the wild-type chief inflow times the time in the Steady
  stage. The colon's 51 was predictable because its progenitor had a second outflow; this is not that
  case, and the registration said so before the run.
- **R5 held.** Moving the one unfetched residence from 50 d to 1000 d moves the fit from 2.5043 to
  2.4969 d. Nothing in the module rests on it.
- **R1, R6, R7, R8 held** as written, and R1 is the reason the table above counts one full agreement and
  not five reproductions.

**What the registration did not foresee, recorded as the colon's and the small intestine's results
record theirs.** The root's self-renewal term was solved for the wild-type outflow, so gating off any
one route leaves it too large and the root rises until its cap stops it, carrying every surviving
lineage up with it: without ESRRG — 3.8% of the daily output — the other three lineages stand at about
2.5 times their wild-type pools, and without ASCL1 — 0.09% of the output — at about 1.1 times. Neither
multiple is a claim about the stomach. It is the formalism and the cap, and it is now pinned in
`tests/test_stomach.py` so the next module meets it as a test and not as a surprise.

## What stays open

- **The gut's grouped population is now fully covered, and that is a bookkeeping statement, not a
  biological one.** `tissues.bio`'s `GutEpithelium` (1.97e11) has mechanism modules beside it for
  1.2e11 + 6.0e10 + 1.7e10 = 1.97e11. Not one of the three is wired in, for the same reason: the
  grouped population already counts those cells.
- **The chief cell's residence was never fetched.** R5 shows nothing rests on it, which is not the same
  as having it. A measured zymogenic turnover time would turn an inference into a citation.
- **The pit lineage has no gate**, and it carries 93% of the stomach's daily output. Every fetched
  candidate (FOXQ1, SPDEF, KLF4) spoils the cell instead of removing it, so no knockout in this module
  touches the lineage that matters most to the turnover.
- **A maturation defect is not a deletion**, and in the stomach that mis-states more than in either gut
  module: the whole XBP1 → MIST1 → chief-cell axis is about a cell that stays where it is and stops
  working. A formalism with a maturity axis, not only a presence axis, is what those five exclusions
  are waiting for.
- **No corpus/antrum axis.** The parietal and chief lineages are corpus-only in life; the G cell is
  antral; this module has one gland and folds them together. The second knockout of the endocrine step
  (ASCL1, antral) and the first (NEUROG3, whose fetched sentence names the G cell) are therefore
  indistinguishable here, and they ran identically.
- **The caveolated (tuft) lineage** waits on a quantity, and the fetched sentence says outright that
  none exists: "Little information is available on the dynamics of caveolated cells."

`bio test data/demo genomeos/std data/organisms`: 56 files and 383/383 checks before, **58 files and
410/410 after** — 5 new checks in the module and 22 in the mutants.

# The epidermis, chosen and registered before anything was built (2026-10-10, lane-module3)

Haematopoiesis, the small intestine, the colon and the stomach leave the human body program with
mechanism in blood and the gut and counts everywhere else. This section chooses the next tissue by
criteria written down first, says why the others were set aside, and registers the closed list, the
knockouts, the one fitted number and eight readings before a line of the module exists. The
authoritative, machine-readable form is `REGISTRATION` in `scripts/epidermis_registration.py`, committed
with this section. Every citation was fetched in this lane (PubMed E-utilities, PMC full text, EBI OLS4
for the Cell Ontology ids, the `milo-lab/cellular_turnover` repository for the tabulated kinetics,
BioNumbers for Hoath & Leahy's Table I), and the sentence each gate rests on is quoted in that file.

Stop list (ROADMAP §5, the list of 2026-09-27 reinstated): a human tissue module is not the worm's fate
ceiling, adds no reader biosample, builds nothing in area J and nothing on the node model.

## The criteria, written before the candidates were compared

| # | Criterion | Why |
| --- | --- | --- |
| C1 | a tissue that renews in the adult, with a published lineage hierarchy | a flow module is a renewing hierarchy |
| C2 | already a counts-only population of `tissues.bio`, with a human count and a daily turnover that Sender & Milo 2021 did NOT extrapolate from rodents | so the one fitted number closes against a published human number, as in every module before it |
| C3 | per-compartment kinetics, measured in human, for at least two compartments | so the fitted number is one knob and not the whole module |
| C4 | at least six published single-gene perturbations with stated phenotypes in that tissue, fetched and read, most of them a LOSS of a compartment | the stomach showed that a list of maturation defects leaves the lineages ungated |
| C5 | outside the gut, and best of all in a germ layer no module covers | three gut modules already exist; blood is mesoderm, the gut endoderm |
| C6 | not on the stop list | ROADMAP §5 |

## The candidates

| Candidate | C1 | C2 | C3 | C4 | C5 | Verdict |
| --- | --- | --- | --- | --- | --- | --- |
| **interfollicular epidermis** | basal, spinous, granular, cornified (Blanpain et al. 2006: "a terminal differentiation program consisting of three distinct stages: spinous, granular layer, and stratum corneum") | `Keratinocyte`, 1.2e11 cells, 3.7e9/day, 31 d, not extrapolated | three compartments with turnovers, Tc 311 h and a 39-day total (Weinstein et al. 1984); stratum corneum 14 d (Hoath & Leahy 2003) | 24 perturbations of 16 genes fetched; 7 genes enter as gates, each with a removal-type phenotype in at least one fetched experiment | ectoderm, the first | **chosen** |
| spermatogenesis | yes | no: no row in the Sender & Milo table and no population in `tissues.bio` | not assessed | not assessed | yes | set aside on C2 |
| bronchial epithelium | yes | 1.04e10 cells, 139 d, not extrapolated | only a whole-population lifespan was looked at; no per-compartment human kinetics fetched | not assessed, nothing fetched | endoderm-derived, outside the gut | set aside on C3 |
| hepatocytes, alveolar epithelium | slow or facultative | no: turnover extrapolated from rodents | — | — | — | set aside on C2 |

Nothing about the genetics of the set-aside candidates is claimed: none was fetched.

## The compartments and the numbers they rest on

Four populations in a chain: `BasalKeratinocyte` (CL:0002187) → `SpinousKeratinocyte` (CL:0000649) →
`GranularKeratinocyte` (CL:0002189) → `Corneocyte` (CL:0002153), shed from the surface. Sender & Milo's
1.2e11 are the NUCLEATED cells ("As corneocytes are essentially dead cells, we don't count them",
`Skin_cells.ipynb`), so the closure is on the first three, and the corneocytes are a fourth population
the count leaves out.

- **Measured, used as given:** 1.2e11 nucleated cells and 3.7e9 per day (Sender & Milo 2021); the basal
  compartment's 311 h, Weinstein et al. 1984's cell-cycle time for the proliferative compartment, which
  Sender & Milo tabulate as the germinative layer's 12.958 d; the stratum corneum's 14 d (Hoath & Leahy
  2003, Table I, from Rothberg 1961, Weinstein & Van Scott 1965, Jansen 1974); one flux through every
  layer, on Weinstein's "very similar 'turnover' rates in each of the compartments", 1246, 1417 and 1490
  cells/day/mm².
- **Fitted, the one knob:** the suprabasal residence (spinous plus granular), so the three nucleated pools
  sum to 1.2e11. R2 states its value, 19.47 d, in advance and the published span it must land in.
- **Inferred, no primary citation:** the granular residence, 2.0 d, which only splits the fitted suprabasal
  pool between two layers (R5); no death in the viable layers; a cap on basal renewal near its pool (R6).
- **Not claimed:** the basal pool, about 4.8e10, checked against no basal density.

## The closed, sourced list: seven gates

| # | Factor | Gates | Published epidermal phenotype of its loss | Source | Registered as |
| --- | --- | --- | --- | --- | --- |
| E1 | TP63 | the basal compartment | the skin "lacks stratification and does not express differentiation markers"; "a striking absence of all squamous epithelia" | Mills et al. 1999, Nature 398:708; Yang et al. 1999, Nature 398:714 | full |
| E2 | DNMT1 | basal divisions only | in regenerated human epidermis, "exit from the progenitor cell compartment, premature differentiation and eventual tissue loss" | Sen et al. 2010, Nature 463:563 | partial: the premature differentiation and the loss of a formed tissue are not expressible; the knockdown left 34% of tissues alive at week 3 |
| E3 | RBPJ | the spinous compartment | "Spinous layers are absent"; "fewer granular layer cells"; basal markers "relatively normal" | Blanpain et al. 2006, Genes Dev 20:3022 | partial: "fewer granular" is overstated as none by a chain |
| E4 | IRF6 | the granular compartment | "the basal layer was present along with a greatly expanded spinous layer, but the outer two layers appeared absent" | Ingraham et al. 2006, Nat Genet 38:1335 (Richardson et al. 2006 agreeing) | full; the spinous proliferation named as not expressible |
| E5 | KDF1 | the granular compartment | granular and cornified markers "are not induced"; basal and spinous markers expanded | Lee, Kong & Weatherbee 2013, Dev Biol 383:201 | full |
| E6 | CHUK (IKKα) | the granular compartment | germline null: no cornified layer, "granular cells ... were not present", basal unchanged, K10+ cells "increased dramatically" | Li et al. 1999, Genes Dev 13:1322 | **contested**: keratinocyte-restricted ablation gives "a normally differentiated stratified epidermis" (Gareus et al. 2007, Nat Cell Biol 9:461) |
| E7 | ZNF750 | the granular compartment | human regenerated epidermis: K1 and K10 "unaffected", late differentiation genes (FLG, LOR) not induced | Sen et al. 2012, Dev Cell 22:669 | **contested**: the mouse Zfp750 null keeps its granular layer and stratum corneum, loricrin "unchanged" (Butera et al. 2023, Sci Adv 9:eadg5423) |

The granular-to-cornified step has no gate. Every fetched candidate for it spoils the layer rather than
removing it, the stomach's finding in a new tissue. Excluded, each with the fetched reason in the
registration: **NOTCH1** (hyperplasia, "an increase rather than decrease of granular marker expression
and no noticeable defect in spinous layers"), **ADAM10** (a reduction with premature differentiation),
**RIPK4** (the full knockout's layer phenotype is in no sentence fetched here, and the epidermis-specific
knockout "was normally stratified"), **KLF4** (barrier lost "without morphological and biochemical
alterations to the well-known structural features of epidermis"), **GRHL3**, **TGM1** and **PRDM1** (a
present but defective stratum corneum or granular layer), **SFN** (the Er allele makes a truncated protein
with a heterozygous phenotype, so it is not a loss of function), **OVOL1** (an increase in
proliferation), and three PAIRS considered without a fetch (CEBPA/B, MAF/MAFB, HDAC1/2), because a gate
of one factor cannot state redundancy.

## The knockouts, registered with their pass rules

Seven `experiment` blocks in `data/organisms/human/epidermis_mutants.bio`, each run to one year: lost
layers exactly 0, kept layers at least 1e6, and for the four granular gates the registered direction that
the spinous pool stands above twice its wild-type design. TP63 leaves the single seed cell in the basal
compartment and nothing above it; DNMT1 leaves every layer below one cell. Three are full (TP63, IRF6,
KDF1), two partial (DNMT1, RBPJ), two contested (CHUK, ZNF750), so the result cannot read "seven of
seven" (R1).

## Registered in advance about the formalism and the sources

- **R2** the fitted suprabasal residence will be 19.47 d (1.2e11 / 3.7e9 − 311 / 24): between Weinstein's
  12 d and Halprin's 26–42 d, inside neither study.
- **R3** the corneocyte pool, not fitted, will be 5.18e10, 1.53 times the 3.38e10 that 17,778 per mm²
  over 1.9 m² gives: exactly the ratio of Sender & Milo's flux to Hoath & Leahy's rate of stratum corneum
  formation, a disagreement between the sources that the module inherits.
- **R4** a blocked granular route leaves the spinous compartment with no outflow, so it grows at the
  wild-type inflow for the whole run: only the direction is claimed, the size (about 18 times) belongs to
  the run's length. For ZNF750 the same growth contradicts the human paper's "unaffected" K1 and K10.
- **R5** the inferred granular residence moves only the spinous/granular split.
- **R6** under RBPJ the basal compartment stands within 10% of its wild-type pool, because of a cap this
  module assumes; the sign agrees with the paper and is not evidence.
- **R7** TP63 and DNMT1 cannot be told apart by what is left. **R8** no appendages, no space, no barrier.

## Wiring into `body.bio` / `tissues.bio`: not wired, and why

`tissues.bio`'s `Keratinocyte` already counts these 1.2e11 cells; importing the module beside it would
count the epidermis twice, which is why no mechanism module is wired.
