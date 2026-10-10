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
