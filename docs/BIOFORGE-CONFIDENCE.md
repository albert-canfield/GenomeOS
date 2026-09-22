# BioForge's predicted confidence: what it is, and why it cannot be calibrated

*2026-09-22. Area H, third open item: "calibrate predicted confidence: how often a design's answer
matches the published outcome". Registration and verdict in `genomeos/forge/calibration.py`; the
facts the verdict rests on are pinned in `tests/test_forge_calibration.py`.*

**The item asked the wrong question, and the answer is worth having anyway.** It presumes the
number beside a design's answer is a probability whose level could be checked. It is not one. It
is a single hand-set constant, and a constant has one reliability bin at every sample size, so no
calibration curve over it can exist at any n. That conclusion is read off the source rather than
off the data, which makes it the rare small-n result that does not depend on n.

## 1. What the confidence is

`genomeos/organism/forge.py:143` writes

```python
# fmt: off
confidence=0.3,
```

as a literal onto the `Experiment` that a `design` block emits. Six further sites cap an existing
confidence at the same value — `organism/forge.py:66,70`, `forge/design.py:59,63,66`, and
`attribution/budget.py:60` for an entirely unrelated quantity in another area. That last one is
the tell: a number that means the same thing for a design search over a worm embryo and for a
phyloP budget over chromosome 21 is not a probability about either. It is a provenance marker
reading *predicted, not validated*.

The caps are defensible. `min(x, 0.3)` can only lower a number, so it cannot inflate a claim. The
literal at line 143 is different in kind, because it creates a number from nothing and attaches it
to a freshly predicted entity, which is where a reader is most likely to take it for a posterior.

**The search computes better quantities and discards all of them.** `run_design` knows each
candidate's `loss`, whether it is `feasible`, how many `evaluations` it took, how many rival
candidates were also feasible, and the margin to the runner-up. None reaches the confidence. Over
the four designs the repository ships, feasible-candidate counts run 2, 6, 2, 22 and evaluations
run 7, 6, 11, 29 — and every one of the four answers is stamped 0.3. The design block's *own*
stated confidence (0.5 and 0.6 in the shipped programs) is parsed, stored on the `Design`, and
then ignored by `to_experiment()`.

So: **a hand-set constant, not a fitted quantity and not an aggregate of other confidences.** The
first item of the task answers itself — this was never a calibration question.

## 2. The outcome census, counted before anything was designed

**n = 2.** Two `design` blocks in this repository state a target and record a published,
citation-backed answer, both in `data/organisms/celegans/designs.bio`: `two_intestinal_founders`
(Lin et al. 1995 → POP-1) and `two_ems_founders_no_germline` (Mello et al. 1992 → PIE-1).

**n = 4** if two more are admitted whose answers are recorded only in prose, without a citation in
the program file: `lose_neutrophils_keep_the_rest` (→ GFI1) and `no_lymphocytes_at_all` (→ IKZF1)
in `data/organisms/human/haematopoiesis_designs.bio`. The two populations are reported separately
and never pooled, per the 2026-09-17 entry in LESSONS.md.

| quantity | count |
|---|---|
| `design` blocks in `data/` | 4 |
| of those, published outcome cited in the file | 2 |
| of those, answer recorded only in prose | 2 |
| `design` blocks in tests, as toys with no published answer | 2 |
| JSON/TSV/markdown fixtures pairing a design with a published outcome | **0** |
| tests that read a published outcome and compare it to a design answer | **0** |

An adjacent and much larger population exists and does not substitute: 31 `experiment` blocks in
`data/`, 29 of them with a published citation, plus 58 `# test:` lines across 21 of 65 `.bio`
files. Those are perturbation/outcome pairs, not design/outcome pairs. An experiment *states* the
perturbation and checks that the simulator reproduces the phenotype; a design must *find* the
perturbation. There is nothing in an experiment for BioForge to get wrong, so none of the 29 can
be borrowed to fill the missing population.

## 3. The registration

Written after the census and after the degeneracy was read off the source, and before any match
rate was compared with 0.3 or any interval computed. Full text in
`genomeos/forge/calibration.py:PREREGISTRATION`.

- **Compared:** for each design with a published outcome, the knockout set BioForge returns as its
  best feasible candidate against the one the cited publication reports. Exact set equality is a
  match; anything else is a miss.
- **Baselines.** *Trivial/modal:* always predict the modal outcome. Every design in the repository
  is solved, so the modal outcome is "the answer is right", and a predictor that says so
  unconditionally scores exactly what BioForge scores. Registered as the baseline BioForge must
  **beat, not tie**. *Uniform from the candidate list:* 1/6, 1/5, 1/10 and 1/28 per design; all
  four by chance is 1/8400. *Constant confidence:* quote 0.3 for everything, which is what the code
  does.
- **Calibrated would require all three:** at least two distinct confidence values so that at least
  two bins are populated; each populated bin's observed rate interval containing its stated value;
  and bin order agreeing with observed-rate order.
- **Uncheckable if any one of:** (a) the predictor is single-valued — one bin, no curve, at any n;
  (b) the outcome is single-valued — no negative case, no rate; (c) the 95% Clopper-Pearson
  interval is wider than 0.5 at the available n; (d) the candidate list or the module rules were
  authored from the publication the answer is scored against.
- **Population clause.** The 0.3 is quoted **unconditionally, for every feasible answer BioForge
  emits, over any organism program, for any target** — nothing in the code narrows it. It could be
  checked against 2 designs, or 4. Those are not the same population: the second is a hand-built
  subset chosen by an author who held the published answer while writing both the candidate list
  and the module rule leading to it.
- **Registered expectation:** refusal.

## 4. The result: uncheckable by construction, on four counts before n

**(a) The predictor is degenerate.** One distinct value over every design, and by inspection of
line 143 over every design that could ever be written. A reliability diagram needs two populated
bins. This is fatal at every sample size and is the whole answer to the item as posed.

**(b) The outcome is degenerate too.** All four designs solve, at loss exactly 0.0. The label
column is constant. A 2×2 with both margins collapsed supports no estimate.

**(c) The interval, if one insisted.** At 4 of 4 the exact 95% interval is **[0.3976, 1.0000]** and
**excludes 0.3** — the arithmetic says the stated confidence is too low and invites raising it. At
2 of 2 it is [0.1581, 1.0000] and includes 0.3. The same quantity therefore rejects or accepts
depending on which of the two populations is used, which is the 2026-09-17 finding arriving again
on new material. The narrowest interval available at n = 4 is 0.6024 wide.

**That "too low" number is not published, and (d) is why.** The population it is computed on is
the population into which the answer was placed before the search ran:

> `designs.bio:14` states its published outcome as *"Lin et al. 1995: without POP-1, MS takes the
> E fate"*. The rule the search must traverse to reach that outcome, `fate_E` at
> `data/organisms/celegans/founders.bio:47`, carries the evidence *"Lin et al. 1995 (without POP-1,
> MS takes the E fate)"* at confidence 0.9.

The same clause, of the same paper, on both sides of the comparison — once as the question's
answer key and once as a hand-written rule inside the model. GFI1 is the same story:
`GMP_to_Myeloblast` requires GFI1 and cites Hock et al. 2003 (GFI1: neutropenia), which is the
phenotype the design asks for. And every `knockout_any_of` contains the published factor, because
the author put it there; no design whose answer is absent from its own candidate list exists here.

**BioForge is not predicting POP-1. It is evaluating a rule authored from POP-1, over a list
containing POP-1.** The 4 of 4 measures that the forward simulator reads its own rules correctly —
which is worth knowing, and is what `tests/test_experiment.py` and `tests/test_design.py` already
establish. It carries no information about a design question whose answer was not written into the
model first.

**Verdict: `UNCHECKABLE_BY_CONSTRUCTION`.** The number has never been a probability. It should be
read, and documented, as a provenance marker.

## 5. What would change it, smallest obstruction last

1. **A negative case.** Every design solves and every answer is right, so the outcome column is
   constant. At least one design whose published answer is "no single perturbation achieves this",
   or whose candidate list does not contain the answer, is needed before the question has two sides.
2. **A blind candidate list.** Fix `knockout_any_of` by a rule written before the design is chosen —
   every factor the module declares, say — rather than by curation around the answer. This is the
   loci-benchmark move: the frame is a function of the data, not of the result.
3. **A design whose answer is not already in the module.** While the rule that yields the outcome
   cites the paper the outcome comes from, the search retrieves. Only a held-out publication, used
   to author no rule, tests prediction.
4. **A confidence that varies.** Derive it from something the search already computes and throws
   away — feasible candidates found, or the margin to the runner-up. Necessary for calibration and
   nowhere near sufficient: the derived number would then need calibrating on a population that
   does not yet exist. Item 1 comes first.
5. **n.** With a varying predictor and both outcomes present, telling a stated 0.3 from a true 0.9
   by exact binomial test at 80% power needs n ≥ 6. The repository holds 2 to 4. **n is the
   smallest of the five obstructions, and the only one that more of the same work would fix** —
   which is why "collect more designs" would have been the wrong response to this item.

## 6. What should change in the code, and what should not

Nothing measured here licenses moving the 0.3 up or down: no honest number exists to move it to.
Two changes are honest and neither needs data.

- **Say what it is where it is read.** `format()` prints the constant directly beneath a truthful,
  loss-based feasibility verdict it bears no relation to. The word *confidence* there should be
  labelled as provenance, not as a rate.
- **Stop discarding the design's own stated confidence, or stop parsing it.** At present a program
  author writes `confidence: 0.6`, the parser stores it, and `to_experiment()` silently overwrites
  it with 0.3. One of those two behaviours is wrong and they should not both remain.

Both are edits to how the number is *described*, not to its value, and are left to whoever takes
area H next so that this entry records a measurement rather than a change.
