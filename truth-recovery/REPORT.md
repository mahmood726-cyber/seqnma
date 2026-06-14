# SeqNMA — Truth-Recovery Report

**Verdict: STRONG VALIDATION** (one minor, honest small-k note — not a bug)

## Method under test
SeqNMA is a browser-based **sequential network meta-analysis** engine
(`index.html`, single-file). The pure numerical core:

- `runNMA` — contrast-based NMA by weighted least squares with a network
  DerSimonian-Laird heterogeneity estimate, then a **random-effects
  refit** with weights `1/(v + tau2)`. CIs use the RE variance-covariance
  `XtWX_inv`, as `effect ± 1.96·se`.
- `computeRIS` — Required Information Size with heterogeneity design
  effect `D^2 = 1 + tau2·(Σ(1/v_i^2)/(Σ(1/v_i))^2·k − 1)`.
- `obrienFlemingBoundary` — `z_k = z_{adjAlpha/2}/sqrt(t)`.
- Bonferroni across `C = T(T-1)/2`: `adjAlpha = alpha/C`.

## Harness (additive only — app code untouched)
- `engine.mjs` — functions extracted **verbatim** from `index.html` (+ the
  matrix utilities they call). DOM-free.
- `dgp-nma.mjs` — seeded (mulberry32) known-truth NMA DGP. Connected
  4-treatment network (ref A), **known basic params**
  `d={A:0,B:0.40,C:−0.25,D:0.65}` so every true relative effect `d_i−d_j`
  is known, **known common tau2**. Contrasts drawn from
  `y ~ N(d_t2−d_t1, tau2+se^2)` — the exact RE contrast model.
- `harness.mjs` — runs the engine on the DGP; measures coverage of the
  true relative effects, bias, tau2-hat, and the boundary identities.

## Results (measured; 2000 reps/cell, coverage over all 6 contrasts, nominal 0.95)

| true tau2 | studies/edge | N studies | overall coverage | mean |bias| | mean tau2-hat |
|----------:|-------------:|----------:|-----------------:|-----------:|--------------:|
| 0.00 | 3 | 18 | 0.962 | 0.0027 | 0.0078 |
| 0.00 | 8 | 48 | 0.960 | 0.0010 | 0.0043 |
| 0.05 | 3 | 18 | 0.925 | 0.0046 | 0.0509 |
| 0.05 | 8 | 48 | 0.938 | 0.0011 | 0.0496 |
| 0.15 | 3 | 18 | 0.926 | 0.0018 | 0.1512 |
| 0.15 | 8 | 48 | 0.943 | 0.0037 | 0.1485 |

Smallest network (k=12, 2 studies/edge, 4000 reps): coverage **0.911**,
mean tau2-hat 0.0518 (true 0.05).

### O'Brien-Fleming + Bonferroni boundary (deterministic)
- alpha=0.05, T=4 ⇒ C=6, adjAlpha=0.008333; z_crit=qnorm(1−adjAlpha/2)=2.6387.
- boundary at t=1 = **2.6387** (= z_crit). boundary at t=.25 = **5.2773** = 2·z_crit.
- OBF ratio b(.25)/b(1) = **2.0000** = 1/sqrt(.25). All match closed form <1e-9.

## Findings
1. **CIs are genuinely random-effects — the LivingNMA fixed-effect-CI bug
   is ABSENT.** Under tau2=0.15 coverage stays 0.926–0.943; it does NOT
   collapse toward ~0.6. `runNMA` refits with `1/(v+tau2)` and reports RE CIs.
2. **Network DL tau2 well recovered**: 0.051/0.050 at true 0.05;
   0.151/0.149 at true 0.15. Point estimates unbiased (|bias|<0.005).
3. **Sequential machinery exactly correct**: OBF `z/sqrt(t)` shape and
   Bonferroni `alpha/C` match closed form to <1e-9.
4. **Minor honest small-k note (NOT a bug):** at k≈12 the 95% CIs
   under-cover slightly (~0.91), the standard normal-quantile limitation
   (z=1.96 with DL-estimated tau2 at small k). tau2-hat itself is unbiased,
   so the cause is the quantile choice, not the estimator. *Informational*
   refit (not shipped into the app): replacing 1.96 with t_{df}, df=n−(T−1),
   lifts the worst cell **0.911→0.945**. Coverage already returns to ~0.94 by k=48.

## Recommendation
Ship as-is — estimator and sequential boundaries are correct and well
calibrated. Optional enhancement: a `t_{n−(T−1)}` (or HKSJ-style) CI for
small sparse networks to close the ~0.91→0.95 gap.
