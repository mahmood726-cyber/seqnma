// node --test  — asserts the MEASURED invariants of SeqNMA's engine.
// No hand-entered statistics: every threshold is checked against a
// fresh harness run on known-truth data.

import { test } from 'node:test';
import assert from 'node:assert/strict';
import { measureCoverage, measureBoundary } from './harness.mjs';
import { runNMA } from './engine.mjs';
import { makeRng, generateNetwork, trueTheta } from './dgp-nma.mjs';

// 1. RE-NMA CIs cover true relative effects near nominal at zero
//    heterogeneity (sanity: the WLS engine + DL is unbiased).
test('coverage ~ nominal at tau2=0 (no heterogeneity)', () => {
  const r = measureCoverage(0, 6, 2000, 101);
  assert.ok(r.overallCoverage > 0.93 && r.overallCoverage < 0.97,
    `tau2=0 coverage ${r.overallCoverage} should be ~0.95`);
  assert.ok(r.meanAbsBias < 0.01, `bias ${r.meanAbsBias} should be ~0`);
});

// 2. THE LivingNMA BUG CHECK: under real heterogeneity the CIs must NOT
//    collapse. A fixed-effect-CI bug would push coverage toward ~0.6.
//    The engine uses 1/(v+tau2) weights, so coverage must stay >= 0.90.
test('CIs are random-effects: NO coverage collapse under heterogeneity', () => {
  const r = measureCoverage(0.15, 8, 2000, 202);
  assert.ok(r.overallCoverage >= 0.90,
    `tau2=0.15 coverage ${r.overallCoverage} — collapse below 0.90 would indicate FE-CI bug`);
});

// 3. Network DerSimonian-Laird tau2 is recovered (unbiased to within 20%).
test('network DL tau2 is recovered', () => {
  const r1 = measureCoverage(0.05, 8, 2500, 303);
  assert.ok(Math.abs(r1.meanTau2Hat - 0.05) < 0.01,
    `tau2-hat ${r1.meanTau2Hat} should be ~0.05`);
  const r2 = measureCoverage(0.15, 8, 2500, 404);
  assert.ok(Math.abs(r2.meanTau2Hat - 0.15) < 0.03,
    `tau2-hat ${r2.meanTau2Hat} should be ~0.15`);
});

// 4. O'Brien-Fleming + Bonferroni boundary formula is exact.
test('O\'Brien-Fleming boundary = z_{adjAlpha/2}/sqrt(t), Bonferroni adjAlpha=alpha/C', () => {
  const b = measureBoundary(0.05, 4);
  assert.equal(b.C, 6);                                   // T(T-1)/2
  assert.ok(Math.abs(b.adjAlpha - 0.05 / 6) < 1e-12);
  assert.ok(Math.abs(b.boundaryAt1 - b.zCrit) < 1e-9,
    'boundary at t=1 must equal critical value');
  assert.ok(Math.abs(b.obfRatio - 2.0) < 1e-6,
    'b(.25)/b(1) must equal 1/sqrt(.25)=2 (OBF shape)');
});

// 5. Direct contrast: a 2-arm-only subnetwork recovers the known theta.
test('engine recovers a known direct contrast point estimate', () => {
  const rng = makeRng(777);
  // many studies on a single edge -> point estimate -> true theta
  let est = 0, n = 0;
  for (let i = 0; i < 1500; i++) {
    const { studies, treatments } = generateNetwork(0.05, 8, rng);
    const fit = runNMA(studies, treatments);
    if (!fit) continue;
    const ab = fit.comparisons.find(c => c.comp === 'A vs B');
    est += ab.effect; n++;
  }
  est /= n;
  assert.ok(Math.abs(est - trueTheta('A', 'B')) < 0.02,
    `mean A-vs-B estimate ${est} should be ~${trueTheta('A', 'B')}`);
});
