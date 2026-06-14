// ============================================================
// harness.mjs — wire SeqNMA's OWN runNMA engine to the known-truth
// NMA DGP and measure the operating characteristics the method targets:
//   (1) coverage of the TRUE relative effects d_i - d_j by the engine's
//       95% CIs (the engine uses 1.96*se), across tau2 levels and
//       studies-per-edge counts;
//   (2) bias of the estimated relative effects;
//   (3) mean estimated network tau2 vs the true tau2 (DL calibration);
//   (4) a deterministic check of the O'Brien-Fleming boundary +
//       Bonferroni adjustment (z_k = z_{alpha/2C}/sqrt(t)).
//
// All numbers are MEASURED by running the engine. No hand-entered
// expected values.
// ============================================================

import { runNMA, obrienFlemingBoundary, qnorm } from './engine.mjs';
import { makeRng, generateNetwork, allTrueContrasts, TRUE_D } from './dgp-nma.mjs';

// Monte-Carlo coverage/bias of the engine's relative-effect CIs.
export function measureCoverage(tau2, studiesPerEdge, nRep, seed) {
  const contrasts = allTrueContrasts();
  const cov = {};         // comp -> count covered
  const biasSum = {};     // comp -> sum(est - true)
  const tau2Sum = { v: 0 };
  let nUsed = 0;
  for (const c of contrasts) { cov[c.comp] = 0; biasSum[c.comp] = 0; }

  const rng = makeRng(seed);
  for (let rep = 0; rep < nRep; rep++) {
    const { studies, treatments } = generateNetwork(tau2, studiesPerEdge, rng);
    const fit = runNMA(studies, treatments);
    if (!fit) continue;
    nUsed++;
    tau2Sum.v += fit.tau2;
    for (const c of contrasts) {
      const est = fit.comparisons.find(x => x.comp === c.comp);
      if (!est) continue;
      if (c.theta >= est.ci_lo && c.theta <= est.ci_hi) cov[c.comp]++;
      biasSum[c.comp] += (est.effect - c.theta);
    }
  }

  const perComp = contrasts.map(c => ({
    comp: c.comp,
    theta: c.theta,
    coverage: cov[c.comp] / nUsed,
    bias: biasSum[c.comp] / nUsed,
  }));
  const overallCov = perComp.reduce((a, b) => a + b.coverage, 0) / perComp.length;
  const meanAbsBias = perComp.reduce((a, b) => a + Math.abs(b.bias), 0) / perComp.length;
  return {
    tau2, studiesPerEdge, nRep: nUsed,
    overallCoverage: overallCov,
    meanAbsBias,
    meanTau2Hat: tau2Sum.v / nUsed,
    perComp,
  };
}

// Deterministic O'Brien-Fleming + Bonferroni boundary check.
// For C comparisons and overall alpha, adjustedAlpha = alpha/C,
// boundary at info fraction t is z_{adjAlpha/2}/sqrt(t), and at t=1
// it must equal the critical value z_{adjAlpha/2}.
export function measureBoundary(alpha, T) {
  const C = T * (T - 1) / 2;
  const adjAlpha = alpha / C;
  const zCrit = qnorm(1 - adjAlpha / 2);
  const b_at_1 = obrienFlemingBoundary(1.0, adjAlpha);
  const b_at_quarter = obrienFlemingBoundary(0.25, adjAlpha);
  return {
    alpha, T, C, adjAlpha,
    zCrit,
    boundaryAt1: b_at_1,                 // should equal zCrit
    boundaryAtQuarter: b_at_quarter,     // should equal zCrit/sqrt(.25)=2*zCrit
    obfRatio: b_at_quarter / b_at_1,     // should be 1/sqrt(.25) = 2
  };
}

function fmt(x, d = 4) { return Number(x).toFixed(d); }

function main() {
  console.log('=== SeqNMA truth-recovery harness ===');
  console.log('True basic params (vs A):', JSON.stringify(TRUE_D));
  console.log('All true contrasts:', allTrueContrasts().map(c => `${c.comp}=${c.theta}`).join('  '));
  console.log('');

  const tau2Levels = [0, 0.05, 0.15];
  const speLevels = [3, 8];
  const NREP = 2000;

  console.log('--- NMA relative-effect 95% CI coverage / bias / tau2-hat ---');
  console.log('tau2   studies/edge  N_studies  overallCov  meanAbsBias  meanTau2Hat');
  const cells = [];
  for (const tau2 of tau2Levels) {
    for (const spe of speLevels) {
      const r = measureCoverage(tau2, spe, NREP, 12345 + Math.round(tau2 * 1000) + spe);
      cells.push(r);
      const nStud = spe * 6;
      console.log(`${fmt(tau2,2)}   ${spe}             ${nStud}         ${fmt(r.overallCoverage,3)}       ${fmt(r.meanAbsBias,4)}       ${fmt(r.meanTau2Hat,4)}`);
    }
  }

  console.log('');
  console.log('--- O\'Brien-Fleming + Bonferroni boundary (deterministic) ---');
  const b = measureBoundary(0.05, 4);
  console.log(`alpha=0.05  T=4  C=${b.C}  adjAlpha=${fmt(b.adjAlpha,5)}`);
  console.log(`z_crit (qnorm(1-adjAlpha/2)) = ${fmt(b.zCrit)}`);
  console.log(`boundary at t=1   = ${fmt(b.boundaryAt1)}  (expect = z_crit)`);
  console.log(`boundary at t=.25 = ${fmt(b.boundaryAtQuarter)}  (expect = 2*z_crit)`);
  console.log(`OBF ratio b(.25)/b(1) = ${fmt(b.obfRatio)}  (expect 2.0)`);

  return { cells, boundary: b };
}

const _invoked = process.argv[1] && import.meta.url.replace(/\\/g, '/').endsWith(
  process.argv[1].replace(/\\/g, '/').split('/').pop()
);
if (_invoked) {
  main();
}

export { main };
