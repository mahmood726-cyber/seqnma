// ============================================================
// dgp-nma.mjs — Known-truth DGP for a contrast-based network
// meta-analysis. Seeded (mulberry32) -> fully reproducible.
//
// Truth model:
//   Treatments A,B,C,D with A as reference.
//   Basic parameters d (relative to A):
//     d_A = 0 (reference), d_B, d_C, d_D KNOWN.
//   Every true relative effect is therefore known:
//     theta_true(x,y) = d_y - d_x.
//   Common between-study heterogeneity variance tau2 (shared across
//   the whole network — the consistency + common-tau model that
//   SeqNMA's runNMA assumes).
//
//   For a study comparing arms (t1,t2) with within-study SE s:
//     study random effect:  delta ~ N(theta_true(t1,t2), tau2)
//     observed contrast:     y ~ N(delta, s^2)
//   so y ~ N(theta_true, tau2 + s^2). This is exactly the
//   random-effects contrast model. A correct RE-NMA estimator must
//   cover theta_true(t1,t2) at its nominal 95% rate.
//
// The network is CONNECTED by construction: a spanning set of
// edges (A-B, A-C, A-D) is always present, plus optional extra
// edges to make it a richer (non-star) network with indirect
// information.
// ============================================================

export function makeRng(seed) {
  let a = seed >>> 0;
  return function () {
    a |= 0; a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

function randn(rng) {
  let u1 = rng(), u2 = rng();
  if (u1 < 1e-12) u1 = 1e-12;
  return Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
}

function drawSe(rng, seLo, seHi) {
  const lo = Math.log(seLo), hi = Math.log(seHi);
  return Math.exp(lo + (hi - lo) * rng());
}

// Known truth. d is keyed by treatment letter; A is reference (0).
export const TRUE_D = { A: 0, B: 0.40, C: -0.25, D: 0.65 };

export function trueTheta(t1, t2) {
  // contrast convention matches runNMA: effect of t2 relative to t1 = d_{t2} - d_{t1}
  return TRUE_D[t2] - TRUE_D[t1];
}

// Edge design: spanning edges keep the network connected; the extra
// edges add direct evidence on additional contrasts (so the network
// is a genuine loop-containing network, not a pure star).
const SPANNING = [['A', 'B'], ['A', 'C'], ['A', 'D']];
const EXTRA    = [['B', 'C'], ['C', 'D'], ['B', 'D']];

/**
 * Generate one network meta-analysis dataset of known truth.
 * @param {number} tau2          true common between-study variance
 * @param {number} studiesPerEdge studies on each edge (>=1)
 * @param {function} rng         seeded PRNG
 * @param {object} opts          { seLo, seHi, year0 }
 * @returns {{studies: Array, treatments: string[]}}
 *   studies: [{study, year, t1, t2, effect, se, variance}]
 *   (t1<t2 alphabetical, matching parseCSV's normalisation)
 */
export function generateNetwork(tau2, studiesPerEdge, rng,
                                { seLo = 0.12, seHi = 0.55, year0 = 2000 } = {}) {
  const sd = Math.sqrt(tau2);
  const edges = [...SPANNING, ...EXTRA];
  const studies = [];
  let idx = 0;
  for (const [a, b] of edges) {
    for (let r = 0; r < studiesPerEdge; r++) {
      const s = drawSe(rng, seLo, seHi);
      // ensure alphabetical t1<t2 (a,b already are by construction)
      const t1 = a, t2 = b;
      const theta = trueTheta(t1, t2);            // true relative effect t2 vs t1
      const delta = theta + sd * randn(rng);       // study-specific true effect
      const y = delta + s * randn(rng);            // observed contrast
      studies.push({
        study: `S${String(idx + 1).padStart(2, '0')}`,
        year: year0 + idx,
        t1, t2,
        effect: y,
        se: s,
        variance: s * s,
      });
      idx++;
    }
  }
  studies.sort((x, y2) => x.year - y2.year || x.study.localeCompare(y2.study));
  return { studies, treatments: ['A', 'B', 'C', 'D'] };
}

// The list of all 6 pairwise contrasts and their known true thetas.
export function allTrueContrasts() {
  const T = ['A', 'B', 'C', 'D'];
  const out = [];
  for (let i = 0; i < T.length; i++)
    for (let j = i + 1; j < T.length; j++)
      out.push({ comp: `${T[i]} vs ${T[j]}`, t1: T[i], t2: T[j], theta: trueTheta(T[i], T[j]) });
  return out;
}
