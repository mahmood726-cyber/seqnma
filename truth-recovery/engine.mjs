// engine.mjs — pure numerical functions extracted VERBATIM from index.html (SeqNMA).
// No DOM. Line ranges from source index.html as of clone.

// --- phi, qnorm (484-504) ---
function phi(x) {
  const a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741;
  const a4 = -1.453152027, a5 = 1.061405429, p = 0.3275911;
  const sign = x < 0 ? -1 : 1;
  x = Math.abs(x) / Math.SQRT2;
  const t = 1.0 / (1.0 + p * x);
  const y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);
  return 0.5 * (1.0 + sign * y);
}

// Inverse standard normal (rational approximation)
function qnorm(p) {
  if (p <= 0) return -Infinity;
  if (p >= 1) return Infinity;
  if (p === 0.5) return 0;
  if (p > 0.5) return -qnorm(1 - p);
  const t = Math.sqrt(-2 * Math.log(p));
  const c0 = 2.515517, c1 = 0.802853, c2 = 0.010328;
  const d1 = 1.432788, d2 = 0.189269, d3 = 0.001308;
  return -(t - (c0 + c1 * t + c2 * t * t) / (1 + d1 * t + d2 * t * t + d3 * t * t * t));
}

// --- matrix utilities (508-606) ---
function matCreate(rows, cols, fill) {
  const m = [];
  for (let i = 0; i < rows; i++) {
    m[i] = new Array(cols).fill(fill !== undefined ? fill : 0);
  }
  return m;
}

function matTranspose(A) {
  const rows = A.length, cols = A[0].length;
  const T = matCreate(cols, rows);
  for (let i = 0; i < rows; i++)
    for (let j = 0; j < cols; j++)
      T[j][i] = A[i][j];
  return T;
}

function matMul(A, B) {
  const rA = A.length, cA = A[0].length, cB = B[0].length;
  const C = matCreate(rA, cB);
  for (let i = 0; i < rA; i++)
    for (let j = 0; j < cB; j++)
      for (let k = 0; k < cA; k++)
        C[i][j] += A[i][k] * B[k][j];
  return C;
}

function matVecMul(A, v) {
  const r = A.length, c = A[0].length;
  const res = new Array(r).fill(0);
  for (let i = 0; i < r; i++)
    for (let j = 0; j < c; j++)
      res[i] += A[i][j] * v[j];
  return res;
}

function matDiag(v) {
  const n = v.length;
  const D = matCreate(n, n);
  for (let i = 0; i < n; i++) D[i][i] = v[i];
  return D;
}

function matTrace(A) {
  let s = 0;
  for (let i = 0; i < Math.min(A.length, A[0].length); i++) s += A[i][i];
  return s;
}

// Solve Ax = b via Gaussian elimination with partial pivoting
function matSolve(A, b) {
  const n = A.length;
  const aug = A.map((row, i) => [...row, b[i]]);
  for (let col = 0; col < n; col++) {
    let maxRow = col;
    for (let row = col + 1; row < n; row++) {
      if (Math.abs(aug[row][col]) > Math.abs(aug[maxRow][col])) maxRow = row;
    }
    [aug[col], aug[maxRow]] = [aug[maxRow], aug[col]];
    if (Math.abs(aug[col][col]) < 1e-14) return null;
    for (let row = col + 1; row < n; row++) {
      const factor = aug[row][col] / aug[col][col];
      for (let j = col; j <= n; j++) aug[row][j] -= factor * aug[col][j];
    }
  }
  const x = new Array(n);
  for (let i = n - 1; i >= 0; i--) {
    x[i] = aug[i][n];
    for (let j = i + 1; j < n; j++) x[i] -= aug[i][j] * x[j];
    x[i] /= aug[i][i];
  }
  return x;
}

// Invert matrix via augmented row reduction
function matInverse(A) {
  const n = A.length;
  const aug = A.map((row, i) => {
    const r = [...row];
    for (let j = 0; j < n; j++) r.push(i === j ? 1 : 0);
    return r;
  });
  for (let col = 0; col < n; col++) {
    let maxRow = col;
    for (let row = col + 1; row < n; row++) {
      if (Math.abs(aug[row][col]) > Math.abs(aug[maxRow][col])) maxRow = row;
    }
    [aug[col], aug[maxRow]] = [aug[maxRow], aug[col]];
    if (Math.abs(aug[col][col]) < 1e-14) return null;
    const pivot = aug[col][col];
    for (let j = 0; j < 2 * n; j++) aug[col][j] /= pivot;
    for (let row = 0; row < n; row++) {
      if (row === col) continue;
      const factor = aug[row][col];
      for (let j = 0; j < 2 * n; j++) aug[row][j] -= factor * aug[col][j];
    }
  }
  return aug.map(row => row.slice(n));
}

// --- network helpers (638-650) ---
function getTreatments(studies) {
  const set = new Set();
  for (const s of studies) { set.add(s.t1); set.add(s.t2); }
  return [...set].sort();
}

function getAllComparisons(treatments) {
  const comps = [];
  for (let i = 0; i < treatments.length; i++)
    for (let j = i + 1; j < treatments.length; j++)
      comps.push(treatments[i] + ' vs ' + treatments[j]);
  return comps;
}

// --- runNMA (663-774) ---
function runNMA(studies, treatments) {
  if (!treatments) treatments = getTreatments(studies);
  const T = treatments.length;
  const n = studies.length;
  if (T < 2 || n < 1) return null;

  const ref = treatments[0]; // reference
  const params = treatments.slice(1); // T-1 basic parameters

  // Build design matrix X (n x T-1)
  const X = matCreate(n, T - 1);
  const y = [];
  const v = [];

  for (let i = 0; i < n; i++) {
    const s = studies[i];
    y.push(s.effect);
    v.push(s.variance);

    // t2 vs t1: positive effect means t2 better (by convention)
    // If t1 = ref: effect = d_{t2,ref}, column for t2 gets +1
    // If t2 = ref: effect = d_{ref,t1} = -d_{t1,ref}, column for t1 gets -1
    // If neither is ref: effect = d_{t2,ref} - d_{t1,ref}

    const idx1 = params.indexOf(s.t1);
    const idx2 = params.indexOf(s.t2);

    if (idx1 >= 0) X[i][idx1] = -1;
    if (idx2 >= 0) X[i][idx2] = 1;
    // If t1 is reference (idx1 = -1): only +1 for t2
    // If t2 is reference (idx2 = -1): only -1 for t1
  }

  // Step 1: Fixed-effect estimate (tau2 = 0)
  const w0 = v.map(vi => 1 / vi);
  const W0 = matDiag(w0);
  const Xt = matTranspose(X);
  const XtW0 = matMul(Xt, W0);
  const XtW0X = matMul(XtW0, X);
  const XtW0y = matVecMul(XtW0, y);
  const XtW0X_inv = matInverse(XtW0X);
  if (!XtW0X_inv) return null;
  const beta0 = matVecMul(XtW0X_inv, XtW0y);

  // Cochran Q for network
  const Xbeta0 = matVecMul(X, beta0);
  let Q = 0;
  for (let i = 0; i < n; i++) Q += w0[i] * (y[i] - Xbeta0[i]) * (y[i] - Xbeta0[i]);

  const df = n - (T - 1);

  // DerSimonian-Laird tau2 for network
  let tau2 = 0;
  if (df > 0 && Q > df) {
    // Denominator: sum(w) - trace(W*X*(X'WX)^{-1}*X'W * diag(w))
    // Simplified: sum(w) - trace(H * diag(w)) where H = X*(X'WX)^-1*X'W
    const H = matMul(matMul(X, XtW0X_inv), XtW0); // n x n hat matrix
    let trHW = 0;
    for (let i = 0; i < n; i++) trHW += H[i][i] * w0[i];
    const sumW = w0.reduce((a, b) => a + b, 0);
    const denom = sumW - trHW;
    if (denom > 0) tau2 = (Q - df) / denom;
  }

  // Step 2: Random-effects estimate with tau2
  const w = v.map(vi => 1 / (vi + tau2));
  const W = matDiag(w);
  const XtW = matMul(Xt, W);
  const XtWX = matMul(XtW, X);
  const XtWy = matVecMul(XtW, y);
  const XtWX_inv = matInverse(XtWX);
  if (!XtWX_inv) return null;
  const beta = matVecMul(XtWX_inv, XtWy);

  // Variance-covariance of basic parameters
  const varBeta = XtWX_inv;

  // All pairwise comparisons from basic parameters
  const comparisons = [];
  for (let i = 0; i < T; i++) {
    for (let j = i + 1; j < T; j++) {
      const compLabel = treatments[i] + ' vs ' + treatments[j];
      let eff, seEff;

      if (i === 0) {
        // treatment[j] vs reference: beta[j-1]
        eff = beta[j - 1];
        seEff = Math.sqrt(varBeta[j - 1][j - 1]);
      } else {
        // treatment[j] vs treatment[i]: beta[j-1] - beta[i-1]
        eff = beta[j - 1] - beta[i - 1];
        const varDiff = varBeta[j - 1][j - 1] + varBeta[i - 1][i - 1] - 2 * varBeta[j - 1][i - 1];
        seEff = Math.sqrt(Math.max(0, varDiff));
      }

      const z = seEff > 0 ? eff / seEff : 0;
      comparisons.push({
        comp: compLabel,
        t1: treatments[i],
        t2: treatments[j],
        effect: eff,
        se: seEff,
        ci_lo: eff - 1.96 * seEff,
        ci_hi: eff + 1.96 * seEff,
        z: z,
        pval: 2 * (1 - phi(Math.abs(z)))
      });
    }
  }

  return { beta, varBeta, tau2, Q, df, treatments, comparisons, n: studies.length };
}

// --- computeRIS (783-811) ---
function computeRIS(alpha, power, delta, tau2, studies, adjustedAlpha) {
  const za = qnorm(1 - adjustedAlpha / 2);
  const zb = qnorm(power);

  // Estimate sigma2 from study variances (median within-study variance)
  let sigma2;
  if (studies.length > 0) {
    const sortedV = studies.map(s => s.variance).sort((a, b) => a - b);
    sigma2 = sortedV[Math.floor(sortedV.length / 2)];
  } else {
    sigma2 = 0.1;
  }

  const RIS_base = Math.pow(za + zb, 2) * 2 * (sigma2 + tau2) / (delta * delta);

  // Design effect for heterogeneity:
  // D^2 = 1 + tau^2 * (sum(1/v_i^2) / (sum(1/v_i))^2 * k - 1)
  // This is NOT the cluster-design effect
  let D2 = 1;
  if (studies.length > 1 && tau2 > 0) {
    const k = studies.length;
    const sumInvV = studies.reduce((s, st) => s + 1 / st.variance, 0);
    const sumInvV2 = studies.reduce((s, st) => s + 1 / (st.variance * st.variance), 0);
    D2 = 1 + tau2 * ((sumInvV2 / (sumInvV * sumInvV)) * k - 1);
    D2 = Math.max(1, D2);
  }

  return { ris: RIS_base * D2, D2, sigma2 };
}

// --- obrienFlemingBoundary (818-822) ---
function obrienFlemingBoundary(t, adjustedAlpha) {
  if (t <= 0) return Infinity;
  const za = qnorm(1 - adjustedAlpha / 2);
  return za / Math.sqrt(t);
}

export { phi, qnorm, matCreate, matTranspose, matMul, matVecMul, matDiag, matTrace, matSolve, matInverse, getTreatments, getAllComparisons, runNMA, computeRIS, obrienFlemingBoundary };
