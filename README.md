# SeqNMA -- Sequential Network Meta-Analysis

Browser-based tool extending Trial Sequential Analysis (TSA) to network meta-analysis. Provides O'Brien-Fleming monitoring boundaries for all pairwise comparisons in a treatment network as new studies accumulate.

## Features

- **Contrast-based NMA engine**: weighted least squares with DerSimonian-Laird heterogeneity estimation
- **Sequential monitoring**: O'Brien-Fleming alpha-spending boundaries with Bonferroni correction for multiple comparisons
- **Required Information Size**: adjusted by design effect for heterogeneity (D^2 formula, not cluster DEFF)
- **Futility boundaries**: non-binding, conditional-power-based
- **5 interactive tabs**:
  1. Sequential Boundary Plot (SVG) with efficacy + futility boundaries
  2. Cumulative NMA Forest showing effect evolution over time
  3. Network Graph Evolution with year slider (force-directed layout)
  4. Information Tracker table + bar chart with status badges
  5. League Table (T x T matrix of all pairwise comparisons)
- **Export**: CSV results and SVG plots
- **Offline**: single-file HTML, no CDN dependencies

## Usage

1. Open `index.html` in a browser
2. Click "Demo Data" or paste CSV with columns: Study, Year, Treat1, Treat2, Effect, SE
3. Set effect scale (SMD/logOR/logRR/logHR), delta, alpha, power
4. Click "Analyze"

## Statistical Methods

- NMA via weighted least squares with design matrix X mapping contrasts to basic parameters
- Heterogeneity: network DerSimonian-Laird (generalized Q-based)
- Boundaries: z_k = z_{alpha/2} / sqrt(t_k) (O'Brien-Fleming alpha-spending)
- Multiple testing: Bonferroni alpha/C where C = T(T-1)/2
- Design effect: D^2 = 1 + tau^2 * (sum(1/v_i^2) / (sum(1/v_i))^2 * k - 1)

## Testing

```bash
cd C:\Models\SeqNMA
python -m pytest test_app.py -v
```

22 tests covering: NMA engine, boundary formulas, Bonferroni, monotonic information, 2-/3-/4-treatment networks, star topology indirect comparisons, league table, export, UI rendering.

## Author

Mahmood Ahmad, Tahir Heart Institute
