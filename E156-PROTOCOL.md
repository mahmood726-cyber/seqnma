# E156-PROTOCOL: SeqNMA

## Project
Sequential Network Meta-Analysis (SeqNMA)

## Dates
- Created: 2026-04-09
- Status: v1.0

## E156 Body (CURRENT)
Can network meta-analysis adapt sequential monitoring boundaries to determine when cumulative evidence across a treatment network is sufficient? We built a browser-based sequential NMA engine processing contrast-level data (study, year, treatment pair, effect size, standard error) for networks of any topology. The tool implements weighted least-squares NMA with DerSimonian-Laird heterogeneity, O'Brien-Fleming alpha-spending boundaries (z_k = z_alpha/sqrt(t_k)), Bonferroni correction across C = T(T-1)/2 comparisons, and heterogeneity-adjusted required information size using the design effect D^2 = 1 + tau^2*(sum(1/v_i^2)/(sum(1/v_i))^2*k - 1). Applied to a 12-study, 4-treatment demo network, the tool tracks information fractions for all 6 pairwise comparisons, rendering sequential boundary plots, cumulative forest plots, evolving network graphs, and league tables. All 22 Selenium tests pass, confirming correct boundary formulas, monotonic information accrual, Bonferroni adjustment (alpha/6), indirect comparison estimation in star networks, and 2-treatment degeneration to standard TSA. The tool provides real-time guidance on which NMA comparisons have reached conclusive, inconclusive, or futile evidence status, but assumes a consistency model without formal inconsistency testing.

## Dashboard
https://mahmood726-cyber.github.io/SeqNMA/
