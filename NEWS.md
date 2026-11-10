# Development changes

- Reduce bandwidth and critical-value root evaluations with safeguarded Illinois interpolation; retain relative bandwidth tolerance `1e-12` and critical-value tolerance `1e-9`.
- Cache angular quadrature nodes and evaluate bounded products without repeated logarithms. Rescale the complement near probability one to resolve narrow endpoint layers.
- Precompute unit-sphere positions for geographic cross-kernels, with accurate nearby-point and antipodal fallbacks. Retain the 256-row memory bound.

- Compute rejection probabilities using an adaptive angular integral; the zero statistic has p-value one, and the independent-error case agrees with Student-t probabilities.
- Solve bandwidth equations to relative accuracy without a fixed upper bound in coordinate units; reject unattainable correlation targets.
- Use great-circle distances throughout the geographic Nyström approximation and partial Fisher-Yates sampling without replacement.
- Retain the full `W'W` diagonal contribution for conditional, nonorthonormal weights.
- Bound kernel scratch memory by projecting 256 rows at a time directly into the stored basis.

Critical-value calibration uses a finite covariance grid. Its numerical tolerance does not bound continuous-grid or Monte Carlo approximation error. Exact upstream parity is not an acceptance criterion.

Upstream comparison suites are opt-in diagnostics (`SCPC_RUN_UPSTREAM_PARITY=true`). Default regression checks use analytic identities, defining matrix formulas, invariance, and Monte Carlo uncertainty.
