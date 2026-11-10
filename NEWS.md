# Development changes

- Compute rejection probabilities using an adaptive angular integral; the zero statistic has p-value one, and the independent-error case agrees with Student-t probabilities.
- Solve bandwidth equations to relative accuracy without a fixed upper bound in coordinate units; reject unattainable correlation targets.
- Use great-circle distances throughout the geographic Nyström approximation and partial Fisher-Yates sampling without replacement.
- Retain the full `W'W` diagonal contribution for conditional, nonorthonormal weights.
- Bound kernel scratch memory by projecting 256 rows at a time directly into the stored basis.

Critical-value calibration uses a finite covariance grid. Its numerical tolerance does not bound continuous-grid or Monte Carlo approximation error. Exact upstream parity is not an acceptance criterion.

- Apply the CRAN two-thread cap after every thread-setting source; return the original option from thread setters so automatic detection can be restored.

Upstream comparison suites are opt-in diagnostics (`SCPC_RUN_UPSTREAM_PARITY=true`). Default regression checks use analytic identities, defining matrix formulas, invariance, and Monte Carlo uncertainty.
