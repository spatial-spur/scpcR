---
name: Bug report
about: Report a reproducible problem with scpcR
title: "[bug] "
labels: bug
---

## Summary

Describe the problem in 1-3 sentences.

## Minimal reproduction

Provide a copy-pasteable example if possible.

```r

```

## Expected behavior

What did you expect to happen?

## Actual behavior

What happened instead? Include the full traceback, warning, or numerical mismatch if relevant.

```text

```

## Environment

- OS:
- R version:
- install method: `install.packages` / `remotes::install_github` / `devtools::load_all` / local source install / other
- scpcR version or commit:
- relevant model type: `lm` / `fixest::feols` / IV / other
- affected function or workflow: `scpc` / `summary.scpc` / `coef.scpc` / `confint.scpc` / coordinates / clustering / docs / CI

## Additional context

Anything else that might help reproduce or explain the issue.

- Does this reproduce on synthetic data, example data, or only your own data?
- Is the issue a crash, wrong result, docs mismatch, install problem, or performance regression?
- If this is a numerical mismatch, does it depend on `method = "exact"` vs `"approx"` or on `uncond = TRUE/FALSE`, and do you see the same issue relative to `scpc-python`?
