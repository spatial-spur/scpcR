# Reference

This page documents the public `scpcR` API.

## Overview

| Function | Description |
|---|---|
| `scpc()` | Run spatial correlation-robust inference for a fitted model |
| `print()` | Print the main SCPC coefficient table |
| `summary()` | Print the coefficient table plus confidence intervals |
| `coef()` | Extract the coefficient estimates |
| `confint()` | Extract confidence intervals for selected coefficients |

## Conventions

### Supported Models

`scpc()` accepts fitted:

- `lm` models
- `fixest::feols()` models, including IV specifications

### Coordinates

Provide exactly one coordinate specification:

- `lon` and `lat` for geographic coordinates
- `coords_euclidean` for planar coordinates

Do not pass both specifications in the same call.

### Clustering and Method Selection

- `cluster` names an optional clustering variable in `data`
- coordinates must be constant within clusters
- when a model uses absorbed fixed effects and external clustering, the
  conditional adjustment is not implemented; use `uncond = TRUE`
- `method = "auto"` selects the exact branch for smaller problems and the
  approximation branch once the large-`n` threshold is reached

## Core Inference

### `scpc()`

Run spatial correlation-robust inference.

**Signature**

```r
scpc(
  model,
  data,
  lon = NULL,
  lat = NULL,
  coords_euclidean = NULL,
  cluster = NULL,
  ncoef = NULL,
  avc = 0.03,
  method = "auto",
  large_n_seed = 1,
  uncond = FALSE,
  cvs = FALSE
)
```

**Parameters**

- `model`: fitted model object. Supported classes are `lm` and
  `fixest::feols()`, including IV models.
- `data` (`data.frame`): data used to fit `model`. It must contain the
  coordinate columns referenced by `lon` / `lat` or `coords_euclidean`.
- `lon`, `lat` (`character | NULL`): longitude and latitude column names for
  geographic distance calculations.
- `coords_euclidean` (`character | NULL`): Euclidean coordinate column names.
  Use instead of `lon` and `lat`.
- `cluster` (`character | NULL`): optional clustering variable in `data`.
- `ncoef` (`integer | NULL`): number of coefficients to report. `NULL` reports
  all coefficients.
- `avc` (`numeric`): upper bound on average pairwise correlation. Must lie in
  `(0.001, 0.99)`. Default: `0.03`.
- `method` (`character`): spatial algorithm.
  - `"auto"`: choose between exact and approximate branches
  - `"exact"`: always use the full distance matrix
  - `"approx"`: always use the large-`n` approximation branch
- `large_n_seed` (`numeric`): integer-valued seed used by the large-`n`
  approximation branch. Ignored when `method = "exact"`.
- `uncond` (`logical`): if `TRUE`, skip the conditional adjustment and report
  unconditional critical values only.
- `cvs` (`logical`): if `TRUE`, store per-coefficient critical values at the
  32%, 10%, 5%, and 1% levels.

**Returns**

- object of class `"scpc"`
  - `scpcstats`: matrix with coefficient estimates, SCPC standard errors,
    t-statistics, p-values, and 95% interval bounds
  - `scpccvs`: matrix of stored critical values, or `NULL` when `cvs = FALSE`
  - `W`: spatial projection matrix
  - `avc`
  - `c0`
  - `cv`
  - `q`
  - `method`
  - `large_n_seed`
  - `call`

## Methods

### `print.scpc()`

Print the main SCPC coefficient table.

**Signature**

```r
print(x, ...)
```

**Parameters**

- `x`: object of class `"scpc"`

**Returns**

- the input object, invisibly

The printed output shows the coefficient table with estimates, SCPC standard
errors, t-statistics, and p-values. If `cvs = TRUE`, stored critical values
are printed as well.

### `summary.scpc()`

Print the fuller SCPC summary.

**Signature**

```r
summary(object, ...)
```

**Parameters**

- `object`: object of class `"scpc"`

**Returns**

- the input object, invisibly

The summary adds the 95% confidence interval table and reports `avc` alongside
the main coefficient output.

### `coef.scpc()`

Extract the coefficient estimates.

**Signature**

```r
coef(object, ...)
```

**Parameters**

- `object`: object of class `"scpc"`

**Returns**

- named numeric vector of coefficient estimates

### `confint.scpc()`

Extract confidence intervals for selected coefficients.

**Signature**

```r
confint(object, parm = NULL, level = 0.95, ...)
```

**Parameters**

- `object`: object of class `"scpc"`
- `parm` (`character | numeric | NULL`): coefficient names or indices to
  include. `NULL` uses all coefficients.
- `level` (`numeric`): confidence level to request.
  - `0.95` is always available
  - `0.68`, `0.90`, and `0.99` are also available when `cvs = TRUE`

**Returns**

- numeric matrix of lower and upper interval bounds

## Return Object

### `"scpc"`

Returned by `scpc()`.

**Fields**

- `scpcstats`
- `scpccvs`
- `W`
- `avc`
- `c0`
- `cv`
- `q`
- `method`
- `large_n_seed`
- `call`

**Methods**

- `print()`
- `summary()`
- `coef()`
- `confint()`

The primary inferential output is stored in `scpcstats`. Its rows correspond to
coefficients, and its columns contain the estimate, standard error,
t-statistic, p-value, and 95% interval endpoints.
