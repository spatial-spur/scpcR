# scpcR

`scpcR` provides the SCPC inference stage of the SPUR workflow for
cross-sectional regressions with spatial dependence. For the diagnostic and
transformation stage in R, see `spuR`.

## Installation

`scpcR` can be installed as standalone package with 

```r
# install.packages("remotes")
remotes::install_github("spatial-spur/scpcR@v0.1.3")
```

If you want to use `scpcR` through the full `spuR` workflow, install both
tagged releases explicitly:
```r
remotes::install_github("spatial-spur/scpcR@v0.1.3")
remotes::install_github("spatial-spur/spuR@v0.1.2")
```

## Example: Chetty Dataset

In this example, we walk you through the SCPC inference workflow step-by-step.
SCPC is the inference stage of the SPUR workflow, so we start from a fitted
regression model. We also provide a one-stop [pipeline wrapper](#pipeline-wrapper)
if you want to run the full SPUR workflow in one step.

### Data preparation

For illustration, we load the Chetty dataset from `spuR`. Of course, the
analysis in principle follows the same logic on any other dataset. In this
specific case, we first omit the non-contiguous US states. We also drop rows
with missing values.

```r
library(spuR)

data(spur_example)

df <- subset(
  spur_example,
  !state %in% c("AK", "HI"),
  select = c(am, gini, fracblack, lat, lon)
)

df <- stats::na.omit(df)
```

### Fitting the regression

SCPC needs a fitted model, the data used to fit that model, and spatial
coordinates supplied either as `lon` / `lat` or as Euclidean coordinates.

```r
library(scpcR)

fit <- stats::lm(am ~ gini + fracblack, data = df)
```

### Running SCPC inference

We suggest applying SCPC inference after estimating the regression:

```r
out <- scpc(
  fit,
  data = df,
  lon = "lon",
  lat = "lat"
)

summary(out)
```

### Interpreting the output

`summary(out)` prints the main SCPC inference table. You can also use
`coef(out)` to extract coefficient estimates and `confint(out)` to extract
confidence intervals.

If your coordinates are Euclidean rather than geographic, use
`coords_euclidean = c(...)` instead of `lon` and `lat`.

### Pipeline wrapper

As a shortcut to implementing the full SPUR workflow manually, use the `spuR`
pipeline wrapper. It runs the diagnostics, transformation, and SCPC inference
steps in one call.

```r
result <- spuR::spur(
  am ~ gini + fracblack,
  data = df,
  lon = "lon",
  lat = "lat"
)
```

The nested SCPC results can be printed directly:

```
summary(result$fits$levels$scpc)
summary(result$fits$transformed$scpc)
```

## Next Step

See [Reference](reference.md) for the full public API, parameter meanings, and
return objects.
