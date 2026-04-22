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

`scpcR` needs:

- a fitted model
- the data used to fit that model
- spatial coordinates supplied either as `lon` / `lat` or as Euclidean
  coordinates

The example below uses `spur_example` from `spuR`.

```r
library(scpcR)
library(spuR)

data(spur_example)

df <- subset(
  spur_example,
  !state %in% c("AK", "HI"),
  select = c(am, gini, fracblack, lat, lon)
)

df <- stats::na.omit(df)

fit <- stats::lm(am ~ gini + fracblack, data = df)

out <- scpc(
  fit,
  data = df,
  lon = "lon",
  lat = "lat"
)

summary(out)
coef(out)
confint(out)
```

- `fit` is the fitted model
- `data` is the underlying data frame
- `lon` and `lat` identify the coordinate columns

If your coordinates are Euclidean rather than geographic, use
`coords_euclidean = c(...)` instead of `lon` and `lat`.

`summary(out)` prints the main SCPC inference table. `coef(out)` extracts the
coefficient estimates, and `confint(out)` extracts confidence intervals.

If you need the diagnostic and transformation stage before inference, the
easiest entry point is `spuR`, which uses `scpcR` internally for the inference
step.

## spuR integration

`spuR` uses `scpcR` internally, so in particular, you can apply `scpc`
as part of the pipeline using:

```r
result <- spuR::spur(
  am ~ gini + fracblack,
  data = df,
  lon = "lon",
  lat = "lat"
)
```

The `scpc` stats in the result object can be accessed using:

```
result$fits$levels$scpc$scpcstats
result$fits$transformed$scpc$scpcstats
```

## Next Step

See [Reference](reference.md) for the full public API, parameter meanings, and
return objects.
