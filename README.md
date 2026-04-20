# scpcR

`scpcR` provides spatial correlation-robust inference for regression coefficients
following Mueller and Watson (2022, 2023).

## Installation

Install the development version from GitHub with:

```r
remotes::install_github("pdavidboll/scpcR")
```

## Example

```r
set.seed(42)
n <- 60
dat <- data.frame(
  y = rnorm(n),
  x = rnorm(n),
  lon = runif(n, -100, -80),
  lat = runif(n, 30, 45)
)

fit <- lm(y ~ x, data = dat)

out <- scpcR::scpc(
  fit,
  data = dat,
  coords_euclidean = c("lon", "lat"),
  avc = 0.1,
  uncond = TRUE
)

summary(out)
```

For the matching Python implementation, see
[`scpc-python`](https://github.com/DGoettlich/scpc-python).
