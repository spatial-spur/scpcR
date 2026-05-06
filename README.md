<p align="center">
  <img src="assets/logo.png" alt="SPUR logo">
</p>

# scpcR

`scpcR` provides spatial correlation-robust inference for regression coefficients following Müller and Watson (2022, 2023), implemented in R based on their [original Stata implementation](https://github.com/ukmueller/SCPC).

**Citation:** If you use this package, please cite [Becker et al. 2026](CITATION.cff), [Müller & Watson 2022](CITATION.cff), and [Müller & Watson 2023](CITATION.cff). See [CITATION.cff](CITATION.cff) for copyable citation metadata.

If you encounter any issues or have any questions, please open an issue on GitHub or contact the authors.

## Installation

Install the latest release version from GitHub with:

```r
remotes::install_github("spatial-spur/scpcR@v0.1.3")
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
  lon = "lon",
  lat = "lat"
)

summary(out)
```

## Documentation

Please refer to [the package documentation](https://spatial-spur.github.io/scpcR/) for detailed information and other (R, Python, Stata) packages.

## References

Becker, Sascha O., P. David Boll and Hans-Joachim Voth "Testing and Correcting for Spatial Unit Roots in Regression Analysis", The Stata Journal, in press.

Müller, Ulrich K. and Mark W. Watson "Spatial Correlation Robust Inference", Econometrica 90(6) (2022), 2901–2935. https://www.princeton.edu/~umueller/SHAR.pdf.

Müller, Ulrich K. and Mark W. Watson "Spatial Correlation Robust Inference in Linear Regression and Panel Models", Journal of Business & Economic Statistics 41(4) (2023), 1050–1064. https://www.princeton.edu/~umueller/SpatialRegression.pdf.
