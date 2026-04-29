<p align="center">
  <img src="assets/logo.png" alt="SPUR logo">
</p>

# scpcR

`scpcR` provides spatial correlation-robust inference for regression coefficients following Müller and Watson (2022, 2023), implemented in R based on their [original Stata implementation](https://github.com/ukmueller/SCPC).

**When using this code, please cite [Becker, Boll and Voth (2026)](https://pauldavidboll.com/SPUR_Stata_Journal_website.pdf):**

```bibtex
@Article{becker2026,
  author    = {Becker, Sascha O. and Boll, P. David and Voth, Hans-Joachim},
  title     = {Testing and Correcting for Spatial Unit Roots in Regression Analysis},
  journal   = {Stata Journal},
  year      = {forthcoming},
  note      = {Forthcoming}
}
```

If you encounter any issues or have any questions, please open an issue on GitHub or contact the authors.

## Installation

Install the development version from GitHub with:

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
  coords_euclidean = c("lon", "lat"),
  avc = 0.1,
  uncond = TRUE
)

summary(out)
```

## Documentation

Please refer to [the package documentation](https://spatial-spur.github.io/scpcR/) for detailed information and other (R, Python, Stata) packages.

## References

Becker, Sascha O., P. David Boll and Hans-Joachim Voth "Testing and Correcting for Spatial Unit Roots in Regression Analysis", Forthcoming at the Stata Journal.

Müller, Ulrich K. and Mark W. Watson "Spatial Correlation Robust Inference", Econometrica 90(6) (2022), 2901–2935. https://www.princeton.edu/~umueller/SHAR.pdf.

Müller, Ulrich K. and Mark W. Watson "Spatial Correlation Robust Inference in Linear Regression and Panel Models", Journal of Business & Economic Statistics 41(4) (2023), 1050–1064. https://www.princeton.edu/~umueller/SpatialRegression.pdf.
