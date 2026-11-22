<p align="center">
  <img src="assets/logo.png" alt="SPUR logo">
</p>

# scpcR

`scpcR` provides spatial correlation-robust inference for regression coefficients following Müller & Watson 2022 and Müller & Watson 2023, implemented in R based on their [original Stata implementation](https://github.com/ukmueller/SCPC).

**Citation:** If you use this package, please cite Becker, Boll and Voth 2026, Müller & Watson 2022, and Müller & Watson 2023. See [CITATION.bib](CITATION.bib) for the BibTeX entries.

If you encounter any issues or have any questions, please open an issue on GitHub or contact the authors.

## Installation

For this development checkout, install from the workspace root with Cargo on your PATH.
Set `SPATIAL_RUST=0` for a pure R build; `backend = "fallback"` selects pure R at
runtime in either build.

```sh
R CMD INSTALL scpcR
```

For a standalone source archive that can build Rust without network access, run
`bash spatial-rust/scripts/build_r_source.sh scpcR ../r-source`, then install the
resulting archive with `R CMD INSTALL`.

On Windows install the Rust target matching Rtools and your R installation.
The build uses `x86_64-pc-windows-gnu` for x86-64 R or
`aarch64-pc-windows-gnullvm` for ARM64 R. Set `WIN_TARGET` to override
the selection. Windows native execution still needs validation.

For the published version:

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

Becker, Sascha O., P. David Boll and Hans-Joachim Voth "Testing and Correcting for Spatial Unit Roots in Regression Analysis", The Stata Journal 26(2): 177–202. [https://doi.org/10.1177/1536867X261449932](https://journals.sagepub.com/doi/10.1177/1536867X261449932).

Müller, Ulrich K. and Mark W. Watson "Spatial Correlation Robust Inference", Econometrica 90(6) (2022), 2901–2935. https://www.princeton.edu/~umueller/SHAR.pdf.

Müller, Ulrich K. and Mark W. Watson "Spatial Correlation Robust Inference in Linear Regression and Panel Models", Journal of Business & Economic Statistics 41(4) (2023), 1050–1064. https://www.princeton.edu/~umueller/SpatialRegression.pdf.

## Shared-core conditional projections

Conditional OLS includes the full fitted design, including absorbed fixed effects,
before cluster aggregation. Conditional IV uses the transpose of its oblique residual
maker, so identical instrument and regressor spaces reduce to OLS. These corrections
apply to both native and fallback calculations and can change intervals and p-values
from archived releases. `backend` and `nthreads` are independent options.

Each interval covers a single coefficient under the corresponding SCPC
assumptions. The correlation bound applies within the exponential benchmark,
and the full-matrix method evaluates a finite decay grid.

### IV inference

The SCPC t-statistic for a 2SLS coefficient requires strong identification
and the spatial score-process conditions in Müller & Watson (2023, §5.3).
Conditional Gaussian size control requires the stated error model given
the endogenous regressors, instruments, and exogenous controls. Instrument
exogeneity alone is insufficient for that conditional claim. Where the
unconditional argument applies, taking the larger critical value preserves
its protection.

For one endogenous regressor `p` and one excluded instrument `z`, Anderson–Rubin
inference can be built from an SCPC test of the instrument coefficient in
the regression of `y - b0 * p` on `z` and the exogenous controls. The example
uses control `w` and coordinate columns `lon` and `lat` in `d`.

```r
b0 <- 0
fit <- lm(I(y - b0 * p) ~ z + w, data = d)
scpcR::scpc(fit, d, lon = "lon", lat = "lat")$scpcstats["z", "P>|t|"]
```

Retain the values of `b0` for which the instrument test fails to reject. The
resulting confidence set may be disconnected or unbounded.
