<p align="center">
  <img src="assets/logo.png" alt="SPUR logo">
</p>

# scpcR

`scpcR` provides spatial correlation-robust inference for regression coefficients following Müller & Watson 2022 and Müller & Watson 2023, implemented in R based on their [original Stata implementation](https://github.com/ukmueller/SCPC).

**Citation:** If you use this package, please cite Becker, Boll and Voth 2026, Müller & Watson 2022, and Müller & Watson 2023. See [CITATION.bib](CITATION.bib) for the BibTeX entries.

If you encounter any issues or have any questions, please open an issue on GitHub or contact the authors.

## Installation

For this development checkout, install from the workspace root with Cargo and
Rust 1.89 or newer on your PATH. The Rust backend uses OxiBLAS for matrix
products and dense symmetric eigendecomposition.

Set `SPATIAL_RUST=0` for a pure R build. `backend = "fallback"` selects pure R at
runtime in either build.

```sh
R CMD INSTALL scpcR
```

For a standalone source archive that can build Rust without network access, run
`bash spatial-rust/scripts/build_r_source.sh scpcR ../r-source`, then install the
resulting archive with `R CMD INSTALL`.

OxiBLAS uses SME2 for eligible matrix products on supported Arm64 macOS
and Linux systems. Set `OXIBLAS_DISABLE_SME2=1` before starting R to
use the other OxiBLAS kernels.

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
