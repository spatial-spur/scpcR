<p align="center">
  <img src="assets/logo.png" alt="SPUR logo">
</p>

# scpcR

`scpcR` provides spatial correlation-robust inference for regression coefficients following Müller & Watson 2022 and Müller & Watson 2023, implemented in R based on their [original Stata implementation](https://github.com/ukmueller/SCPC).

**Citation:** If you use this package, please cite Becker, Boll and Voth 2026, Müller & Watson 2022, and Müller & Watson 2023. See [CITATION.bib](CITATION.bib) for the BibTeX entries.

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

## Multithreading and Performance

`scpcR` uses `RcppParallel` for multi-core acceleration in computationally intensive steps (distance matrix calculations, subsampling, and Omega grid evaluations). By default, it automatically uses available cores.

You can configure or restrict multithreading using:

```r
# Restrict to a specific thread count via helper function:
scpc_set_num_threads(2)

# Or set via standard R package option:
options(scpcR.threads = 1)

# Or set via environment variable:
Sys.setenv(SCPCR_NUM_THREADS = "1")

# Check currently active thread count:
scpc_get_num_threads()

# Restore automatic multi-core threading:
scpc_set_num_threads("auto")
```

> **Tip for cluster & batch workflows:** If you are running `scpc()` inside parallel loops (e.g. `parallel::mclapply`, `future`, or a shared cluster node), set `scpc_set_num_threads(1)` or `options(scpcR.threads = 1)` to prevent CPU oversubscription. Setting thread count to `1` completely bypasses background task submission, ensuring pure single-threaded execution without overhead or thread creation issues.

## Documentation

Please refer to [the package documentation](https://spatial-spur.github.io/scpcR/) for detailed information and other (R, Python, Stata) packages.

## References

Becker, Sascha O., P. David Boll and Hans-Joachim Voth "Testing and Correcting for Spatial Unit Roots in Regression Analysis", The Stata Journal 26(2): 177–202. [https://doi.org/10.1177/1536867X261449932](https://journals.sagepub.com/doi/10.1177/1536867X261449932).

Müller, Ulrich K. and Mark W. Watson "Spatial Correlation Robust Inference", Econometrica 90(6) (2022), 2901–2935. https://www.princeton.edu/~umueller/SHAR.pdf.

Müller, Ulrich K. and Mark W. Watson "Spatial Correlation Robust Inference in Linear Regression and Panel Models", Journal of Business & Economic Statistics 41(4) (2023), 1050–1064. https://www.princeton.edu/~umueller/SpatialRegression.pdf.
