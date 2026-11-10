test_that("SCPC probabilities and critical values satisfy the Student-t identity", {
  for (q in c(1L, 2L, 3L, 10L, 60L)) {
    for (cv in c(0, 1e-5, 0.01, 0.1, 1, 10)) {
      expect_equal(.getrp(diag(q + 1L), cv), 2 * pt(-sqrt(q) * cv, q), tolerance = 1e-9)
    }
  }
  for (q in c(1L, 2L, 10L)) {
    for (level in c(0.32, 0.05, 0.01)) {
      cv <- .getcv(list(diag(q + 1L)), q, level)
      expect_equal(cv, qt(1 - level / 2, q), tolerance = 1e-8)
      expect_lte(.maxrp(list(diag(q + 1L)), q, cv / sqrt(q))$max, level + 1e-10)
    }
  }
})

test_that("bandwidth solving is invariant to coordinate units", {
  distances <- c(0.2, 0.7, 1, 1.3, 2, 4)
  roots <- vapply(c(1e-12, 1, 1e12), function(scale) {
    c0 <- .getc0fromavc(distances * scale, 0.03)
    expect_equal(mean(exp(-c0 * distances * scale)), 0.03, tolerance = 1e-10)
    c0 * scale
  }, 0.0)
  expect_equal(roots, rep(roots[2], 3), tolerance = 1e-9)
  expect_error(.getc0fromavc(c(0, 0, 1), 0.1), "co-located")
})

test_that("sampled covariance retains the diagonal of a nonorthonormal basis", {
  for (n in c(50L, 1300L)) {
    S <- cbind(seq_len(n), 2 * seq_len(n))
    W <- cbind(rep(2, n), seq_len(n) / n, sin(seq_len(n)))
    out <- .lnget_Oms(S, 1e6, 2e6, W, 1.2, capM = 2000L, random_t = 7)
    for (om in out$Oms) expect_equal(om, crossprod(W), tolerance = 1e-12)
  }
})

test_that("a full geographic subsample spans the exact spatial eigenspace", {
  set.seed(581)
  S <- cbind(runif(40, -179, 179), runif(40, -70, 70))
  basis <- .lnset_wc0(S, 0.1, 3L, 1e-4, TRUE, capN = 1L, m = 40L, random_t = 7)
  exact <- .getW(.getdistmat(S, TRUE), basis$c0, 3L)[, -1, drop = FALSE]
  approximate <- basis$W[, -1, drop = FALSE]
  expect_equal(svd(crossprod(exact, approximate))$d, rep(1, 3), tolerance = 1e-8)
})

test_that("partial shuffling samples locations uniformly without replacement", {
  S <- cbind(seq_len(4), seq_len(4))
  seeds <- seq(0, 2^32 - 1, length.out = 4000)
  counts <- integer(4)
  for (seed in seeds) {
    out <- .jumble_s(S, 2L, floor(seed))$coords
    if (seed %in% seeds[c(1, 1000, 4000)]) expect_equal(sort(out[, 1]), seq_len(4))
    counts <- counts + tabulate(out[1:2, 1], nbins = 4)
  }
  expect_lt(max(abs(counts / length(seeds) - 0.5)), 0.025)
})


test_that("correlated Gaussian rejection probabilities agree with simulation uncertainty", {
  set.seed(7401)
  for (q in c(1L, 2L, 5L)) {
    A <- matrix(rnorm((q + 1L)^2), q + 1L)
    Om <- crossprod(A) + diag(q + 1L)
    draws <- t(chol(Om)) %*% matrix(rnorm((q + 1L) * 100000L), q + 1L)
    denominator <- sqrt(colSums(draws[-1L, , drop = FALSE]^2))
    for (cv in c(0.3, 1, 2)) {
      estimated <- mean(abs(draws[1L, ]) > cv * denominator)
      probability <- .getrp(Om, cv)
      mc_error <- sqrt(probability * (1 - probability) / ncol(draws))
      expect_lt(abs(estimated - probability), 6 * mc_error + 1e-4)
    }
  }
})


test_that("angular probabilities match independent adaptive integration", {
  reference <- function(Om, cv) {
    q <- nrow(Om) - 1L
    if (cv == 0) return(1)
    U <- chol(Om)
    roots <- eigen(U %*% diag(c(1, rep(-cv^2, q))) %*% t(U), symmetric = TRUE)$values
    ratios <- -roots[roots < 0] / max(roots)
    integrand <- function(theta) vapply(theta, function(angle) {
      exp(-0.5 * sum(log1p(ratios / cos(angle)^2)))
    }, 0.0)
    # Explicit dyadic panels keep the independent Gauss-Kronrod oracle from
    # overlooking a narrow endpoint layer when the cutoff is close to zero.
    end <- pi/2
    depth <- max(0L, min(40L, ceiling(log2(end / sqrt(sum(ratios))))))
    breaks <- if (depth == 0L) c(0, end) else c(0, end - end / 2^seq_len(depth), end)
    values <- vapply(seq_len(length(breaks) - 1L), function(i) {
      integrate(integrand, breaks[i], breaks[i + 1L],
        abs.tol = 1e-12 / length(breaks), rel.tol = 1e-11)$value
    }, 0.0)
    2 / pi * sum(values)
  }
  set.seed(825)
  for (q in c(2L, 5L, 12L)) {
    A <- matrix(rnorm((q + 1L)^2), q + 1L)
    Om <- crossprod(A) + diag(q + 1L)
    for (cv in c(1e-5, 0.01, 0.3, 1, 5)) {
      expect_lt(abs(.getrp(Om, cv) - reference(Om, cv)), 2e-10)
    }
    if (q <= 5L) {
      Oms <- list(diag(q + 1L), Om)
      for (level in c(0.05, 0.01)) {
        worst <- function(cv) max(vapply(Oms, reference, 0.0, cv = cv))
        upper <- 1
        while (worst(upper) > level) upper <- 2 * upper
        root <- uniroot(function(cv) worst(cv) - level, c(0, upper), tol = 1e-11)$root
        actual <- .getcv(Oms, q, level)
        expect_equal(actual, sqrt(q) * root, tolerance = 1e-8)
        expect_lte(worst(actual / sqrt(q)), level + 2e-10)
      }
    }
  }
})

test_that("quadrature remains accurate for extreme tails and covariance scales", {
  for (q in c(2L, 30L, 120L, 300L)) {
    for (cv in c(1e-8, 0.1, 1, 10)) {
      expect_lt(abs(.getrp(diag(q + 1L), cv) - 2 * pt(-sqrt(q) * cv, q)), 2e-10)
    }
  }
  Om <- matrix(c(4, 1.2, 1.2, 2), 2L)
  location <- Om[1, 2] / Om[2, 2]
  scale <- sqrt(det(Om)) / Om[2, 2]
  for (cv in c(0, 1e-5, 0.1, 1, 10, 1e6)) {
    expected <- pcauchy(-cv, location, scale) + pcauchy(cv, location, scale, lower.tail = FALSE)
    for (units in c(1e-100, 1, 1e100)) {
      expect_lt(abs(.getrp(Om * units, cv) - expected), 1e-10)
    }
  }
})

test_that("safeguarded bandwidth interpolation respects attainable correlation targets", {
  distances <- c(rep(0, 4), 0.01, 0.02, 0.1, 0.5, 1, 2, 10, 100)
  floor <- mean(distances == 0)
  targets <- floor + (1 - floor) * c(0.99, 0.03, 1e-5)
  for (units in c(1e-12, 1, 1e12)) {
    roots <- .getc0fromavc(distances * units, targets)
    achieved <- vapply(roots, function(c) mean(exp(-c * distances * units)), 0.0)
    expect_equal(achieved, targets, tolerance = 1e-11)
    expect_true(all(diff(roots) > 0))
  }
})

test_that("partial geographic extension agrees with a spherical distance oracle", {
  n <- 40L; m <- 16L; q <- 3L; seed <- 7
  set.seed(826)
  S <- cbind(runif(n, -179, 179), runif(n, -70, 70))
  selected <- .jumble_s(cbind(seq_len(n), seq_len(n)), m, seed)$coords[seq_len(m), 1]
  outside <- setdiff(seq_len(n), selected)
  S[selected[1:4], ] <- rbind(c(20, 40), c(179.999999, 45), c(10, 20), c(80, 89.99999))
  S[outside[1:4], ] <- rbind(c(-160, -40), c(-179.999999, 45), c(10.0000001, 20.0000001), c(-100, -89.99999))
  subsample <- S[selected, , drop = FALSE]
  basis <- .lnset_wc0(S, 0.1, q, 1e-4, TRUE, capN = 1L, m = m, random_t = seed)
  W0 <- .ln_subset_evecs(.getdistmat(subsample, TRUE), basis$c0, q)
  # Spherical Vincenty formula, independent of the optimized chord calculation.
  cross_distance <- vapply(seq_len(m), function(j) {
    lat1 <- S[, 2] * pi/180
    lat2 <- subsample[j, 2] * pi/180
    dl <- (S[, 1] - subsample[j, 1]) * pi/180
    y <- sqrt((cos(lat2) * sin(dl))^2 +
      (cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(dl))^2)
    x <- sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(dl)
    atan2(y, x) / (2*pi)
  }, numeric(n))
  expected <- exp(-basis$c0 * cross_distance) %*% W0
  expected <- qr.Q(qr(sweep(expected, 2L, colMeans(expected), "-")))
  actual <- basis$W[, -1, drop = FALSE]
  expect_lt(norm(tcrossprod(expected) - tcrossprod(actual), "F"), 1e-8)
})
