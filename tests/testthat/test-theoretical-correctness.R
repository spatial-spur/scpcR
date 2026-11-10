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
