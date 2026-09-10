context("Uncertainty quantification")

# Sampling from a precision matrix -----------------------------------------

test_that("rmvn_prec draws from N(mu, Q^-1)", {
  set.seed(1)
  # A dense precision and a sparse one, since the point of the function is
  # that the two go the same way
  dense <- crossprod(matrix(stats::rnorm(100), 10, 10)) + diag(10)
  sparse <- Matrix::bandSparse(60, k = c(0, 1), symmetric = TRUE,
                               diagonals = list(rep(2, 60), rep(-0.4, 59)))
  for (Q in list(dense = dense, sparse = sparse)) {
    q <- nrow(Q)
    mu <- seq_len(q) / q
    V <- as.matrix(Matrix::solve(Q))
    set.seed(2)
    S <- rmvn_prec(2e5, mu, Q)

    expect_equal(dim(S), c(2e5L, q))
    # Monte Carlo error at 2e5 draws is well under 2%
    expect_lt(max(abs(colMeans(S) - mu) / sqrt(diag(V))), 0.02)
    expect_lt(max(abs(apply(S, 2, stats::sd) / sqrt(diag(V)) - 1)), 0.02)
    # Off-diagonal structure too, not just the marginals
    expect_lt(max(abs(stats::cov(S[, 1:5]) - V[1:5, 1:5])) /
                max(abs(V[1:5, 1:5])), 0.05)
  }
})

test_that("rmvn_prec falls back when the precision is not positive definite", {
  set.seed(1)
  Q <- crossprod(matrix(stats::rnorm(100), 10, 10)) + diag(10)
  Q[10, ] <- Q[, 10] <- 0     # singular, so there is no Cholesky factor
  expect_warning(S <- rmvn_prec(5, rep(0, 10), Q), "generalised inverse|ginv")
  expect_equal(dim(S), c(5L, 10L))
  expect_true(all(is.finite(S)))
})

# The two branches of the model -------------------------------------------

## A smooth covariate effect on the transition probabilities, so that the
## random-effects model below is actually identified. A smooth fitted to a
## covariate that does nothing is estimated as flat, its smoothing parameter
## runs off, and the joint precision becomes computationally singular -- which
## exercises the fallback rather than the code these tests are about.
sim_data <- function(n = 1500, seed = 1) {
  set.seed(seed)
  x <- seq(-3, 3, length.out = n)
  s <- numeric(n)
  s[1] <- 1
  for (t in 2:n) {
    p <- stats::plogis(-1.5 + 1.5 * sin(x[t]))
    s[t] <- if (stats::runif(1) < p) 3 - s[t - 1] else s[t - 1]
  }
  data.frame(ID = 1, x = x, z = stats::rnorm(n, c(0, 4)[s], 1))
}

fit_hmm <- function(data, formula = ~1) {
  hmm <- HMM$new(hid = MarkovChain$new(data = data, n_states = 2,
                                       formula = formula, initial_state = 1),
                 obs = Observation$new(data = data, n_states = 2,
                                       dists = list(z = "norm"),
                                       par = list(z = list(mean = c(0, 4),
                                                           sd = c(1, 1)))))
  hmm$fit(silent = TRUE)
  hmm
}

test_that("a model without random effects samples from its Hessian", {
  hmm <- fit_hmm(sim_data())
  rep <- hmm$tmb_rep()
  expect_null(rep$jointPrecision)

  # The precision is the Hessian, not the joint precision
  expect_equal(as.matrix(hmm$post_prec()),
               hmm$tmb_obj()$he(rep$par.fixed), tolerance = 1e-10)

  set.seed(3)
  post <- hmm$post_coeff(2e4)
  free <- which(!is.na(hmm$coeff_array()[, "fixed"]))
  expect_lt(max(abs(apply(post[, free], 2, stats::sd) /
                      sqrt(diag(rep$cov.fixed)) - 1)), 0.1)
})

test_that("a model with random effects samples from the joint precision", {
  hmm <- fit_hmm(sim_data(), formula = ~ s(x, k = 6, bs = "cs"))
  rep <- hmm$tmb_rep()
  expect_false(is.null(rep$jointPrecision))
  expect_identical(hmm$post_prec(), rep$jointPrecision)

  V <- prec_to_cov(rep$jointPrecision)
  set.seed(3)
  post <- hmm$post_coeff(2e4)
  free <- which(!is.na(hmm$coeff_array()[, "fixed"]))
  expect_equal(ncol(post[, free, drop = FALSE]), ncol(V))
  expect_lt(max(abs(apply(post[, free], 2, stats::sd) / sqrt(diag(V)) - 1)), 0.1)
})

test_that("the precision is cached and dropped when the model is refitted", {
  hmm <- fit_hmm(sim_data())
  expect_identical(hmm$post_prec(), hmm$post_prec())   # second call is cached

  # Refitting replaces the report, and the precision belonged to it. A stale
  # cache would still hold the Hessian at the old estimates, which differs.
  hmm$fit(silent = TRUE)
  expect_equal(as.matrix(hmm$post_prec()),
               hmm$tmb_obj()$he(hmm$tmb_rep()$par.fixed), tolerance = 1e-10)
})

test_that("confint needs no inversion of the joint precision", {
  hmm <- fit_hmm(sim_data(), formula = ~ s(x, k = 6, bs = "cs"))
  rep <- hmm$tmb_rep()
  # The covariance of the fixed effects is the leading block of the inverse of
  # the joint precision. That is the identity confint() relies on to use
  # cov.fixed directly, and it is worth pinning down.
  V <- prec_to_cov(rep$jointPrecision)
  np <- length(rep$par.fixed)
  expect_equal(unname(diag(V)[1:np]), unname(diag(rep$cov.fixed)),
               tolerance = 1e-8)

  ci <- hmm$confint()
  se <- sqrt(diag(rep$cov.fixed))
  expect_equal(unname(ci$coeff_fe$hid[, "se"]),
               unname(se[names(rep$par.fixed) == "coeff_fe_hid"]),
               tolerance = 1e-8)
})
