context("SPDE smooths and the generalised penalty")

# Bookkeeping in make_matrices --------------------------------------------

test_that("an ordinary smooth still gets one penalty and one parameter", {
  set.seed(1)
  data <- data.frame(x = stats::runif(200), z = stats::runif(200))
  mats <- hmmTMB:::make_matrices(list(a = ~ s(x, k = 6, bs = "cs")), data = data)

  expect_equal(ncol(mats$ncol_re), 1)
  expect_equal(mats$L, matrix(1, 1, 1))
  expect_equal(mats$gmrf, 0L)
  expect_equal(mats$sp_gmrf, 0L)
  expect_equal(mats$theta_start, 0)
  expect_length(mats$log_det_S, 1)
})

test_that("two smooths in one linear predictor get one determinant each", {
  set.seed(1)
  data <- data.frame(x = stats::runif(200), z = stats::runif(200))
  mats <- hmmTMB:::make_matrices(list(a = ~ s(x, k = 6, bs = "cs") + s(z, k = 5, bs = "cs")),
                        data = data)

  expect_equal(ncol(mats$ncol_re), 2)
  expect_equal(mats$gmrf, c(0L, 0L))
  # Each penalty's own generalised determinant, not the determinant of the
  # block diagonal of both. Getting this wrong is silent, and biases both
  # smoothing parameters.
  expect_length(mats$log_det_S, 2)
  q1 <- diff(mats$ncol_re[, 1]) + 1
  S1 <- as.matrix(mats$S[1:q1, 1:q1])
  expect_equal(mats$log_det_S[1], hmmTMB:::gdeterminant(S1))

  # The two blocks sit side by side, and each penalty covers its own
  # coefficients only
  expect_equal(unname(mats$ncol_re[1, 2]), unname(mats$ncol_re[2, 1]) + 1)
})

test_that("penalties are stacked across formulas", {
  set.seed(1)
  data <- data.frame(x = stats::runif(200), z = stats::runif(200))
  mats <- hmmTMB:::make_matrices(list(a = ~ s(x, k = 6, bs = "cs"),
                             b = ~ s(z, k = 5, bs = "cs")), data = data)
  expect_equal(as.numeric(mats$L), c(1, 0, 0, 1))
  expect_equal(unname(mats$ncol_re[1, 2]), unname(mats$ncol_re[2, 1]) + 1)
})

# The smoother itself ------------------------------------------------------

test_that("the SPDE smoother builds three penalties and two parameters", {
  skip_if_not_installed("fmesher")
  set.seed(1)
  data <- data.frame(x = stats::runif(300), y = stats::runif(300))
  mesh <- fmesher::fm_mesh_2d(loc = cbind(data$x, data$y),
                              max.edge = c(0.2, 0.5), cutoff = 0.1,
                              offset = c(0.1, 0.3))
  sm <- mgcv::smoothCon(mgcv::s(x, y, bs = "spde", xt = list(mesh = mesh)),
                        data = data, absorb.cons = FALSE)[[1]]

  expect_s3_class(sm, "spde.smooth")
  expect_length(sm$S, 3)
  expect_equal(dim(sm$L), c(3L, 2L))
  expect_equal(sm$theta.names, c("sd", "range"))
  expect_equal(sm$null.space.dim, 0L)     # a proper field needs no constraint
  expect_equal(nrow(sm$C), 0L)
  expect_true(sm$no.rescale)
  expect_equal(ncol(sm$X), mesh$n)
  expect_true(inherits(sm$X, "Matrix"))   # the design matrix stays sparse

  # The (sd, range) parameterisation is exactly the (tau, kappa) one: at
  # theta = (log sigma, log rho), sum_i exp((L theta)_i) S_i must reproduce
  # tau^2 (kappa^4 C + 2 kappa^2 G1 + G2)
  sigma <- 0.7
  rho <- 0.4
  kappa <- 2 * sqrt(2) / rho
  tau <- 1 / (sigma * kappa * sqrt(4 * pi))
  lambda <- exp(sm$L %*% c(log(sigma), log(rho)))
  Q1 <- lambda[1] * sm$S[[1]] + lambda[2] * sm$S[[2]] + lambda[3] * sm$S[[3]]
  fem <- fmesher::fm_fem(mesh)
  Q2 <- tau^2 * (kappa^4 * fem$c0 + 2 * kappa^2 * fem$g1 + fem$g2)
  expect_lt(max(abs(as.matrix(Q1 - Q2))) / max(abs(as.matrix(Q2))), 1e-12)
})

test_that("a one-dimensional SPDE smooth builds its own mesh", {
  skip_if_not_installed("fmesher")
  set.seed(1)
  data <- data.frame(x = stats::runif(200))
  mats <- hmmTMB:::make_matrices(list(a = ~ s(x, k = 15, bs = "spde")), data = data)

  q <- diff(mats$ncol_re[, 1]) + 1
  expect_equal(ncol(mats$ncol_re), 3)          # three penalties
  expect_equal(mats$gmrf, rep(1L, 3))
  expect_equal(dim(mats$L), c(3L, 2L))         # combined by two parameters
  expect_equal(mats$sp_gmrf, c(1L, 1L))
  expect_equal(mats$sp_names, c("a.s(x).sd", "a.s(x).range"))
  expect_equal(mats$theta_start[1], 0)         # unit marginal sd
  # All three penalties apply to the same coefficients, which is what groups
  # them into one smooth, and their blocks of S are stacked in order
  expect_equal(unname(mats$ncol_re[1, ]), rep(1, 3))
  expect_equal(unname(mats$ncol_re[2, ]), rep(q, 3))
  expect_equal(ncol(mats$S), 3 * q)
})

test_that("SPDE smooths are rejected where they cannot work", {
  skip_if_not_installed("fmesher")
  set.seed(1)
  data <- data.frame(x = stats::runif(50), y = stats::runif(50),
                     z = stats::runif(50))
  expect_error(
    mgcv::smoothCon(mgcv::s(x, y, bs = "spde"), data = data),
    "needs a mesh")
  expect_error(
    mgcv::smoothCon(mgcv::s(x, y, z, bs = "spde"), data = data),
    "one- or two-dimensional")
  expect_error(
    mgcv::smoothCon(mgcv::s(x, bs = "spde", xt = list(mesh = "not a mesh")),
                    data = data),
    "must be an fmesher mesh")
})

# End to end ---------------------------------------------------------------

## A 2-state HMM whose transition 2 -> 1 is driven by a smooth spatial field,
## observed through a Gaussian state-dependent distribution. Small and coarse
## on purpose, so that the whole thing fits in a few seconds.
sim_field <- function(n = 2000, seed = 3) {
  set.seed(seed)
  x <- y <- numeric(n)
  x[1] <- y[1] <- 0.5
  for (t in 2:n) {
    x[t] <- min(max(x[t - 1] + stats::rnorm(1, 0, 0.1), 0), 1)
    y[t] <- min(max(y[t - 1] + stats::rnorm(1, 0, 0.1), 0), 1)
  }
  u <- 1.5 * sin(4 * x) * cos(4 * y)
  s <- numeric(n)
  s[1] <- 1
  for (t in 2:n) {
    p <- if (s[t - 1] == 1) stats::plogis(-1.2) else stats::plogis(-0.5 + u[t])
    s[t] <- if (stats::runif(1) < p) 3 - s[t - 1] else s[t - 1]
  }
  data.frame(ID = 1, x = x, y = y, z = stats::rnorm(n, c(0, 4)[s], 1))
}

test_that("an SPDE field in the transition probabilities can be fitted", {
  skip_if_not_installed("fmesher")
  skip_on_cran()
  data <- sim_field()
  mesh <- fmesher::fm_mesh_2d(loc = cbind(data$x, data$y),
                              max.edge = c(0.15, 0.5), cutoff = 0.08,
                              offset = c(0.1, 0.3))
  form <- matrix(c(".", "~1",
                   "~ s(x, y, bs = 'spde', xt = list(mesh = mesh))", "."),
                 2, byrow = TRUE)
  hid <- MarkovChain$new(data = data, n_states = 2, formula = form,
                         initial_state = 1)
  obs <- Observation$new(data = data, n_states = 2, dists = list(z = "norm"),
                         par = list(z = list(mean = c(0, 4), sd = c(1, 1))))
  hmm <- HMM$new(hid = hid, obs = obs)

  # A field switches the banded forward algorithm on by itself
  expect_equal(hmm$bw(), 15L)

  # The two smoothing parameters are named, and start where the smoother says
  expect_equal(rownames(hid$lambda()),
               c("S2>S1.s(x,y).sd", "S2>S1.s(x,y).range"))
  expect_equal(unname(hid$lambda()[1, 1]), 1)

  hmm$fit(silent = TRUE)

  sp <- hmm$lambda()$hid
  expect_true(all(is.finite(sp)))
  expect_true(all(sp > 0))
  # Neither parameter has run off: a range approaching the size of the domain
  # would mean the field is unidentified
  expect_lt(sp[2, 1], 2)

  # 1/sqrt(lambda) is not a standard deviation for these two
  expect_true(all(is.na(hid$sd_re())))

  # The intercepts recover the transition probabilities that carry no field
  expect_equal(unname(hmm$hid()$coeff_fe()[1, 1]), -1.2, tolerance = 0.25)

  # And the field itself is recovered on a grid the model never saw
  grid <- expand.grid(x = seq(0.05, 0.95, length = 20),
                      y = seq(0.05, 0.95, length = 20))
  tpm <- hmm$predict("tpm", newdata = grid)
  fitted <- stats::qlogis(tpm[2, 1, ])
  truth <- -0.5 + 1.5 * sin(4 * grid$x) * cos(4 * grid$y)
  expect_gt(stats::cor(fitted, truth), 0.8)
})

test_that("a model without a field is left on the exact algorithm", {
  data <- sim_field(n = 300)
  hid <- MarkovChain$new(data = data, n_states = 2, formula = ~ s(x, k = 5, bs = "cs"),
                         initial_state = 1)
  obs <- Observation$new(data = data, n_states = 2, dists = list(z = "norm"),
                         par = list(z = list(mean = c(0, 4), sd = c(1, 1))))
  hmm <- HMM$new(hid = hid, obs = obs)
  expect_equal(hmm$bw(), 0L)
  expect_false(is.na(hid$sd_re()[1, 1]))
})
