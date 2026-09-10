context("Banded forward algorithm")

## A reference forward algorithm in plain R, exact and banded, written the way
## LaMa::forward_g() writes it. Having the recursion twice, once in R and once
## in C++, is the point: the C++ version is the one that has to be right, and
## it has nothing to be checked against otherwise. Note the indexing
## convention: hmmTMB's array of transition probability matrices is indexed by
## the time step being left, so the transition into t uses slice t - 1.
r_forward <- function(delta, Gamma, allprobs, bw = 0) {
  n <- nrow(allprobs)
  N <- ncol(allprobs)
  step <- function(phi, t) (phi %*% Gamma[, , t - 1]) * allprobs[t, ]

  if (bw < 2) {
    phi <- delta * allprobs[1, ]
    llk <- log(sum(phi))
    phi <- phi / sum(phi)
    for (t in seq_len(n)[-1]) {
      phi <- step(phi, t)
      llk <- llk + log(sum(phi))
      phi <- phi / sum(phi)
    }
    return(llk)
  }

  # First block: exact, from the model's own initial distribution
  phi <- delta * allprobs[1, ]
  llk <- log(sum(phi))
  phi <- phi / sum(phi)
  for (t in seq_len(min(bw, n))[-1]) {
    phi <- step(phi, t)
    llk <- llk + log(sum(phi))
    phi <- phi / sum(phi)
  }

  # Remaining blocks: warmed up from a uniform vector over the previous block
  rho <- rep(1 / N, N)
  b_start <- bw + 1
  while (b_start <= n) {
    phi <- rho * allprobs[b_start - bw, ]
    phi <- phi / sum(phi)
    for (t in seq_len(b_start - 1)[-seq_len(b_start - bw)]) {
      phi <- step(phi, t)
      phi <- phi / sum(phi)
    }
    for (t in b_start:min(b_start + bw - 1, n)) {
      phi <- step(phi, t)
      llk <- llk + log(sum(phi))
      phi <- phi / sum(phi)
    }
    b_start <- b_start + bw
  }
  llk
}

## A 2-state Gaussian HMM with a covariate on the transition probabilities.
## Two things about this design are deliberate. The chain is persistent, and
## the two state-dependent distributions overlap heavily. Either a fast-mixing
## chain or well-separated states would let the forward variable forget its
## initial condition within a couple of time steps, the banded approximation
## would be exact to machine precision at any useful bandwidth, and the tests
## below would pass without testing anything.
## 'beta_tv' is the coefficient vector hmmTMB orders as
## (1>2 intercept, 1>2 slope, 2>1 intercept, 2>1 slope).
beta_tv <- c(-3.5, 0.8, -3.5, -0.8)
mu_state <- c(0, 1.2)

sim_hmm <- function(n = 300, seed = 1) {
  set.seed(seed)
  x <- cumsum(stats::rnorm(n, 0, 0.1))
  eta <- cbind(beta_tv[1] + beta_tv[2] * x, beta_tv[3] + beta_tv[4] * x)
  s <- numeric(n)
  s[1] <- 1
  for (t in 2:n) {
    p <- stats::plogis(eta[t, s[t - 1]])
    s[t] <- ifelse(stats::runif(1) < p, 3 - s[t - 1], s[t - 1])
  }
  # Column must not be called 'state': hmmTMB reads such a column as known states
  data.frame(ID = 1, x = x, y = stats::rnorm(n, mu_state[s], 1), true_state = s)
}

## Model set-up at fixed, time-varying parameter values
make_hid <- function(data) {
  hid <- MarkovChain$new(data = data, n_states = 2, formula = ~x,
                         initial_state = 1)
  hid$update_coeff_fe(beta_tv)
  hid
}

fit_setup <- function(data, bw = NULL) {
  obs <- Observation$new(data = data, n_states = 2,
                         dists = list(y = "norm"),
                         par = list(y = list(mean = mu_state, sd = c(1, 1))))
  hmm <- HMM$new(hid = make_hid(data), obs = obs, bw = bw)
  hmm$setup()
  hmm
}

allprobs_of <- function(data) {
  cbind(stats::dnorm(data$y, mu_state[1], 1), stats::dnorm(data$y, mu_state[2], 1))
}

test_that("the exact forward algorithm is unchanged", {
  data <- sim_hmm()
  hmm <- fit_setup(data)
  expect_equal(hmm$bw(), 0L)

  # Rebuild the same likelihood in R from the model's own quantities
  Gamma <- make_hid(data)$tpm(t = "all")
  expect_gt(stats::sd(Gamma[1, 2, ]), 0.01)   # the tpm really does vary

  expect_equal(hmm$llk(), r_forward(c(1, 0), Gamma, allprobs_of(data)),
               tolerance = 1e-8)
})

test_that("the banded forward algorithm matches its R reference", {
  data <- sim_hmm()
  Gamma <- make_hid(data)$tpm(t = "all")
  allprobs <- allprobs_of(data)

  for (bw in c(2, 5, 17, 100, 400)) {
    hmm <- fit_setup(data, bw = bw)
    expect_equal(hmm$bw(), as.integer(bw))
    expect_equal(hmm$llk(), r_forward(c(1, 0), Gamma, allprobs, bw = bw),
                 tolerance = 1e-8,
                 label = paste("banded log-likelihood at bw =", bw))
  }
})

test_that("the banded log-likelihood converges to the exact one", {
  data <- sim_hmm(n = 400)
  exact <- fit_setup(data)$llk()
  bws <- c(2, 4, 8, 16, 32)
  err <- sapply(bws, function(b) abs(fit_setup(data, bw = b)$llk() - exact))

  # Geometric decay in the bandwidth
  expect_true(all(diff(err) < 0))
  expect_gt(err[1], 1)                     # the approximation really does bite
  expect_lt(err[length(err)], 1e-4)
  expect_lt(err[length(err)], err[1] / 1e4)
})

test_that("banding respects time series boundaries", {
  # Two short series, each shorter than the bandwidth, so the banded algorithm
  # should fall back to the exact one on both and reinitialise in between
  data <- rbind(sim_hmm(n = 30, seed = 2), sim_hmm(n = 30, seed = 3))
  data$ID <- rep(1:2, each = 30)

  expect_equal(fit_setup(data, bw = 50)$llk(), fit_setup(data)$llk(),
               tolerance = 1e-10)

  # And a bandwidth that does bite still gets close, on both series
  expect_equal(fit_setup(data, bw = 10)$llk(), fit_setup(data)$llk(),
               tolerance = 1e-1)
})

test_that("the banded forward algorithm agrees with LaMa", {
  skip_if_not_installed("LaMa")
  skip_if_not_installed("RTMB")
  data <- sim_hmm(n = 300)
  Gamma <- make_hid(data)$tpm(t = "all")
  allprobs <- allprobs_of(data)
  # LaMa indexes its transition matrices by the time step being entered, so
  # slice t of its array is slice t - 1 of hmmTMB's
  Gamma_lama <- Gamma[, , c(1, seq_len(dim(Gamma)[3] - 1)), drop = FALSE]

  for (bw in c(15, 30)) {
    lama <- LaMa::forward_g(c(1, 0), Gamma_lama, allprobs, bw = bw, ad = TRUE,
                            report = FALSE)
    expect_equal(fit_setup(data, bw = bw)$llk(), as.numeric(lama),
                 tolerance = 1e-8,
                 label = paste("hmmTMB vs LaMa at bw =", bw))
  }
})

test_that("update_bw validates its argument and invalidates the setup", {
  data <- sim_hmm(n = 100)
  hmm <- fit_setup(data)
  expect_error(hmm$update_bw(1), "integer >= 2")
  expect_error(hmm$update_bw(-3), "integer >= 2")
  expect_error(hmm$update_bw(2.5), "integer >= 2")

  hmm$update_bw(10)
  expect_equal(hmm$bw(), 10L)
  # The bandwidth is TMB data, so the old object cannot be reused
  expect_error(hmm$tmb_obj(), "Setup or fit model first")
  hmm$setup()
  expect_equal(hmm$bw(), 10L)

  # NULL means "decide at setup", and without a field that means exact
  hmm$update_bw(NULL)
  expect_equal(hmm$bw(), 0L)
})

test_that("check_bw profiles the log-likelihood", {
  data <- sim_hmm(n = 300)
  hmm <- fit_setup(data, bw = 15)
  prof <- hmm$check_bw(bws = c(2, 5, 20))

  expect_equal(nrow(prof), 4)
  expect_true(is.infinite(prof$bw[4]))
  # The last row is the exact algorithm, which the profile should approach
  expect_lt(abs(prof$llk[3] - prof$llk[4]), abs(prof$llk[1] - prof$llk[4]))
  expect_equal(prof$llk[4], fit_setup(data)$llk(), tolerance = 1e-10)
})
