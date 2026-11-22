context("REML")

# Two-state HMM with a smooth in the observation model and in the hidden
# state model, used by the tests below
make_reml_hmm <- function(n = 1000) {
  set.seed(1)
  data <- data.frame(x = seq(-3, 3, length.out = n), z = NA)
  s <- 1
  for(i in 1:n) {
    data$z[i] <- stats::rnorm(1, c(0, 4)[s] + sin(data$x[i]), 0.7)
    p_stay <- stats::plogis(2 + cos(data$x[i]))
    s <- sample(c(s, 3 - s), size = 1, prob = c(p_stay, 1 - p_stay))
  }
  obs <- Observation$new(data = data, dists = list(z = "norm"), n_states = 2,
                         par = list(z = list(mean = c(0, 4), sd = c(1, 1))),
                         formulas = list(z = list(mean = ~ s(x, k = 6, bs = "tp"),
                                                  sd = ~1)))
  hid <- MarkovChain$new(n_states = 2, data = data,
                         formula = ~ s(x, k = 5, bs = "cr"))
  HMM$new(obs = obs, hid = hid)
}

test_that("REML smoothing parameters match mgcv", {
  # With every state known and the Markov chain fixed, the likelihood splits
  # into one Gaussian smooth per state, for which the Laplace approximation
  # is exact. With the sd fixed at mgcv's estimate, REML must then give
  # lambda = sp / sig2 for each state. Both states share the same x, so that
  # the basis built from all rows is the one mgcv builds from each state's.
  set.seed(2)
  x0 <- stats::runif(200, -3, 3)
  data <- data.frame(x = c(x0, x0), state = rep(1:2, each = 200))
  mu <- ifelse(data$state == 1, sin(data$x) + 0.3 * data$x, 3 + cos(data$x))
  data$z <- mu + stats::rnorm(400, 0, 0.5)

  g <- lapply(1:2, function(s) {
    mgcv::gam(z ~ s(x, k = 10, bs = "tp"), data = data[data$state == s, ],
              method = "REML")
  })
  sig2 <- sapply(g, function(m) m$sig2)

  obs <- Observation$new(data = data, dists = list(z = "norm"), n_states = 2,
                         par = list(z = list(mean = c(0, 3), sd = sqrt(sig2))),
                         formulas = list(z = list(mean = ~ s(x, k = 10, bs = "tp"),
                                                  sd = ~1)))
  hid <- MarkovChain$new(n_states = 2, data = data, initial_state = 1)
  fix_hid <- stats::setNames(rep(NA, 2), rownames(hid$coeff_fe()))
  fix_obs <- c("z.sd.state1.(Intercept)" = NA, "z.sd.state2.(Intercept)" = NA)
  hmm <- suppressMessages(
    HMM$new(obs = obs, hid = hid, fixpar = list(obs = fix_obs, hid = fix_hid)))
  hmm$fit(silent = TRUE, method = "REML")

  expect_equal(unname(as.vector(hmm$lambda()$obs)),
               unname(sapply(g, function(m) m$sp / m$sig2)), tolerance = 1e-4)

  # Coefficients of each state: intercept, then the smooth
  fe <- hmm$coeff_fe()$obs[, 1]
  re <- hmm$coeff_re()$obs[, 1]
  expect_equal(unname(c(fe[1], re[1:9])), unname(coef(g[[1]])), tolerance = 1e-4)
  expect_equal(unname(c(fe[2], re[10:18])), unname(coef(g[[2]])), tolerance = 1e-4)
})

test_that("REML fits and every downstream method runs", {
  hmm <- make_reml_hmm()
  hmm$fit(silent = TRUE, method = "REML")
  expect_equal(hmm$method(), "REML")
  expect_equal(hmm$out()$convergence, 0)

  # Only the smoothing parameters and delta0 are outer parameters
  rep <- hmm$tmb_rep()
  expect_setequal(unique(names(rep$par.fixed)),
                  c("log_lambda_obs", "log_lambda_hid", "log_delta0"))

  # par_all() puts everything back in template order, which is the order of
  # jointPrecision and of coeff_array(). (coeff_array() stores lambda rather
  # than log(lambda), so only the coefficients are compared by value.)
  par <- par_all(rep)
  expect_equal(names(par), rownames(rep$jointPrecision))
  arr <- hmm$coeff_array()[!is.na(hmm$coeff_array()[, "fixed"]), ]
  expect_equal(names(par), rownames(arr))
  is_coeff <- grepl("^coeff", names(par))
  expect_equal(unname(par[is_coeff]), unname(arr[is_coeff, "value"]))

  # Confidence intervals are finite and positive-width
  ci <- hmm$confint()
  expect_true(all(is.finite(ci$coeff_fe$obs)))
  expect_true(all(ci$coeff_fe$obs[, "se"] > 0))

  # Posterior draws are centred on the estimates
  set.seed(3)
  post <- hmm$post_coeff(2000)
  fe_obs <- post[, colnames(post) == "coeff_fe_obs"]
  z <- (colMeans(fe_obs) - hmm$coeff_fe()$obs[, 1]) / ci$coeff_fe$obs[, "se"]
  expect_true(all(abs(z) < 0.2))

  # Predictions with intervals bracket the point estimate
  pr <- hmm$predict("tpm", t = c(1, 500, 1000), n_post = 200)
  expect_true(all(pr$lcl <= pr$mle + 1e-8 & pr$mle <= pr$ucl + 1e-8))

  # Model comparison: conditional quantities are defined, the marginal AIC isn't
  edf <- suppressWarnings(hmm$edf())
  expect_true(is.finite(edf) && edf > 0)
  expect_true(is.finite(suppressWarnings(hmm$AIC_conditional())))
  expect_true(is.finite(suppressWarnings(as.numeric(logLik(hmm)))))
  expect_error(hmm$AIC_marginal(), "REML")
})

test_that("REML and ML give similar fits, and fit() switches between them", {
  hmm <- make_reml_hmm()
  hmm$fit(silent = TRUE, method = "REML")
  fe_reml <- hmm$coeff_fe()$obs[, 1]

  hmm$fit(silent = TRUE, method = "ML")
  expect_equal(hmm$method(), "ML")
  expect_true("coeff_fe_obs" %in% names(hmm$tmb_rep()$par.fixed))
  expect_equal(fe_reml, hmm$coeff_fe()$obs[, 1], tolerance = 0.05)
})

test_that("REML falls back to ML without smooths", {
  set.seed(1)
  data <- data.frame(z = c(stats::rnorm(100, 0, 1), stats::rnorm(100, 4, 1)))
  obs <- Observation$new(data = data, dists = list(z = "norm"), n_states = 2,
                         par = list(z = list(mean = c(0, 4), sd = c(1, 1))))
  hid <- MarkovChain$new(n_states = 2, data = data)
  hmm <- HMM$new(obs = obs, hid = hid)
  expect_message(hmm$fit(silent = TRUE, method = "REML"), "fitting by ML")
  expect_equal(hmm$method(), "ML")
})
