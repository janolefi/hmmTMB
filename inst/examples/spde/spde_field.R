### A latent Gaussian field in the transition probabilities of an HMM
###
### Simulated version of the movement-ecology example in Fischer (2026): an
### animal switches between two behaviours, and the probability of switching
### out of the active behaviour varies smoothly over space. The spatial effect
### is modelled as a Matern field, represented on a triangulated mesh through
### the SPDE approach of Lindgren, Rue and Lindstrom (2011), which turns it
### into a Gaussian Markov random field with a sparse precision matrix.
###
### The field enters the model as an ordinary mgcv smooth,
### s(x, y, bs = "spde"), so nothing about the modelling interface changes.
### What does change is how the likelihood is computed: hmmTMB switches to the
### banded forward algorithm, because the exact one would make the Hessian
### with respect to the several hundred field weights dense and the Laplace
### approximation unaffordable.
###
### Runs in well under a minute.

library(hmmTMB)
library(fmesher)


# Simulate ----------------------------------------------------------------

set.seed(3)
n <- 2000

## A track wandering over the unit square
x <- y <- numeric(n)
x[1] <- y[1] <- 0.5
for (t in 2:n) {
  x[t] <- min(max(x[t - 1] + rnorm(1, 0, 0.1), 0), 1)
  y[t] <- min(max(y[t - 1] + rnorm(1, 0, 0.1), 0), 1)
}

## The true field, on the linear predictor of Pr(state 2 -> state 1)
field <- function(x, y) 1.5 * sin(4 * x) * cos(4 * y)

## Two states: state 1 is "resting" (low mean), state 2 is "active"
states <- numeric(n)
states[1] <- 1
for (t in 2:n) {
  p <- if (states[t - 1] == 1) plogis(-1.2) else plogis(-0.5 + field(x[t], y[t]))
  states[t] <- if (runif(1) < p) 3 - states[t - 1] else states[t - 1]
}

data <- data.frame(ID = 1, x = x, y = y, z = rnorm(n, c(0, 4)[states], 1))


# Mesh --------------------------------------------------------------------

## The mesh is the one part of an SPDE model with no sensible default: it has
## to resolve the features you expect to see. A rule of thumb is a maximum
## edge length of about a third of the range you expect, plus an outer
## extension of about one range, without which the field's variance is
## inflated at the boundary.
mesh <- fm_mesh_2d(loc = cbind(data$x, data$y),
                   max.edge = c(0.15, 0.5),
                   cutoff = 0.08,
                   offset = c(0.1, 0.3))
plot(mesh, asp = 1)
points(data$x, data$y, pch = 20, cex = 0.4, col = "#00798c40")
mesh$n   # number of field weights, integrated out by the Laplace approximation


# Fit ---------------------------------------------------------------------

## The field goes on the 2 -> 1 transition only; the 1 -> 2 transition is
## intercept-only. Transition-specific formulas are given as a matrix, with
## "." on the reference (diagonal) entries.
form <- matrix(c(".", "~ 1",
                 "~ s(x, y, bs = 'spde', xt = list(mesh = mesh))", "."),
               nrow = 2, byrow = TRUE)

hid <- MarkovChain$new(data = data, n_states = 2, formula = form,
                       initial_state = "stationary")

obs <- Observation$new(data = data, n_states = 2,
                       dists = list(z = "norm"),
                       par = list(z = list(mean = c(0, 4), sd = c(1, 1))))

hmm <- HMM$new(hid = hid, obs = obs)

## The field switched the banded forward algorithm on, with the default
## bandwidth of 15
hmm$bw()

system.time(hmm$fit(silent = TRUE))


# Check the bandwidth -----------------------------------------------------

## The bandwidth controls how far back in time each log-likelihood
## contribution is allowed to look. Too small and the approximation bites;
## large enough and the log-likelihood stops moving. Profile it at the
## estimates: if the curve has flattened well before the bandwidth in use,
## the approximation is fine.
hmm$check_bw(bws = seq(5, 30, by = 5))

## If it has not, raise it and refit:
# hmm$update_bw(30)
# hmm$fit(silent = TRUE)


# Results -----------------------------------------------------------------

## The field's two parameters are a marginal standard deviation and a range,
## on the scale of the linear predictor and of the coordinates respectively.
## They are not variances, so sd_re() reports NA for them.
hmm$lambda()$hid

## Everything else works as usual
hmm$coeff_fe()$hid
hmm$obs()$par()[, , 1]

## The fitted field, on a grid the model never saw
grid <- expand.grid(x = seq(0.02, 0.98, length = 150),
                    y = seq(0.02, 0.98, length = 150))
tpm <- hmm$predict("tpm", newdata = grid)

par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
image(unique(grid$x), unique(grid$y),
      matrix(plogis(-0.5 + field(grid$x, grid$y)), 150, 150),
      col = hcl.colors(50), asp = 1, xlab = "x", ylab = "y",
      main = "True Pr(active -> resting)", zlim = c(0, 1))
image(unique(grid$x), unique(grid$y),
      matrix(tpm[2, 1, ], 150, 150),
      col = hcl.colors(50), asp = 1, xlab = "x", ylab = "y",
      main = "Estimated", zlim = c(0, 1))
par(mfrow = c(1, 1))

## Correlation between the fitted and the true surface, on the linear
## predictor scale
cor(qlogis(tpm[2, 1, ]), -0.5 + field(grid$x, grid$y))

## Decoded states
table(hmm$viterbi(), states)
