
# hmmTMB 1.1.3

- Add latent Gaussian fields, as a Matern SPDE smoother usable anywhere an
  mgcv smooth is: `s(x, y, bs = "spde", xt = list(mesh = mesh))`. See
  `?smooth.construct.spde.smooth.spec` and `inst/examples/spde/spde_field.R`.
- Add the banded forward algorithm of Fischer (2026), which is what makes a
  high-dimensional field affordable. New `bw` argument of `HMM$new()`, with
  `HMM$update_bw()` and `HMM$check_bw()`; on by default with `bw = 15` when
  the model contains a field.
- Generalise the smoothing penalty: a smooth may combine several penalty
  matrices through fewer smoothing parameters (mgcv's `L`), and may have its
  log-determinant computed inside the likelihood.
- Compile with TMB's TMBad framework by default; set
  `HMMTMB_AD_FRAMEWORK=CppAD` before installing to fall back.
- Fix the log-determinant of penalty matrices when one linear predictor
  contains several smooths, which biased their smoothing parameters.
- Build prediction matrices without fitting a throwaway `mgcv::gam()`.
- Sample from the precision matrix rather than from its inverse in
  `HMM$post_coeff()`, through a sparse Cholesky factorisation. Drawing 1000
  posterior samples from a model with a 3400-node field takes 0.6s instead of
  26s, and no longer forms a dense matrix of the same size. The precision
  itself -- TMB's joint precision with random effects, the Hessian without --
  is exposed as `HMM$post_prec()` and cached.
- `HMM$confint()` no longer inverts the joint precision. It needs only the
  covariance of the fixed effects, which is the leading block of that inverse
  and is what `TMB::sdreport()` already returns as `cov.fixed`.
- New vignette on (semi-)supervised learning
- Fix parameter counts for models with constraints
- Use safe Hessian inversion even for models without random effects in `HMM$post_coeff()` and `HMM$confint()`

# hmmTMB 1.1.2

- Add hurdle negative binomial distribution
- Fix bug in `HMM$fit_stan()`
- Fix bug in `HMM$pseudores()`

# hmmTMB 1.1.1

- Add error message for `laplace = TRUE` in `HMM$fit_stan()`
- Add reference to JSS paper

# hmmTMB 1.1.0

- Add standard error to output of `HMM$confint()`
- Use ID-specific `delta0` when `initial_state = "stationary"`
- Add `par_alt()` functions for pretty display of cat and mvn parameters.
- Improve mvn distribution (automatically detect dimension), and fix mvnorm for dim > 2, using a Cholesky decomposition to get unconstrained parameters of covariance matrix.
- Add zero-one-inflated beta distribution
- Replace optimx by nlminb for model fitting
- Allow empty models for simulation
- Use generalized determinant for penalty matrices 
- Add nbinom(mean, shape) distribution
- Allow for user-specified arguments to be passed directly to mgcv::gam(); e.g., knots.
- Add built-in functions dvm, rvm, dwrpcauchy, rwrpcauchy, to remove dependence on CircStats (as requested by CRAN).

# hmmTMB 1.0.2

- Fix bug for categorical distribution
- Fix bug when using tibble input data
- Fix bugs for models with shared parameters
- Add `HMM$suggest_initial()`
