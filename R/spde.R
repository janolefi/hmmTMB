## A Matern SPDE field as an mgcv smooth. The construction is Lindgren, Rue
## and Lindstrom (2011): a Gaussian field with a Matern covariance is the
## solution of a stochastic PDE, and a finite element approximation on a
## triangulation turns that into a GMRF whose precision is sparse and known in
## closed form. 'fmesher' builds the mesh and the finite element matrices;
## nothing here needs INLA.
##
## Formally the field is just another smoother, which is why it can be written
## as an mgcv smooth and dropped into any hmmTMB formula. It differs from the
## smooths hmmTMB already handles in two ways: it carries three penalty
## matrices combined through two parameters rather than a single lambda * S,
## and its precision is proper and parameter-dependent, so its log-determinant
## is computed inside the likelihood rather than precomputed. Both are handled
## by make_matrices() and by smooth_penalty() in src/likelihood.hpp.

#' Matern SPDE smooth
#'
#' A smooth term for a spatially (or temporally) continuous Gaussian random
#' field, usable as \code{s(x, y, bs = "spde", xt = list(mesh = mesh))} in any
#' hmmTMB formula, and so on any transition probability or any parameter of any
#' observation distribution.
#'
#' @section Building a mesh:
#' The mesh is yours to build and pass in, because a good one depends on the
#' domain, the data and the range you expect. \code{fmesher::fm_mesh_2d()}
#' takes locations, a maximum triangle edge length inside the domain and in
#' the outer extension, and a cutoff below which nearby points are merged:
#'
#' \preformatted{
#' mesh <- fmesher::fm_mesh_2d(loc = cbind(data$x, data$y),
#'                             max.edge = c(0.1, 0.3),
#'                             cutoff = 0.05, offset = c(0.1, 0.3))
#' f <- ~ s(x, y, bs = "spde", xt = list(mesh = mesh))
#' }
#'
#' A rule of thumb is \code{max.edge} no larger than a third of the range you
#' expect, and an offset about equal to the range; without the outer extension
#' the field's variance is inflated at the boundary. For one dimension pass a
#' \code{fmesher::fm_mesh_1d()}, or leave \code{xt} out and a mesh with
#' \code{k} evenly spaced knots is built over the range of the covariate.
#'
#' The field is evaluated at every time step, so missing covariates are
#' replaced by the last non-missing value as elsewhere in hmmTMB. If that is
#' too crude for locations, interpolate them before fitting.
#'
#' @section Parameterisation:
#' The precision is
#' \deqn{Q(\tau, \kappa) = \tau^2(\kappa^4 C + 2\kappa^2 G_1 + G_2)}
#' with \eqn{C}, \eqn{G_1}, \eqn{G_2} the finite element matrices from
#' \code{fmesher::fm_fem()}. For \code{alpha = 2} in two dimensions the Matern
#' smoothness is \eqn{\nu = 1}, the range is
#' \eqn{\rho = \sqrt{8\nu}/\kappa = 2\sqrt2/\kappa} and the marginal standard
#' deviation is \eqn{\sigma = 1/(\tau\kappa\sqrt{4\pi})}.
#'
#' The two estimated parameters are \eqn{\log\sigma} and \eqn{\log\rho}, not
#' \eqn{\log\tau} and \eqn{\log\kappa}, and appear as \code{sd} and
#' \code{range} in \code{lambda()}. The change is exact -- a linear map of the
#' log parameters, absorbed into the constants in front of the three matrices
#' -- and worth making: the parameters mean something on their own, and the
#' likelihood's long ridge along constant marginal variance runs diagonally in
#' \eqn{(\tau, \kappa)} but along an axis in \eqn{(\sigma, \rho)}, which the
#' optimiser finds much easier.
#'
#' Both are estimated freely, with no prior. A single realisation identifies
#' them only loosely, and a field whose range approaches the size of the domain
#' is close to improper: the range then runs off and the standard deviation
#' follows it, with the fitted surface barely changing. The fit is still usable
#' if that happens, but the two parameters are not, and a finer mesh will not
#' help -- it is the data that are uninformative.
#'
#' The field is proper for any positive \eqn{\kappa}, so no identifiability
#' constraint is imposed, as in INLA. It is not centred, and its level is only
#' weakly separated from the intercept when the range is large.
#'
#' @section Bandwidth:
#' A model containing this smooth is fitted with the banded forward algorithm,
#' because the exact one makes the Hessian with respect to the field weights
#' dense and so defeats the sparsity the SPDE representation exists to provide.
#' See \code{HMM$update_bw()} and \code{HMM$check_bw()}.
#'
#' @param object,data,knots As \code{mgcv::smooth.construct()}; for
#'   \code{Predict.matrix}, \code{knots} is absent and the rest are as
#'   \code{mgcv::Predict.matrix()}.
#'
#' @return A \code{smoothCon} object of class \code{spde.smooth}, with sparse
#'   \code{X} and sparse penalties.
#'
#' @references
#' Lindgren, F., Rue, H. and Lindstrom, J. (2011). An explicit link between
#' Gaussian fields and Gaussian Markov random fields. \emph{JRSS-B} 73, 423-498.
#'
#' Fischer, J.-O. (2026). Fast and scalable inference in hidden Markov models
#' with Gaussian fields.
#'
#' @examples
#' # A one-dimensional field over a covariate, with a default mesh
#' if(requireNamespace("fmesher", quietly = TRUE)) {
#'   d <- data.frame(x = runif(100))
#'   sm <- mgcv::smoothCon(mgcv::s(x, bs = "spde", k = 20), data = d)[[1]]
#'   sm$theta.names
#' }
#'
#' @exportS3Method mgcv::smooth.construct
smooth.construct.spde.smooth.spec <- function(object, data, knots) {
  if(!requireNamespace("fmesher", quietly = TRUE)) {
    stop("bs = \"spde\" needs the fmesher package", call. = FALSE)
  }
  if(length(object$term) > 2 | length(object$term) < 1) {
    stop("An SPDE smooth is one- or two-dimensional, but this one has ",
         length(object$term), " terms.", call. = FALSE)
  }

  mesh <- object$xt$mesh
  if(is.null(mesh)) {
    if(length(object$term) == 2) {
      stop("A two-dimensional SPDE smooth needs a mesh: build one with ",
           "fmesher::fm_mesh_2d() and pass it as ",
           "s(x, y, bs = \"spde\", xt = list(mesh = mesh)). The right mesh ",
           "depends on the domain and on the range you expect, so there is ",
           "no useful default.", call. = FALSE)
    }
    k <- ifelse(object$bs.dim < 0, 20, object$bs.dim)
    x <- data[[object$term]]
    mesh <- fmesher::fm_mesh_1d(seq(min(x), max(x), length.out = k),
                                degree = 2, boundary = "free")
  }
  if(!inherits(mesh, c("fm_mesh_1d", "fm_mesh_2d", "inla.mesh", "inla.mesh.1d"))) {
    stop("xt$mesh must be an fmesher mesh, from fm_mesh_2d() or fm_mesh_1d().",
         call. = FALSE)
  }

  object$X <- fmesher::fm_basis(mesh, spde_loc(object$term, data))
  fem <- fmesher::fm_fem(mesh)
  # The mass matrix is the lumped one, as in INLA: fm_fem already forms
  # g2 = g1 C0^-1 g1 with it, and using the consistent c1 in the kappa^4 term
  # would make the three matrices inconsistent with each other.
  # lambda = exp(L theta) with theta = (log sigma, log rho). Substituting
  # kappa = 2 sqrt(2) / rho and tau = 1 / (sigma kappa sqrt(4 pi)) into
  # tau^2 (kappa^4 C + 2 kappa^2 G1 + G2) leaves these constants in front of
  # the three matrices, and these powers of sigma and rho.
  object$S <- list((2/pi) * fem$c0, fem$g1 / (2*pi), fem$g2 / (32*pi))
  object$L <- matrix(c(-2, -2, -2, -2, 0, 2), ncol = 2)
  object$theta.names <- c("sd", "range")
  # A range of a fifth of the domain and unit marginal standard deviation:
  # smooth without being flat, on the scale of a linear predictor whose other
  # terms are of order one.
  object$theta.start <- c(0, log(spde_extent(mesh, object$X)/5))

  object$rank <- rep(ncol(object$X), length(object$S))
  object$null.space.dim <- 0    # proper for any positive kappa
  # Two flags that keep smoothCon out of the way. 'no.rescale' suppresses its
  # penalty rescaling, which would divide the finite element matrices by a norm
  # of X and so break the meaning of tau and kappa -- and which is also the only
  # place it takes a matrix norm of X, so the design matrix can stay sparse. A
  # zero-row C says there is no centring constraint, which is right for a proper
  # field and is what INLA does.
  object$no.rescale <- TRUE
  object$C <- matrix(0, 0, ncol(object$X))
  object$df <- ncol(object$X)
  object$mesh <- mesh
  object$te.ok <- 0
  class(object) <- "spde.smooth"
  return(object)
}

#' @rdname smooth.construct.spde.smooth.spec
#' @exportS3Method mgcv::Predict.matrix
Predict.matrix.spde.smooth <- function(object, data) {
  # Dense, unlike the design matrix built above: this is only reached through
  # predict.gam(), which splices the result into a dense matrix spanning every
  # coefficient in the model and so would densify it anyway.
  as.matrix(fmesher::fm_basis(object$mesh, spde_loc(object$term, data)))
}

#' Observation locations for an SPDE smooth
#'
#' @param term Character vector of one or two covariate names
#' @param data Data frame containing those covariates
#'
#' @return Numeric vector (1d) or two-column matrix (2d) of locations
spde_loc <- function(term, data) {
  if(length(term) == 1) {
    return(as.numeric(data[[term]]))
  }
  return(cbind(as.numeric(data[[term[1]]]), as.numeric(data[[term[2]]])))
}

#' A length scale for the meshed domain, for starting values
#'
#' The larger side of the bounding box of the mesh nodes that carry basis
#' weight, which is the region the data occupy rather than the outer extension.
#' A one-dimensional mesh stores its knots as a plain vector, and its degree-2
#' basis has one more function than knots, so both the shape and the length of
#' 'loc' have to be taken as they come.
#'
#' @param mesh An fmesher mesh
#' @param X The smooth's design matrix
#'
#' @return A positive number
spde_extent <- function(mesh, X) {
  loc <- as.matrix(if(is.null(mesh$loc)) mesh$mid else mesh$loc)
  loc <- loc[, seq_len(min(2, ncol(loc))), drop = FALSE]
  used <- which(Matrix::colSums(abs(X)) > 0)
  if(length(used) > 1 & max(used) <= nrow(loc)) {
    loc <- loc[used, , drop = FALSE]
  }
  return(max(apply(loc, 2, function(z) diff(range(z))), .Machine$double.eps))
}
