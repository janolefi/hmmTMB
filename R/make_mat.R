
#' Create model matrices
#' 
#' @param formulas List of formulas (possibly nested, e.g. for use within Observation)
#' @param data Data frame including covariates
#' @param new_data Optional new data set, including covariates for which
#' the design matrices should be created. This needs to be passed in addition
#' to the argument '\code{data}', for cases where smooth terms or factor
#' covariates are included, and the original data set is needed to determine
#' the full range of covariate values.
#' @param gam_args Named list of additional arguments for \code{mgcv::gam()},
#' such as knots.
#' 
#' @return A list of
#' \itemize{
#'   \item X_fe Design matrix for fixed effects
#'   \item X_re Design matrix for random effects
#'   \item S Smoothness matrix
#'   \item log_det_S Vector of log-determinants of smoothness matrices
#'   \item ncol_fe Number of columns of X_fe for each parameter
#'   \item ncol_re Number of columns of X_re and S for each random effect
#'   \item L Matrix with log(lambda) = L * theta, mapping the smoothing
#'   parameters of each smooth to the weights of its penalties (mgcv's
#'   convention). It is the identity unless a smooth combines several
#'   penalties through fewer parameters, as an SPDE field does.
#'   \item gmrf For each penalty, 1 if its log-determinant depends on the
#'   parameters and so has to be computed inside the likelihood, as for an
#'   SPDE field, and 0 if the penalty is a fixed matrix scaled by exp(theta)
#'   \item sp_names Name of each smoothing parameter
#'   \item theta_start Starting value of each log smoothing parameter
#' }
#'
#' @details
#' Note the two indices: \code{ncol_re}, \code{log_det_S} and \code{gmrf} have
#' one entry per \emph{penalty}, while \code{L}, \code{sp_names} and
#' \code{theta_start} have one per \emph{smoothing parameter}. The two coincide
#' for an ordinary smooth, which has one of each, but an SPDE field has three
#' penalties and two parameters.
#' 
#' @importFrom stats update predict
#' @importFrom mgcv initial.sp
make_matrices = function(formulas, data, new_data = NULL, gam_args = NULL) {
  # Initialise lists of matrices
  X_list_fe <- list()
  X_list_re <- list()
  S_list <- list()
  ncol_fe <- NULL
  ncol_re <- NULL
  names_fe <- NULL
  names_re <- NULL
  names_ncol_re <- NULL
  L_list <- list()
  gmrf <- NULL
  sp_gmrf <- NULL
  sp_names <- NULL
  theta_start <- NULL
  log_det_S <- NULL
  start <- 1
  
  # Unlist formulas so that this function works both for Observation and MarkovChain
  forms <- unlist(formulas)
  names <- names(forms)
  
  # Loop over formulas
  for(k in seq_along(forms)) {
    form <- forms[[k]]
    
    # Check that tensor products aren't used (not supported)
    check_smooths(form)
    
    # Check that random effect variables are factors
    var_re <- find_re(form)
    for(var in var_re) {
      if(!inherits(data[[var]], "factor")) {
        data[[var]] <- factor(data[[var]])
        warning(paste0("'", var, "' is included as a random effect but is ",
                       "not a factor - changing to factor."))
      }
    }
    
    # Prepare gam() arguments
    gam_args_list <- c(list(formula = update(form, dummy_response ~ .), 
                            data = cbind(dummy_response = 1, data)),
                       gam_args)
    
    # Create matrices based on this formula
    gam_setup <- do.call(what = gam,
                         args = c(gam_args_list, list(fit = FALSE)))
    # Extract column names for design matrices
    term_names <- gam_setup$term.names
    # Starting smoothing parameter for each of this formula's penalties
    ini_sp <- initial_lambda(gam_setup)
    if(is.null(new_data)) {
      Xmat <- gam_setup$X
    } else {
      # Get design matrix for new data set. predict.gam() uses none of the
      # fitted quantities of a gam for type = "lpmatrix", so it can be given
      # the unfitted setup wrapped in a shell object. Fitting a gam to the
      # dummy response, as this used to do, is wasted work for an ordinary
      # smooth and infeasible for a mesh smooth with thousands of columns.
      Xmat <- predict(gam_shell(gam_setup), newdata = new_data,
                      type = "lpmatrix")
    }
    
    # Fixed effects design matrix
    X_list_fe[[k]] <- Xmat[, 1:gam_setup$nsdf, drop = FALSE]
    subnames_fe <- paste0(names[k], ".", term_names[1:gam_setup$nsdf])
    names_fe <- c(names_fe, subnames_fe)
    
    # Random effects design matrix
    X_list_re[[k]] <- Xmat[, -(1:gam_setup$nsdf), drop = FALSE]
    if(ncol(X_list_re[[k]]) > 0) {
      subnames_re <- paste0(names[k], ".", term_names[-(1:gam_setup$nsdf)])
      names_re <- c(names_re, subnames_re)                    
    }
    
    # Smoothing matrix
    S_list[[k]] <- bdiag_check(gam_setup$S)
    
    # Number of columns for fixed effects
    ncol_fe <- c(ncol_fe, gam_setup$nsdf)
    
    if(length(gam_setup$smooth) > 0) {
      sub_ncol_re <- matrix(1, nrow = 2, ncol = length(gam_setup$S))
      colnames(sub_ncol_re) <- 1:ncol(sub_ncol_re)
      start_s <- 1
      for (s in 1:length(gam_setup$smooth)) {
        sm <- gam_setup$smooth[[s]]
        # how many penalties for this smooth?
        npen <- length(sm$S)
        # how many parameters for this smooth? 
        npar <- ncol(sm$S[[1]])
        # where does this smooth's parameters start and end?
        sub_ncol_re[, (start_s:(start_s + npen - 1))] <- c(start, start + npar - 1)
        colnames(sub_ncol_re)[start_s:(start_s + npen - 1)] <- rep(sm$label, npen)
        # get names of smooth terms
        # regex from datascience.stackexchange.com/questions/8922
        s_terms <- gsub("(.*)\\..*", "\\1", names_re[sub_ncol_re[1, s]:sub_ncol_re[2, s]])
        s_label <- unique(s_terms)
        names_ncol_re <- c(names_ncol_re, rep(s_label, npen))
        
        # An SPDE field is proper and its precision depends on both of its
        # parameters, so its log-determinant cannot be precomputed here
        is_gmrf <- inherits(sm, "spde.smooth")
        gmrf <- c(gmrf, rep(as.integer(is_gmrf), npen))
        # Log generalised determinant of each penalty. Taking one per penalty,
        # rather than one for the block diagonal of the whole formula, is what
        # makes a linear predictor with several smooths come out right.
        log_det_S <- c(log_det_S, if(is_gmrf) rep(0, npen)
                       else unname(sapply(sm$S, gdeterminant)))
        # mgcv's L convention: several penalties may share fewer parameters
        L_sm <- if(is.null(sm$L)) diag(npen) else as.matrix(sm$L)
        L_list <- c(L_list, list(L_sm))
        ntheta <- ncol(L_sm)
        sp_gmrf <- c(sp_gmrf, rep(as.integer(is_gmrf), ntheta))
        sp_names <- c(sp_names, if(is.null(sm$theta.names)) rep(s_label, ntheta)
                      else paste0(s_label, ".", sm$theta.names))
        theta_start <- c(theta_start,
                         if(is.null(sm$theta.start))
                           log(ini_sp[start_s:(start_s + ntheta - 1)])
                         else rep(sm$theta.start, length = ntheta))
        
        start <- start + npar
        start_s <- start_s + npen
      }
      ncol_re <- cbind(ncol_re, sub_ncol_re)
    }    
  }
  colnames(ncol_re) <- names_ncol_re
  
  # Store as block diagonal matrices
  X_fe <- bdiag_check(X_list_fe)
  colnames(X_fe) <- names_fe
  X_re <- bdiag_check(X_list_re)
  colnames(X_re) <- names_re
  S <- bdiag_check(S_list)
  L <- bdiag_check(L_list)
  if(!is.null(L)) L <- as.matrix(L)
  
  return(list(X_fe = X_fe, 
              X_re = X_re, 
              S = S,
              log_det_S = log_det_S,
              X_list_fe = X_list_fe, 
              X_list_re = X_list_re, 
              S_list = S_list, 
              ncol_fe = ncol_fe, 
              ncol_re = ncol_re,
              L = L,
              # Zero-length rather than NULL when there are no smooths, so
              # that callers can use these without checking
              gmrf = as.integer(gmrf),
              sp_gmrf = as.integer(sp_gmrf),
              sp_names = as.character(sp_names),
              theta_start = as.numeric(theta_start)))
}

#' Shell gam object for building prediction matrices
#' 
#' \code{mgcv::predict.gam()} with \code{type = "lpmatrix"} uses only the model
#' frame, the terms objects, the factor levels and contrasts, and the smooth
#' objects -- all of which \code{mgcv::gam(fit = FALSE)} already produces. This
#' wraps that unfitted setup in something \code{predict.gam()} accepts, so that
#' a prediction matrix can be built without fitting anything.
#' 
#' @param G Output of \code{mgcv::gam()} called with \code{fit = FALSE}
#' 
#' @return An object of class "gam", only usable for
#' \code{predict(type = "lpmatrix")}
gam_shell <- function(G) {
  shell <- G[c("pterms", "terms", "smooth", "nsdf", "assign", "xlevels", 
               "contrasts", "pred.formula")]
  shell$model <- G$mf
  shell$na.action <- attr(G$mf, "na.action")
  shell$coefficients <- stats::setNames(rep(0, ncol(G$X)), G$term.names)
  class(shell) <- c("gam", "glm", "lm")
  return(shell)
}

#' Starting smoothing parameters for the smooths of one formula
#' 
#' Every smooth used to start at \code{lambda = 1}, whatever its basis, its
#' covariate or the size of the data. \code{mgcv::initial.sp()} instead works a
#' starting value out from the scaling of the design matrix against each
#' penalty, which is both much larger -- a few hundred for a typical smooth --
#' and basis-dependent, giving roughly 480 for a cubic regression spline, 290
#' for a P-spline and 36 for a thin plate spline on the same covariate.
#' 
#' mgcv's value is then scaled by \code{factor}, and the scaling is downwards.
#' \code{initial.sp()} is calibrated against a Gaussian penalized least squares
#' data term, and against an HMM likelihood it errs a long way to the smooth
#' side: over the models used to calibrate this, the optima sat at lambda of
#' roughly 15 (cubic regression spline), 3 (P-spline) and 0.04 (thin plate),
#' against initial.sp values of 478, 287 and 36.
#' 
#' \code{factor} is 1, that is, mgcv's value is taken as it stands, and the
#' evidence points both ways on moving it. Scaling it up helps where an
#' over-flexible start is the problem. On a two-state Gaussian model of the
#' \code{MSwM} energy data with \code{s(Oil)} on both the mean and the standard
#' deviation, starting at lambda = 1 gave a false convergence from
#' \code{nlminb} after 384 seconds, while mgcv's value converged cleanly in
#' 176 -- and scaling it by 3 or by 10 then changed nothing at all, reaching
#' the same optimum and the same smoothing parameters.
#' 
#' Scaling it up hurts elsewhere. For a basis with a null space the gradient of
#' the marginal likelihood in log(lambda) shrinks as lambda grows -- for a thin
#' plate spline it fell from 6.1 at lambda = 1 to 0.28 at ten times initial.sp
#' -- so an over-smoothed start can leave the outer optimiser on a flat surface
#' far from the optimum, where it stalls or runs lambda off to infinity. Over
#' 30 simulated models the number failing to converge was 0 at factor 0.1 and
#' 0.3, 1 at factor 1, 3 at factor 3 and 7 at factor 10, and every one of those
#' failures was a smooth on the transition probabilities rather than on a
#' state-dependent parameter.
#' 
#' So a factor of 1 is where the two sets of evidence meet: enough to fix the
#' model that needed fixing, and no further, since going further bought nothing
#' there and cost convergence elsewhere. Use
#' \code{MarkovChain$update_lambda()} or \code{Observation$update_lambda()} to
#' override it.
#' 
#' Smooths that carry their own \code{theta.start} are skipped, and so are
#' their penalties. An SPDE field is the case in hand: its parameters are a
#' marginal standard deviation and a range rather than a smoothing parameter,
#' so mgcv's value would not mean anything, and its penalties are sparse, which
#' \code{initial.sp()} does not handle.
#' 
#' @param G Output of \code{mgcv::gam()} called with \code{fit = FALSE}
#' @param factor Multiplier on mgcv's values; see 'Details' for the calibration
#' 
#' @return One starting smoothing parameter for each penalty in \code{G$S}
initial_lambda <- function(G, factor = 1) {
  sp <- rep(1, length(G$S))
  if(!length(sp)) {
    return(sp)
  }
  
  # Which penalties belong to a smooth that supplies its own starting values?
  keep <- rep(TRUE, length(sp))
  i <- 1
  for(sm in G$smooth) {
    if(!is.null(sm$theta.start)) {
      keep[i:(i + length(sm$S) - 1)] <- FALSE
    }
    i <- i + length(sm$S)
  }
  if(!any(keep)) {
    return(sp)
  }
  
  ini <- try(initial.sp(G$X, G$S[keep], G$off[keep]), silent = TRUE)
  if(!inherits(ini, "try-error")) {
    sp[keep] <- factor * ini
  }
  return(sp)
}
