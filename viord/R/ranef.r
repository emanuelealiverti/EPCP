#' Random-effect posterior means
#'
#' Extracts the posterior means of the random-effect coefficients from a fitted
#' \code{viord} model.
#'
#' @param object A fitted object of class \code{"viord"}.
#' @param ... Currently unused.
#'
#' @return A named list with one component per variance-component group.
#'   Each component is a named numeric vector of posterior means for the
#'   random-effect coefficients belonging to that group.
#'   Returns \code{NULL} (invisibly) if the model contains no random effects.
#'
#' @seealso \code{\link{postVar.viord}} for the posterior covariance matrices,
#'   \code{\link{viord}} for model fitting.
#' @name ranef
#' @importFrom nlme ranef
#' @export
NULL

# The generic comes from nlme (a recommended package, and the same object that
# lme4 re-exports), so that ranef(fit) keeps working whichever of these packages
# is attached last. Registering the method on that generic avoids the clash a
# package-local generic would cause.
#' @exportS3Method nlme::ranef
ranef.viord <- function(object, ...) {
  if (is.null(object$u.names) || length(object$est$m_u) == 0) {
    message("No random effects in this model.")
    return(invisible(NULL))
  }

  u_hat        <- object$est$m_u
  names(u_hat) <- object$u.names
  groups       <- as.character(object$Z_group)
  group_labels <- unique(groups)

  out <- vector("list", length(group_labels))
  names(out) <- group_labels
  for (g in group_labels) {
    idx      <- which(groups == g)
    out[[g]] <- u_hat[idx]
  }
  out
}

#' Posterior covariance matrices of random effects
#'
#' Extracts the posterior covariance matrices of the random-effect coefficients
#' from a fitted \code{viord} model.
#'
#' @param object A fitted object of class \code{"viord"}.
#' @param ... Currently unused.
#'
#' @return A named list with one component per variance-component group.
#'   Each component is the posterior covariance submatrix (with row and column
#'   names from \code{colnames(Z)}) for the random-effect coefficients belonging
#'   to that group.
#'   Returns \code{NULL} (invisibly) if the model contains no random effects.
#'
#' @seealso \code{\link{ranef.viord}} for the posterior means,
#'   \code{\link{viord}} for model fitting.
#' @export
postVar <- function(object, ...) UseMethod("postVar")

#' @export
postVar.viord <- function(object, ...) {
  if (is.null(object$u.names) || length(object$est$m_u) == 0) {
    message("No random effects in this model.")
    return(invisible(NULL))
  }

  S_u                      <- object$est$S_u
  rownames(S_u) <- colnames(S_u) <- object$u.names
  groups                   <- as.character(object$Z_group)
  group_labels             <- unique(groups)

  out <- vector("list", length(group_labels))
  names(out) <- group_labels
  for (g in group_labels) {
    idx      <- which(groups == g)
    out[[g]] <- S_u[idx, idx, drop = FALSE]
  }
  out
}
