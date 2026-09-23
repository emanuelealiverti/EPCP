#' Control parameters for \code{viord}
#'
#' Collects the options that govern the two nested loops of \code{\link{viord}}:
#' the inner loop of the variational (or EP) algorithm, run at fixed thresholds,
#' and the outer loop that re-estimates the thresholds.
#'
#' @param maxit_inner Integer. Maximum number of iterations of the inner
#'   algorithm, at fixed thresholds.
#' @param maxit_outer Integer. Maximum number of threshold updates. Ignored when
#'   \code{fix_alpha = TRUE}, in which case a single inner fit is performed.
#' @param tol_inner Numeric. Convergence tolerance of the inner algorithm.
#' @param tol_outer Numeric. Convergence tolerance of the threshold loop.
#'   Both tolerances are relative: convergence is declared when the change in
#'   the objective is smaller than \code{tol * (1 + abs(objective))}, so that
#'   the same value is meaningful at any sample size.
#' @param min_iter Integer. Minimum number of iterations before convergence may
#'   be declared, in both loops. Useful when two consecutive iterations may
#'   happen to agree before the algorithm has really settled.
#' @param conv_crit Character. \code{"elbo"} monitors the objective (the ELBO,
#'   or the EP log marginal likelihood); \code{"coef"} monitors the variational
#'   parameters themselves — posterior means, thresholds and log variance
#'   components — through the largest relative change
#'   \eqn{\max_j |\Delta \theta_j| / (1 + |\theta_j|)}. The latter is useful
#'   when the objective is nearly flat along directions in which the variance
#'   components still move, which is common with large \eqn{n}. Not used by
#'   \code{algorithm = "EP"}, whose inner loop always monitors its own
#'   objective.
#' @param warm_start Logical. If \code{TRUE}, each threshold iteration starts
#'   the inner algorithm from the state reached by the previous one instead of
#'   restarting from scratch. The fixed point is unchanged, but far fewer inner
#'   iterations are needed. Not used by \code{algorithm = "EP"}.
#' @param alpha_init Numeric vector of \code{K - 1} starting thresholds, or
#'   \code{NULL} to obtain them from a covariate-free fit. Required when
#'   \code{fix_alpha = TRUE}.
#' @param fix_alpha Logical. If \code{TRUE}, the thresholds are held fixed at
#'   \code{alpha_init} and the outer loop is skipped entirely.
#' @param verbose Integer. \code{0} prints nothing, \code{1} reports the
#'   progress of the threshold loop, \code{2} also reports the inner iterations.
#' @param full_out Logical. If \code{TRUE}, the returned object also carries the
#'   quantities used for diagnostics and posterior sampling (\code{xiZ},
#'   \code{sigmaZ}, the ELBO sequence, the convergence flag).
#'
#' @return A list of validated control parameters, to be passed as the
#'   \code{control} argument of \code{\link{viord}}.
#'
#' @seealso \code{\link{viord}}
#'
#' @examples
#' \dontrun{
#' # looser threshold loop, chatty output
#' viord(Y, X, prior, algorithm = "PMF",
#'       control = viord.control(maxit_outer = 200, tol_outer = 1e-4, verbose = 1))
#'
#' # thresholds held fixed at known values
#' viord(Y, X, prior, algorithm = "PMF",
#'       control = viord.control(alpha_init = c(-0.5, 0.5), fix_alpha = TRUE))
#' }
#' @export
viord.control <- function(maxit_inner = 100,
                          maxit_outer = 100,
                          tol_inner   = 1e-6,
                          tol_outer   = 1e-6,
                          min_iter    = 1,
                          conv_crit   = c("elbo", "coef"),
                          warm_start  = TRUE,
                          alpha_init  = NULL,
                          fix_alpha   = FALSE,
                          verbose     = 0,
                          full_out    = TRUE) {

  conv_crit <- match.arg(conv_crit)

  check_count <- function(x, name, min = 1) {
    if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x < min || x != as.integer(x))
      stop(sprintf("%s must be a single integer >= %d.", name, min))
    as.integer(x)
  }
  check_tol <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x <= 0)
      stop(sprintf("%s must be a single strictly positive number.", name))
    as.numeric(x)
  }
  check_flag <- function(x, name) {
    if (!is.logical(x) || length(x) != 1 || is.na(x))
      stop(sprintf("%s must be TRUE or FALSE.", name))
    x
  }

  maxit_inner <- check_count(maxit_inner, "maxit_inner")
  maxit_outer <- check_count(maxit_outer, "maxit_outer")
  tol_inner   <- check_tol(tol_inner,   "tol_inner")
  tol_outer   <- check_tol(tol_outer,   "tol_outer")
  min_iter    <- check_count(min_iter,  "min_iter")
  warm_start  <- check_flag(warm_start, "warm_start")
  fix_alpha   <- check_flag(fix_alpha,  "fix_alpha")
  full_out    <- check_flag(full_out,   "full_out")

  verbose <- check_count(verbose, "verbose", min = 0)
  if (verbose > 2) stop("verbose must be 0, 1 or 2.")

  if (!is.null(alpha_init)) {
    if (!is.numeric(alpha_init) || any(!is.finite(alpha_init)))
      stop("alpha_init must be a numeric vector of finite values.")
    if (is.unsorted(alpha_init, strictly = TRUE))
      stop("alpha_init must be strictly increasing.")
  }
  if (fix_alpha && is.null(alpha_init))
    stop("fix_alpha = TRUE requires alpha_init.")

  list(maxit_inner = maxit_inner,
       maxit_outer = maxit_outer,
       tol_inner   = tol_inner,
       tol_outer   = tol_outer,
       min_iter    = min_iter,
       conv_crit   = conv_crit,
       warm_start  = warm_start,
       alpha_init  = alpha_init,
       fix_alpha   = fix_alpha,
       verbose     = verbose,
       full_out    = full_out)
}
