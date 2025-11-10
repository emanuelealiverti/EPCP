#' Variational inference for ordinal models
#'
#' Unified interface based for approximate Bayesian inference under the cumulative probit based on Expectation Propagation (EP),
#' Mean-Field (MF), and Partially Mean-Field (PMF) algorithms.
#'
#' @param Y Response (ordinal factor)
#' @param X Design matrix
#' @param prior List with prior parameters (mu0, S0, Q0)
#' @param algorithm Character, one of "EP", "MF", or "PMF"
#' @param maxit Maximum number of iterations
#' @param conv_tr Convergence tolerance
#' @param method Optimizer for the thresholds ("grad" or "polr")
#'
#' @return A list containing the estimates and convergence information
#' @export
#'
#' @examples
#' \dontrun{
#' res <- viord(Y, X, prior, algorithm = "MF")
#' }
viord = function(Y, X, prior,
                 algorithm = c("EP", "MF", "PMF"),
                 maxit = 100, conv_tr = 1e-6) {

  algorithm <- match.arg(algorithm)

  if (algorithm == "EP") {
    out <- optim_ep_ml(Y = Y, X = X, prior = prior,
                       maxit = maxit, conv_tr = conv_tr)

  } else {
    vb_factor <- ifelse(algorithm == "MF", "MF", "PMF")
    out <- optim_vb_ml(Y = Y, X = X, prior = prior,
                       maxit = maxit, conv_tr = conv_tr,
                       vb_factor = vb_factor)
  }

  out$coef.names = colnames(X)
  out$algorithm = algorithm
  class(out) = 'viord'
  return(out)
}
