#' Approximate Bayesian Inference for Cumulative Probit Models
#'
#' Provides a unified interface for scalable approximate Bayesian inference under
#' the cumulative probit model using one of three algorithms:
#' \emph{Expectation Propagation (EP)}, \emph{Mean-Field Variational Bayes (MF)}, or
#' \emph{Partially Factorized Mean-Field (PMF)}.
#'
#' Threshold (cutoff) parameters are estimated via approximate marginal likelihood,
#' alternating between the optimization of the thresholds (via Newton–Raphson steps)
#' and the internal optimization based on EP, MF, or PMF.  
#' The core computational routines are implemented in C++ for efficiency,
#' while this function provides a high-level R interface.
#'
#' @param Y Ordinal response variable (an ordered factor).
#' @param X Design matrix of covariates.
#' @param prior A list containing prior parameters:
#'   \code{mu0} (prior mean), \code{S0} (prior covariance), and \code{Q0} (prior precision matrix).
#' @param algorithm Character string specifying the inference algorithm to use:
#'   one of \code{"EP"}, \code{"MF"}, or \code{"PMF"}.
#' @param maxit Integer specifying the maximum number of iterations used in both
#'   the alternating optimization of the thresholds and the internal optimization
#'   based on the selected approximation algorithm.
#' @param conv_tr Numeric value giving the convergence tolerance for both the
#'   threshold optimization and the internal inference loop.
#' @return A list containing the estimated model quantities and convergence information.  
#'   The output includes:
#'   \itemize{
#'     \item \code{est}: posterior (approximate) mean vector \eqn{m} and covariance matrix \eqn{S};
#'     \item \code{alpha}: estimated thresholds (cutoffs);
#'     \item \code{algorithm}: the selected inference algorithm;
#'     \item \code{prior}: the prior provided as input
#'     \item additional fields used for convergence diagnostics and summaries.
#'   }
#' 
#' @details
#' Consider an ordinal probit model for outcomes \eqn{y_i \in \{1, \dots, K\}} with
#' a \eqn{p}-dimensional covariate vector \eqn{x_i = (x_{i1}, \dots, x_{ip})}.
#' The model is specified as
#' \deqn{
#'   \Pr(y_i \leq k) = \Phi(\alpha_k - x_i^\top \beta),
#' }
#' where \eqn{\Phi(\cdot)} denotes the standard normal cumulative distribution function.
#' A Gaussian prior \eqn{\beta \sim N_p(\mu_0, S_0)} is assumed for the regression coefficients.
#' The thresholds \eqn{\alpha_1 < \dots < \alpha_{K-1}} are treated as nuisance parameters
#' and estimated by maximizing the (approximate) marginal likelihood.
#'
#' The \code{algorithm} argument controls the internal approximation used
#' for inference conditional on the current thresholds:
#' \itemize{
#'   \item \code{"EP"} – Expectation Propagation;
#'   \item \code{"MF"} – Mean-Field Variational Bayes;
#'   \item \code{"PMF"} – Partially Factorized Mean-Field Variational Bayes.
#' }
#'
#' @references
#' Aliverti, E. (2025).
#' *Approximate Bayesian Inference for Cumulative Probit Regression Models*.
#' \url{https://arxiv.org/abs/2511.06967}
#' @seealso
#' \code{\link{summary.viord}} for model summaries,
#' \code{\link{predict.viord}} for posterior predictive probabilities, and
#' \code{\link{simulate.viord}} for posterior sampling from the approximate posterior
#'
#' @examples
#' \dontrun{
#' # See the tutorials folder at github.com/emanuelealiverti/EPCP for more examples
#' library(viord)
#'
#' # Simulate a small ordinal dataset
#' set.seed(1)
#' n <- 200
#' X <- cbind(1, matrix(rnorm(n * 2), n, 2))
#' beta <- c(0.5, -1, 1.2)
#' alpha <- c(-0.5, 0.5)
#' z <- X %*% beta + rnorm(n)
#' Y <- cut(z, breaks = c(-Inf, alpha, Inf), labels = FALSE)
#' Y <- factor(Y,ordered=T)
#'
#' # Define Gaussian prior
#' prior <- list(mu0 = rep(0, 3), S0 = diag(3), Q0 = diag(3))
#'
#' # Fit the model using Expectation Propagation
#' fit <- viord(Y, X, prior = prior, algorithm = "EP")
#'
#' # Summarize results
#' summary(fit)
#'
#' # Predictive probabilities for new data
#' pred <- predict(fit, Xn = X)
#'
#' # Posterior simulation from approximate posterior
#' sims <- simulate(fit, nMC = 100)
#' }
#' @export
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
  out$prior = prior
  class(out) = 'viord'
  return(out)
}
