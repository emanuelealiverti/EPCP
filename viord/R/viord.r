#' Approximate Bayesian Inference for Cumulative Probit Models
#'
#' Provides a unified interface for scalable approximate Bayesian inference under
#' the cumulative probit model using one of five algorithms:
#' \emph{Expectation Propagation (EP)}, \emph{Mean-Field Variational Bayes (MF)},
#' \emph{Partially Factorized Mean-Field (PMF)}, a mean-field VB method with
#' an inverse-gamma prior on the Gaussian prior variance (\code{VB_prior}), or a
#' PMF mixed model with half-Cauchy priors on the random-effect standard
#' deviations (\code{PMF_mixed}).
#'
#' Threshold (cutoff) parameters are estimated via approximate marginal likelihood,
#' alternating between the optimization of the thresholds (via Newton–Raphson steps)
#' and the internal optimization based on EP, MF, or PMF.  
#' The core computational routines are implemented in C++ for efficiency,
#' while this function provides a high-level R interface.
#'
#' @param Y Ordinal response variable (an ordered factor).
#' @param X Design matrix of covariates.
#' @param Z Optional random-effects design matrix. Random-effect structures can
#'   be column-bound in R before calling \code{viord}; columns belonging to the
#'   same variance component are identified by \code{Z_group}.
#' @param Z_group Optional vector of length \code{ncol(Z)} indicating the
#'   variance-component group for each column of \code{Z}. Currently used only
#'   with \code{algorithm = "VB_prior"} and \code{algorithm = "PMF_mixed"}.
#' @param prior A list containing prior parameters. For \code{"EP"}, \code{"MF"},
#'   and \code{"PMF"}, provide \code{mu0} (prior mean), \code{S0} (prior
#'   covariance), and \code{Q0} (prior precision matrix). For \code{"VB_prior"},
#'   provide \code{mu0}, \code{a0}, and \code{b0},
#'   corresponding to \eqn{\beta \mid \sigma_b^2 \sim N(\mu_0, \sigma_b^2 I_p)}
#'   and \eqn{\sigma_b^2 \sim IG(a_0, b_0)}. For \code{"VB_prior"} with random
#'   effects, also provide \code{au0} and \code{bu0}. For \code{"PMF_mixed"},
#'   provide \code{mu0} and \code{Q0} (known prior mean and precision of the
#'   fixed effects) and \code{s_sigma}, the scale of the half-Cauchy prior
#'   shared by all random-effect standard deviations.
#' @param algorithm Character string specifying the inference algorithm to use:
#'   one of \code{"EP"}, \code{"MF"}, \code{"PMF"}, \code{"VB_prior"}, or
#'   \code{"PMF_mixed"}.
#' @param control A list of fitting options, as returned by
#'   \code{\link{viord.control}}. It collects the maximum number of iterations
#'   and the convergence tolerance of the inner algorithm and of the threshold
#'   loop, the convergence criterion, the warm start, starting or fixed
#'   thresholds, and how much information is printed.
#' @return A list containing the estimated model quantities and convergence information.  
#'   The output includes:
#'   \itemize{
#'     \item \code{est}: posterior (approximate) mean vector \eqn{m} and covariance matrix \eqn{S};
#'     \item \code{alpha}: estimated thresholds (cutoffs);
#'     \item \code{algorithm}: the selected inference algorithm;
#'     \item \code{prior}: the prior provided as input
#'     \item \code{it_outer} and \code{conv_outer}: number of threshold updates
#'       performed and whether the threshold loop converged;
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
#' A Gaussian prior \eqn{\beta \sim N_p(\mu_0, S_0)} is assumed for the regression coefficients
#' under \code{"EP"}, \code{"MF"}, and \code{"PMF"}. Under \code{"VB_prior"},
#' the fixed-effect prior is
#' \eqn{\beta \mid \sigma_b^2 \sim N_p(\mu_0, \sigma_b^2 I_p)}
#' with \eqn{\sigma_b^2 \sim IG(a_0, b_0)}. If \code{Z} is supplied, the
#' random-effect coefficients have independent group-specific priors
#' \eqn{u_j \mid \sigma_{u,g(j)}^2 \sim N(0, \sigma_{u,g(j)}^2)}.
#' Under \code{"PMF_mixed"}, \eqn{\beta \sim N_p(\mu_0, Q_0^{-1})} with known
#' \eqn{Q_0}, \eqn{u_j \mid \sigma_{u,g(j)}^2 \sim N(0, \sigma_{u,g(j)}^2)} and
#' \eqn{\sigma_{u,g} \sim C^+(0, s_\sigma)}, represented as
#' \eqn{\sigma_{u,g}^2 \mid a_g \sim IG(1/2, 1/a_g)},
#' \eqn{a_g \sim IG(1/2, 1/s_\sigma^2)}.
#' The thresholds \eqn{\alpha_1 < \dots < \alpha_{K-1}} are treated as nuisance parameters
#' and estimated by maximizing the (approximate) marginal likelihood.
#'
#' The \code{algorithm} argument controls the internal approximation used
#' for inference conditional on the current thresholds:
#' \itemize{
#'   \item \code{"EP"} – Expectation Propagation;
#'   \item \code{"MF"} – Mean-Field Variational Bayes;
#'   \item \code{"PMF"} – Partially Factorized Mean-Field Variational Bayes;
#'   \item \code{"VB_prior"} – Mean-Field Variational Bayes with an
#'     inverse-gamma update for the prior variances;
#'   \item \code{"PMF_mixed"} – Partially Factorized Mean-Field for mixed
#'     models, with half-Cauchy priors on the random-effect standard deviations.
#' }
#'
#' @references
#' Aliverti, E. (2025).
#' *Approximate Bayesian Inference for Cumulative Probit Regression Models*.
#' \url{https://arxiv.org/abs/2511.06967}
#' @seealso
#' \code{\link{viord.control}} for the fitting options,
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
                 Z = NULL, Z_group = NULL,
                 algorithm = c("EP", "MF", "PMF", "VB_prior", "PMF_mixed"),
                 control = viord.control()) {

  algorithm <- match.arg(algorithm)
  if (!(algorithm %in% c("VB_prior", "PMF_mixed")) && !is.null(Z) && NCOL(Z) > 1) {
    stop("Z and Z_group are currently supported only with algorithm = 'VB_prior' or 'PMF_mixed'.")
  }

  if (!is.list(control) || is.null(control$maxit_inner))
    stop("control must be a list as returned by viord.control().")

  if (!is.null(control$alpha_init) &&
      length(control$alpha_init) != nlevels(factor(Y)) - 1)
    stop("control$alpha_init must have length equal to the number of categories minus one.")

  if (algorithm == "EP") {
    out <- optim_ep_ml(Y = Y, X = X, prior = prior, control = control)

  } else {
    out <- optim_vb_ml(Y = Y, X = X, prior = prior,
                       vb_factor = algorithm,
                       Z = Z, Z_group = Z_group,
                       control = control)
  }

  out$coef.names = colnames(X)
  if (!is.null(Z) && NCOL(Z) > 1) {
    out$u.names = colnames(Z)
    if (is.null(out$u.names))
      out$u.names = paste0("u[", seq_len(NCOL(Z)), "]")
    out$Z_group = Z_group
  }
  out$algorithm = algorithm
  out$prior = prior
  out$control = control
  class(out) = 'viord'
  return(out)
}
