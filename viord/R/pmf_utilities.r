#' Posterior Sampling from a Fitted \code{viord} Model
#'
#' Generates posterior samples from the approximate posterior for the regression coefficients from a fitted
#' cumulative probit model obtained via \code{\link{viord}}.
#'
#' When the model is estimated using the \emph{Partially Factorized Mean-Field (PMF)}
#' algorithm, independent Monte Carlo sampling is used to
#' chatacterize the marginal approximate posterior of \eqn{\beta}. In this case, the user must supply the
#' original data (\code{Y}, \code{X}) and the prior specification used for model fitting.
#' The sampling relies on independent truncated normal draws from the optimal density \eqn{q^\star_{\tiny PMF}(z)}
#'
#' For models estimated using \emph{Expectation Propagation (EP)} or
#' \emph{Mean-Field (MF)}, the approximate posterior of \eqn{\beta}
#' is Gaussian with mean \code{object$est$m} and covariance \code{object$est$S};
#' hence, posterior samples are drawn directly from this multivariate normal
#' distribution.
#'
#' @param object A fitted object of class \code{"viord"}.
#' @param Y Optional. Ordinal response vector (factor or integer) required
#'   only if \code{object$algorithm == "PMF"}.
#' @param X Optional. Design matrix used in the model, required only if
#'   \code{object$algorithm == "PMF"}.
#' @param prior Optional. A list containing prior quantities
#'   (\code{mu0}, \code{S0}, \code{Q0}), required only if
#'   \code{object$algorithm == "PMF"}.
#' @param nMC Integer. Number of Monte Carlo samples to draw (default \code{1e3}).
#' @param ... Additional arguments (ignored).
#'
#' @return
#' A numeric matrix of posterior samples of the regression coefficients,
#' with one row per sample.
#'
#'
#' @seealso
#' \code{\link{viord}} for model fitting,
#' \code{\link{summary.viord}} for summaries, and
#' \code{\link{predict.viord}} for predictive probabilities.
#'
#' @references
#' Aliverti, E. (2025).  
#' *Approximate Bayesian Inference for Cumulative Probit Regression Models*.  
#' \url{https://arxiv.org/abs/2511.06967}
#'
#' @examples
#' \dontrun{
#' # Fit model with PMF approximation
#' fit <- viord(Y, X, prior = prior, algorithm = "PMF")
#'
#' # Posterior sampling via Monte Carlo
#' samples <- simulate(fit, Y = Y, X = X, prior = prior, nMC = 1000)
#'
#' # Fit model with EP (Gaussian posterior)
#' fit_ep <- viord(Y, X, prior = prior, algorithm = "EP")
#'
#' # Sampling directly from Gaussian approximation
#' samples_ep <- simulate(fit_ep, nMC = 1000)
#' }
#'
#' @export
simulate.viord = function(object, Y, X, prior, nMC = 1e3) {
	method = toupper(object$algorithm)

	if (method == "PMF") {
		# --- PMF: sample from truncated normals ---
		require(truncnorm)
		tresh = object$alpha
		xiZ = object$est$xiZ
		sigmaZ = object$est$sigmaZ

		li = tresh[Y]
		ui = tresh[as.numeric(Y) + 1]
		V = solve(crossprod(X) + prior$Q0)
		L = chol(V)
		XV = X %*% V
		Vprior = V %*% prior$Q0 %*% prior$mu

		ZMC = matrix(0, nrow(X), nMC)
		for (i in seq_len(nrow(X))) {
			ZMC[i, ] = rtruncnorm(
					       nMC,
					       a = li[i],
					       b = ui[i],
					       mean = xiZ[i],
					       sd = sigmaZ[i]
			)
		}

		b0 = matrix(rnorm(nMC * ncol(X)), ncol = nMC)
		m = apply(ZMC, 2, function(zs) Vprior + t(XV) %*% zs)
		beta_samp = m + L %*% b0

		rm(ZMC, b0)
		gc()
		out = t(beta_samp)

	} else if (method %in% c("MF", "EP")) {
		# --- MF and EP: Gaussian sampling from N(m, S) ---
		m = object$est$m
		S = object$est$S
		L = chol(S)
		p = length(m)

		beta_samp = matrix(0, nMC, p)
		for (i in seq_len(nMC)) {
			beta_samp[i, ] = m + L %*% rnorm(p)
		}

		out = beta_samp

	} else {
		stop("Unknown method: must be one of 'PMF', 'MF', or 'EP'.")
	}

	class(out) = c("simulate_viord", class(out))
	attr(out, "method") = method
	attr(out, "nMC") = nMC
	return(out)
}
