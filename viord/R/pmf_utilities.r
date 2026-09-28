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
#' For models estimated using \emph{Expectation Propagation (EP)},
#' \emph{Mean-Field (MF)}, or \code{VB_prior}, the approximate posterior of \eqn{\beta}
#' is Gaussian with mean \code{object$est$m} and covariance \code{object$est$S};
#' hence, posterior samples are drawn directly from this multivariate normal
#' distribution.
#'
#' For \code{PMF_mixed} fits, the joint approximate posterior of
#' \eqn{(\beta, u)} is sampled; the returned matrix contains the fixed effects,
#' and the random-effect draws are stored in the \code{"ranef"} attribute.
#'
#' @param object A fitted object of class \code{"viord"}.
#' @param nsim Integer. Number of posterior samples requested by the
#'   \code{\link[stats]{simulate}} generic.
#' @param seed Optional random seed.
#' @param Y Optional. Ordinal response vector (factor or integer) required
#'   only for PMF and PMF_mixed fits.
#' @param X Optional. Design matrix used in the model, required only for
#'   PMF and PMF_mixed fits.
#' @param prior Optional. A list containing prior quantities
#'   (\code{mu0}, \code{S0}, \code{Q0}), required only if
#'   \code{object$algorithm == "PMF"}.
#' @param nMC Integer. Number of Monte Carlo samples to draw. Defaults to
#'   \code{nsim}.
#' @param Z Optional. Random-effects design matrix used in the model, required
#'   only if \code{object$algorithm == "PMF_mixed"}.
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
simulate.viord = function(object, nsim = 1, seed = NULL, Y = NULL, X = NULL,
			  prior = NULL, nMC = nsim, Z = NULL, ...) {
	if(!is.null(seed)){
		set.seed(seed)
	}
	method = toupper(object$algorithm)

	u_samp = NULL
	if (method %in% c("PMF", "PMF_MIXED")) {
		# --- PMF / PMF_mixed: sample from truncated normals ---
		tresh  = object$alpha
		xiZ    = object$est$xiZ
		sigmaZ = object$est$sigmaZ

		li = tresh[Y]
		ui = tresh[as.numeric(Y) + 1]

		p_fix = ncol(X)
		if (method == "PMF_MIXED") {
			if (is.null(Z))
				stop("For PMF_mixed simulation, you must specify Z.")
			X      = cbind(X, Z)
			jp     = pmf_mixed_prior(object)
			Q0_eff = jp$Q0
			mu0    = jp$mu0
		} else {
			Q0_eff = prior$Q0
			mu0    = prior$mu0
		}
		V      = solve(crossprod(X) + Q0_eff)
		L      = chol(V)
		XV     = X %*% V
		Vprior = V %*% (Q0_eff %*% mu0)

		ZMC = matrix(0, nrow(X), nMC)
		for (i in seq_len(nrow(X))) {
			ZMC[i, ] = truncnorm::rtruncnorm(
					       nMC,
					       a = li[i],
					       b = ui[i],
					       mean = xiZ[i],
					       sd = sigmaZ[i]
			)
		}

		b0_samp = matrix(rnorm(nMC * ncol(X)), ncol = nMC)
		m = apply(ZMC, 2, function(zs) Vprior + t(XV) %*% zs)
		beta_samp = m + L %*% b0_samp

		rm(ZMC, b0_samp)
		gc()
		out = t(beta_samp)
		if (method == "PMF_MIXED") {
			u_samp = out[, -seq_len(p_fix), drop = FALSE]
			colnames(u_samp) = object$u.names
			out = out[, seq_len(p_fix), drop = FALSE]
		}

	} else if (method %in% c("MF", "EP", "VB_PRIOR")) {
		# --- MF, VB_prior and EP: Gaussian sampling from N(m, S) ---
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
		stop("Unknown method: must be one of 'PMF', 'PMF_mixed', 'MF', 'VB_prior', or 'EP'.")
	}

	class(out) = c("simulate_viord", class(out))
	attr(out, "ranef") = u_samp
	attr(out, "method") = method
	attr(out, "nMC") = nMC
	return(out)
}


# Joint Gaussian prior on (beta, u) implied by a PMF_mixed fit, at the
# E[1/sigma2_g] used to build the returned q(z)
#' @noRd
pmf_mixed_prior = function(object) {
	groups = as.character(object$Z_group)
	g_id   = match(groups, unique(groups))
	tau_u  = object$est$sigma_u2_inv_mean[g_id]
	p = length(object$prior$mu0)
	q = length(tau_u)
	Q = matrix(0, p + q, p + q)
	Q[seq_len(p), seq_len(p)] = object$prior$Q0
	diag(Q)[p + seq_len(q)] = tau_u
	list(mu0 = c(object$prior$mu0, rep(0, q)), Q0 = Q)
}
