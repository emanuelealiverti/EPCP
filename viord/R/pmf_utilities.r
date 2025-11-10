#' Sample from a fitted PMF viord model
#'
#' Draws posterior samples of regression coefficients for a \code{viord}
#' object estimated with the PMF algorithm.
#'
#' @param object A fitted object of class \code{"viord"} obtained with method = "PMF".
#' @param Y A factor or numeric vector of ordinal responses.
#' @param X The design matrix used in the model.
#' @param prior A list containing prior quantities (\code{mu0}, \code{S0}, \code{Q0}).
#' @param nMC Number of Monte Carlo samples to draw (default 1e3).
#' @param ... Ignored, for compatibility.
#'
#' @return A matrix of posterior samples of regression coefficients,
#'         with one row per sample.
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
