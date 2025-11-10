#' Posterior predictive probabilities for a fitted \code{viord} model
#'
#' @description
#' Computes posterior predictive probabilities or predicted classes for new
#' observations from a fitted variational ordinal regression model.
#'
#' Depending on the inference algorithm used to fit the model:
#' \itemize{
#'   \item For \strong{EP} and \strong{MF} methods, predictive probabilities are
#'   obtained using a Gaussian approximation of the posterior (analytical integration).
#'   \item For the \strong{PMF} method, predictive probabilities are computed
#'   using a Monte Carlo procedure based on truncated normal sampling.
#' }
#'
#' @param object A fitted model of class \code{"viord"}.
#' @param Xn A numeric matrix of new covariate values (one row per observation).
#' @param Y (PMF only) The ordinal response used in fitting the model.
#' @param X (PMF only) The design matrix used in fitting the model.
#' @param prior (PMF only) A list containing the prior specification
#'   (\code{mu0}, \code{S0}, and \code{Q0}).
#' @param nMC Number of Monte Carlo samples for the PMF predictive distribution
#'   (default \code{1e3}).
#' @param type A character string specifying the type of prediction:
#'   \code{"prob"} returns the full matrix of posterior predictive probabilities,
#'   while \code{"class"} returns only the most probable categories
#'   (default \code{"class"}).
#' @param ... Further arguments (currently ignored).
#'
#' @return
#' Depending on \code{type}:
#' \itemize{
#'   \item If \code{type = "prob"}: a numeric matrix of posterior predictive probabilities.
#'   \item If \code{type = "class"}: a numeric vector of predicted category indices.
#' }
#'
#' @examples
#' \dontrun{
#' fit_mf <- viord(Y, X, prior, method = "MF")
#' pred_class <- predict(fit_mf, Xn = X_new, type = "class")
#' pred_prob  <- predict(fit_mf, Xn = X_new, type = "prob")
#'
#' fit_pmf <- viord(Y, X, prior, method = "PMF")
#' pred_pmf <- predict(fit_pmf, Y = Y, X = X, prior = prior, Xn = X_new, nMC = 1000)
#' }
#'
#' @export
predict.viord <- function(object, Xn, Y = NULL, X = NULL, prior = NULL,
                          nMC = 1e3, type = c("class", "prob"), ...) {
  type <- match.arg(type)
  method <- toupper(object$algorithm)
  tresh <- object$alpha

  if (method == "PMF") {
    if (is.null(Y) || is.null(X) || is.null(prior))
      stop("For PMF prediction, you must specify Y, X, and prior.")

    xiZ <- object$est$xiZ
    sigmaZ <- object$est$sigmaZ

    pred <- pred_pmf(
      Y = Y,
      X = X,
      tresh = tresh,
      Xn = Xn,
      prior = prior,
      xiZ = xiZ,
      sigmaZ = sigmaZ,
      nMC = nMC
    )

  } else if (method %in% c("MF", "EP")) {
    m <- object$est$m
    S <- object$est$S
    pred <- pred_gauss(
      Xn = Xn,
      tresh = tresh,
      m = m,
      R = S
    )

  } else {
    stop("Unknown method: must be one of 'PMF', 'MF', or 'EP'.")
  }

  if (type == "class") {
    return(pred$y_hat)
  } else {
    return(pred$y_pr)
  }
}



## Compute predicitve probabilities for the ordinal probit regression
# from a gaussian approximation with variance R and mean m
# From EP or variational gaussian approx (compute the integral analytically)
pred_gauss = function(Xn, tresh, m, R) {
	loc = Xn %*% m
	sc = rowSums((Xn %*% R)  * Xn)
	pr_s = sapply(tresh, function(x) pnorm( (x-loc) / sqrt(1+sc)))
	y_pr = t(apply(pr_s, 1, diff))
	y_hat = apply(y_pr, 1, which.max)
	return(list(y_pr = y_pr, y_hat = y_hat))
}

# from PMF (monte carlo independent via truncated normals)
pred_pmf = function(Y, X, tresh, Xn, prior,  xiZ, sigmaZ, nMC = 1e3) {
	# simulate from independent truncated normals
	li = tresh[Y]
	ui = tresh[as.numeric(Y) + 1]
	V = solve(crossprod(X) + prior$Q0)
	Hn = Xn %*% V %*% t(X)
	Hnp = Xn %*% V %*% (prior$Q0 %*% prior$mu0)
	ZMC = matrix(0, NROW(X), nMC)
	for(i in 1:NROW(X)){
		ZMC[i,] = rtruncnorm(nMC,a = li[i],b =  ui[i], mean = xiZ[i], sd = sigmaZ[i])
	}
	sc = sqrt(1 + rowSums((Xn %*% V)  * Xn))

	# vectorized across n, single MC sample
	y_prMC = array(0, dim = c(NROW(Xn), length(tresh)-1, nMC)) 
	for(tt in 1:nMC){
		pr_s = sapply(tresh, function(x) pnorm( (x - Hn %*% ZMC[,tt] - Hnp) / sc))
		y_prMC[,,tt] = t(apply(pr_s,1,diff))
	}
	y_pr = apply(y_prMC,c(1,2), mean)
	y_hat = apply(y_pr, 1, which.max)
	# free up memory before closing function
	rm(y_prMC, ZMC,Hn); gc();
	return(list(y_pr = y_pr, y_hat = y_hat))
}
