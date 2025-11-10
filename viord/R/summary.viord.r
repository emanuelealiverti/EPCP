#' Summary method for viord objects
#'
#' Summarizes the fitted variational ordinal model, showing posterior estimates,
#' uncertainty, and approximate log marginal likelihood.
#'
#' @param object An object of class \code{"viord"}.
#' @param level Confidence level for credible intervals (default 0.95).
#' @param ... Additional arguments (ignored).
#'
#' @return An object of class \code{"summary.viord"}.
#' @export
summary.viord <- function(object, ci = FALSE, level = 0.95, ...) {
  if (!inherits(object, "viord"))
    stop("object must be of class 'viord'")

  est <- object$est
  m <- est$m
  S <- est$S

  se <- sqrt(diag(S))
  z <- qnorm(1 - (1 - level) / 2)
  lower <- m - z * se
  upper <- m + z * se

  coef_names <- object$coef.names
  if (is.null(coef_names))
    coef_names <- paste0("beta[", seq_along(m), "]")

  coef_table <- cbind(m, se, lower, upper)
  rownames(coef_table) <- coef_names
  colnames(coef_table) <- c("Estimate", "Std. Error", "Lower", "Upper")
  if(!ci) coef_table = coef_table[,1:2] # drop intervals


  # Thresholds (cutpoints)
  alpha <- object$alpha
  thresholds <- NULL
  if (!is.null(alpha)) {
    alpha <- alpha[is.finite(alpha)]
    thresholds <- matrix(alpha, ncol = 1)
    rownames(thresholds) <- paste0("alpha[", seq_along(alpha), "]")
    colnames(thresholds) <- "Estimate"
  }

  # Determine algorithm and log marginal likelihood
  if (!is.null(est$logZ)) {
    log_marginal_lik <- est$logZ
    algorithm <- "EP"
  } else if (!is.null(est$elbo)) {
    log_marginal_lik <- est$elbo
    if (!is.null(est$sigmaZ)) {
      algorithm <- "PMF"
    } else {
      algorithm <- "MF"
    }
  } else {
    log_marginal_lik <- NA_real_
    algorithm <- "Unknown"
  }

  out <- list(
    coefficients = coef_table,
    thresholds = thresholds,
    log_marginal_lik = log_marginal_lik,
    algorithm = algorithm,
    n_it = est$it
  )
  class(out) <- "summary.viord"
  out
}
#' @export
print.summary.viord <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nSummary of VI Ordinal Model\n")
  cat("Inference algorithm:", x$algorithm, "\n\n")

  cat("Posterior estimates:\n")
  print(noquote(format(round(x$coefficients, digits = digits), nsmall = 3, justify = "right")))
  
  if (!is.null(x$thresholds)) {
    cat("\nThreshold parameters (cutpoints):\n")
    print(noquote(format(round(x$thresholds, digits = digits), nsmall = 3, justify = "right")))
  }

  cat(sprintf("\nConverged in %s iterations. ",x$n_it))
  cat("Approx. log marginal likelihood:",
      format(x$log_marginal_lik, digits = digits), "\n")
  invisible(x)
}
