## ---- Setup ------------------------------------------------------------------
rm(list = ls())
library(viord)
data(mcycle, package = "MASS")

x <- mcycle$times
y <- mcycle$accel

## Ordinal response: discretise acceleration into K = 5 ordered categories
K  <- 5
Yt <- cut(y, breaks = K, include.lowest = TRUE)

plot(x, as.numeric(Yt), xlab = "times (ms)", ylab = "category",
     main = "Ordinal response (mcycle)")

n <- length(x)

## ---- O'Sullivan penalized-spline basis ------------------------------------
## ZOSull (Wand 2015) constructs an O'Sullivan basis whose columns can be used
## either as fixed effects (with an IG prior on the common variance) or as
## random effects (the classic penalized-spline mixed-model parameterisation).

Xz    <- ZOSull(x)
xgrid <- seq(min(x), max(x), length.out = 300)
Xn    <- ZOSull(xgrid,
                range.x  = attr(Xz, "range.x"),
                intKnots = attr(Xz, "intKnots"))
p <- ncol(Xz)

z95 <- qnorm(0.975)

## ============================================================
## Approach 1 — ZOSull as fixed effects with IG prior on beta
## ============================================================
## The IG prior (a0, b0) on sigma_b2 acts as an automatic penalty on the
## spline coefficients; its scale is estimated from the data.

prior_ig <- list(mu0 = rep(0, p), a0 = 1, b0 = 2)

fit_vb <- viord(Y = Yt, X = Xz, prior = prior_ig, algorithm = "VB_prior")
summary(fit_vb)

## Posterior mean and pointwise 95% credible interval
f_vb  <- drop(Xn %*% coef(fit_vb))
se_vb <- sqrt(rowSums((Xn %*% vcov(fit_vb)) * Xn))

cat("VB_prior (fixed) — sigma_b2: mean =", round(fit_vb$est$sigma_b2_mean, 4),
    "  E[1/sigma2] =", round(fit_vb$est$sigma_b2_inv_mean, 4), "\n")

## ============================================================
## Approach 2 — ZOSull as random effects (mixed-model spline)
## ============================================================
## In the mixed-model parameterisation of a penalized spline (Wand 2003),
## the smooth deviations from the linear trend are treated as random effects
## u ~ N(0, sigma_u^2 I).  The estimated sigma_u^2 plays the role of the
## smoothing parameter.
##
## Fixed part: the linear trend x (centred); the intercept is absorbed by
## viord's thresholds.
## Random part: ZOSull basis columns, all in a single variance component.
##
## Two algorithms are compared:
##   * VB_prior  — mean-field VB, IG prior on sigma_b2 and sigma_u2;
##   * PMF_mixed — partially factorized VB, known Gaussian prior on the fixed
##                 effects and a half-Cauchy prior sigma_u ~ C+(0, s_sigma).

x_c     <- x - mean(x)                           # centre x for identifiability
X_re    <- matrix(x_c, n, 1, dimnames = list(NULL, "times_c"))
xgrid_c <- xgrid - mean(x)
Wn      <- cbind(xgrid_c, Xn)                    # joint design on the grid

prior_re  <- list(mu0 = 0, a0 = 1, b0 = 2, au0 = 1, bu0 = 1)
prior_hc  <- list(mu0 = 0, Q0 = matrix(1 / 100), s_sigma = 1)

fit_re <- viord(Y       = Yt,
                X       = X_re,
                Z       = Xz,
                Z_group = rep(0, ncol(Xz)),
                prior   = prior_re,
                algorithm = "VB_prior")

fit_hc <- viord(Y       = Yt,
                X       = X_re,
                Z       = Xz,
                Z_group = rep(0, ncol(Xz)),
                prior   = prior_hc,
                algorithm = "PMF_mixed")

summary(fit_re)
summary(fit_hc)

## Fitted smooth on the grid: linear fixed part + random-effect part, with
## pointwise 95% credible intervals from the joint posterior of (beta, u)
smooth_fit <- function(fit) {
  f  <- drop(Wn %*% fit$est$m_joint)
  se <- sqrt(rowSums((Wn %*% fit$est$S_joint) * Wn))
  list(f = f, se = se)
}
sm_re <- smooth_fit(fit_re)
sm_hc <- smooth_fit(fit_hc)

cat("VB_prior (random)  — sigma_u2 mean =", round(fit_re$est$sigma_u2_mean, 4), "\n")
cat("PMF_mixed (random) — sigma_u2 mean =", round(fit_hc$est$sigma_u2_mean, 4), "\n")

## ---- Credible bands --------------------------------------------------------
plot_band <- function(f, se, col, main) {
  plot(xgrid, f, type = "n", ylim = ylim,
       xlab = "times (ms)", ylab = "linear predictor", main = main)
  polygon(c(xgrid, rev(xgrid)), c(f + z95 * se, rev(f - z95 * se)),
          col = adjustcolor(col, 0.20), border = NA)
  lines(xgrid, f, col = col, lwd = 2)
  rug(x, col = col)
}

ylim <- range(c(f_vb - z95 * se_vb, f_vb + z95 * se_vb,
                sm_re$f - z95 * sm_re$se, sm_re$f + z95 * sm_re$se,
                sm_hc$f - z95 * sm_hc$se, sm_hc$f + z95 * sm_hc$se))

par(mfrow = c(1, 3))
plot_band(f_vb,    se_vb,    "steelblue", "VB_prior — fixed effects")
plot_band(sm_re$f, sm_re$se, "darkgreen", "VB_prior — random effects")
plot_band(sm_hc$f, sm_hc$se, "tomato3",   "PMF_mixed — random effects")
par(mfrow = c(1, 1))

## ---- Overlay all three estimates (centred for comparability) ---------------
f_list <- list("VB_prior fixed"   = f_vb    - mean(f_vb),
               "VB_prior random"  = sm_re$f - mean(sm_re$f),
               "PMF_mixed random" = sm_hc$f - mean(sm_hc$f))

cols <- c("steelblue", "darkgreen", "tomato3")
ltys <- c(1, 3, 2)

ylim2 <- range(unlist(f_list))
plot(xgrid, f_list[[1]], type = "l", col = cols[1], lwd = 2,
     ylim = ylim2, xlab = "times (ms)", ylab = "centred linear predictor",
     main = "O'Sullivan penalized splines: algorithm comparison")
for (k in 2:3)
    lines(xgrid, f_list[[k]], col = cols[k], lwd = 2, lty = ltys[k])
rug(x)
legend("topright", legend = names(f_list),
       col = cols, lwd = 2, lty = ltys, bty = "n")
