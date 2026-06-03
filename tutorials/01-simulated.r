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

## ============================================================
## Approach 1 — ZOSull as fixed effects with IG prior on beta
## ============================================================
## Both VB_prior and PMF_prior share the same prior specification.
## The IG prior (a0, b0) on sigma_b2 acts as an automatic penalty on the
## spline coefficients; its scale is estimated from the data.

prior_ig <- list(mu0 = rep(0, p), a0 = 1, b0 = 2)

fit_vb  <- viord(Y = Yt, X = Xz, prior = prior_ig, algorithm = "VB_prior")
fit_pmf <- viord(Y = Yt, X = Xz, prior = prior_ig, algorithm = "PMF_prior")

summary(fit_vb)
summary(fit_pmf)

## Posterior means and pointwise 95% credible intervals
f_vb  <- drop(Xn %*% coef(fit_vb))
f_pmf <- drop(Xn %*% coef(fit_pmf))

se_vb  <- sqrt(rowSums((Xn %*% vcov(fit_vb))  * Xn))
se_pmf <- sqrt(rowSums((Xn %*% vcov(fit_pmf)) * Xn))

z95  <- qnorm(0.975)
ylim <- range(c(f_vb  - z95 * se_vb,  f_vb  + z95 * se_vb,
                f_pmf - z95 * se_pmf, f_pmf + z95 * se_pmf))

par(mfrow = c(1, 2))

# VB_prior
plot(xgrid, f_vb, type = "n", ylim = ylim,
     xlab = "times (ms)", ylab = "linear predictor",
     main = "VB_prior — O'Sullivan fixed effects")
polygon(c(xgrid, rev(xgrid)),
        c(f_vb + z95 * se_vb, rev(f_vb - z95 * se_vb)),
        col = adjustcolor("steelblue", 0.20), border = NA)
lines(xgrid, f_vb,  col = "steelblue", lwd = 2)
rug(x, col = "dodgerblue")

# PMF_prior
plot(xgrid, f_pmf, type = "n", ylim = ylim,
     xlab = "times (ms)", ylab = "linear predictor",
     main = "PMF_prior — O'Sullivan fixed effects")
polygon(c(xgrid, rev(xgrid)),
        c(f_pmf + z95 * se_pmf, rev(f_pmf - z95 * se_pmf)),
        col = adjustcolor("tomato3", 0.20), border = NA)
lines(xgrid, f_pmf, col = "tomato3", lwd = 2)
rug(x, col = "tomato")

par(mfrow = c(1, 1))

cat("VB_prior  — sigma_b2: mean =", round(fit_vb$est$sigma_b2_mean,  4),
    "  E[1/sigma2] =", round(fit_vb$est$sigma_b2_inv_mean,  4), "\n")
cat("PMF_prior — sigma_b2: mean =", round(fit_pmf$est$sigma_b2_mean, 4),
    "  E[1/sigma2] =", round(fit_pmf$est$sigma_b2_inv_mean, 4), "\n")

## ============================================================
## Approach 2 — ZOSull as random effects (mixed-model spline)
## ============================================================
## In the mixed-model parameterisation of a penalized spline (Wand 2003),
## the smooth deviations from the linear trend are treated as random effects
## u ~ N(0, sigma_u^2 I).  The estimated sigma_u^2 plays the role of the
## smoothing parameter.  Only VB_prior supports random effects.
##
## Fixed part: the linear trend x (centred), absorbed by viord's thresholds.
## Random part: ZOSull basis columns, all in a single variance component.

x_c  <- x - mean(x)                              # centre x for identifiability
X_re <- matrix(x_c, n, 1)

prior_re <- list(mu0 = 0, a0 = 1, b0 = 2, au0 = 1, bu0 = 1)

fit_re <- viord(Y       = Yt,
                X       = X_re,
                Z       = Xz,
                Z_group = rep(0, ncol(Xz)),
                prior   = prior_re,
                algorithm = "VB_prior")

summary(fit_re)

## Reconstruct fitted smooth on the grid
xgrid_c <- xgrid - mean(x)
f_re_lin <- xgrid_c * coef(fit_re)              # linear fixed-effect part
f_re_smo <- drop(Xn %*% ranef(fit_re)[["0"]])   # smooth random-effect part
f_re     <- f_re_lin + f_re_smo

cat("VB_prior (random) — sigma_u2 mean =",
    round(fit_re$est$sigma_u2_mean, 4),
    "  sigma_b2 mean =", round(fit_re$est$sigma_b2_mean, 4), "\n")

## ---- Overlay all three estimates (centred for comparability) ---------------
f_list <- list("VB_prior fixed"  = f_vb  - mean(f_vb),
               "PMF_prior fixed" = f_pmf - mean(f_pmf),
               "VB_prior random" = f_re  - mean(f_re))

cols <- c("steelblue", "tomato3", "darkgreen")
ltys <- c(1, 2, 3)

ylim2 <- range(unlist(f_list))
plot(xgrid, f_list[[1]], type = "l", col = cols[1], lwd = 2,
     ylim = ylim2, xlab = "times (ms)", ylab = "centred linear predictor",
     main = "O'Sullivan penalized splines: algorithm comparison")
for (k in 2:3)
    lines(xgrid, f_list[[k]], col = cols[k], lwd = 2, lty = ltys[k])
rug(x)
legend("topright", legend = names(f_list),
       col = cols, lwd = 2, lty = ltys, bty = "n")
