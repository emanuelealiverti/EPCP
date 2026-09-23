viord: Climate Policy Support Case Study
================
Emanuele Aliverti

This document illustrates the mixed-effects extension of the `viord`
package, using the TISP dataset (Mede and Cologna 2023) — an
international survey on public attitudes towards climate change
collected across 66 countries. We model ordinal support for climate
policies (`CLIM_POLSUPPORT`) as a function of demographic predictors,
with nonlinear effects of age and log-income represented via O’Sullivan
penalized splines (Wand and Ormerod 2008) and country-level random
intercepts capturing cross-national heterogeneity. Both the spline
coefficients and the country intercepts enter the model as random
effects, each block with its own variance. Inference uses the
`"PMF_mixed"` algorithm, a partially factorized mean-field approximation
that places a Gaussian prior with known variance on the fixed effects
and half-Cauchy priors (Wand et al. 2011) on the random-effect standard
deviations, estimated jointly from the data.

------------------------------------------------------------------------

First, we load the required packages.

``` r
library(viord)
```

## Data Loading

We download the data directly from the OSF repository of Mede and
Cologna (2023). The dataset is stored in a semicolon-separated CSV file
with European decimal notation, so we use `read.csv2` after downloading.
The parsed data are cached locally in `tisp.rds`, so that re-knitting
this document does not download the file again.

``` r
url_data   <- "https://osf.io/download/xjc4p"
cache_file <- "tisp.rds"   # local cache, not tracked by git

if (file.exists(cache_file)) {
  dd <- readRDS(cache_file)
} else {
  tmp_file <- tempfile(fileext = ".csv")
  curl::curl_download(url_data, tmp_file)
  dd <- read.csv2(tmp_file)
  unlink(tmp_file)
  saveRDS(dd, cache_file)
}
```

    ## Dataset dimensions: 69534 x 141

## Response Variable

Five items (`CLIM_POLSUPPORT_*`) measure support for specific climate
policies (fuel taxes, public transport, sustainable energy,
environmental protection, food taxes) on a 1–3 ordinal scale. Value 4
(“Not applicable”) is treated as missing following the codebook. We sum
the five items and discretise into `K = 5` ordered categories.

``` r
risp <- grep("CLIM_POLSUPPORT", names(dd), value = TRUE)
for (v in risp) dd[[v]][dd[[v]] == 4] <- NA
dd$y_sum <- rowSums(dd[, risp])
K <- 5
Yt <- cut(dd$y_sum, breaks = K, include.lowest = TRUE)
table(Yt, useNA = "always")
```

    ## Yt
    ## [4.99,7]    (7,9]   (9,11]  (11,13]  (13,15]     <NA> 
    ##     1476     4308    15574    22178    13017    12981

## Covariate Preparation

We retain demographic predictors: age, log-income, gender, education,
political orientation, religiosity, and urban/rural residence. After
listwise deletion of incomplete cases the analysis sample has 41,976
individuals from 66 countries.

``` r
vars_keep <- c("DEM_AGE", "DEM_INCOME_USD_log",
               "COUNTRY_CODE",
               "DEM_GENDER_male",
               "DEM_EDU",
               "DEM_POL_right",
               "DEM_RELIGIOUS",
               "DEM_RESIDENCE")
df <- dd[, vars_keep]
df$Yt <- Yt
df <- na.omit(df)
cat("Sample size after listwise deletion:", nrow(df), "\n")
```

    ## Sample size after listwise deletion: 41976

## Nonlinear Effects via Penalized Splines

We use the mixed-model representation of penalized splines (Wand and
Ormerod 2008). Each smooth is decomposed as
$f(x) = \beta x + \sum_k z_k(x) u_k$: the null space of the penalty (the
linear term) enters as a fixed effect, while the coefficients $u_k$ of
the penalized part are random effects $u_k \sim N(0, \sigma^2_f)$. The
variance $\sigma^2_f$ acts as the smoothing parameter and is estimated
from the data. The intercept is absorbed by the ordinal thresholds.

The bases are built with `mgcv`: `smoothCon()` constructs a thin-plate
regression spline with the sum-to-zero identifiability constraint
already absorbed, and `smooth2random()` returns the fixed/random
decomposition above. With `k = 15` each smooth contributes one fixed
column and 13 random columns, considerably fewer than an unreduced
basis.

Age and log-income are standardized first, so that the half-Cauchy scale
used below is meaningful for both smooths.

``` r
library(mgcv)

age    <- df$DEM_AGE
income <- df$DEM_INCOME_USD_log

std <- function(x, ref = x) (x - mean(ref)) / sd(ref)
age_s    <- std(age)
income_s <- std(income)

# penalized basis in mixed-model form: Xf spans the penalty null space
# (fixed effect), Xr the penalized part (random effects)
make_smooth <- function(x, k = 15) {
  sm <- smoothCon(s(x, bs = "tp", k = k), data = data.frame(x = x),
                  absorb.cons = TRUE)[[1]]
  re <- smooth2random(sm, "", type = 2)
  list(sm = sm, trans = re$trans.U %*% diag(re$trans.D),
       n_rand = ncol(re$rand$Xr))
}

# evaluate a smooth at new points, with the same transformation
eval_smooth <- function(S, x) {
  M <- PredictMat(S$sm, data.frame(x = x)) %*% S$trans
  list(Xr = M[, seq_len(S$n_rand), drop = FALSE],
       Xf = M[, -seq_len(S$n_rand), drop = FALSE])
}

S_age <- make_smooth(age_s)
S_inc <- make_smooth(income_s)

B_age <- eval_smooth(S_age, age_s)
B_inc <- eval_smooth(S_inc, income_s)
```

Prediction grids and the corresponding bases, transformed in the same
way, are prepared for plotting.

``` r
age_grid    <- seq(min(age),    max(age),    length.out = 200)
income_grid <- seq(min(income), max(income), length.out = 200)

Bn_age <- eval_smooth(S_age, std(age_grid,    age))
Bn_inc <- eval_smooth(S_inc, std(income_grid, income))
```

## Fixed-Effects Design Matrix

The fixed effects are the linear parts of the two smooths and the
remaining predictors (gender, education, political orientation,
religiosity, residence). No intercept is included; the ordinal
thresholds play its role.

``` r
X_other <- model.matrix(~ DEM_GENDER_male + DEM_EDU + DEM_POL_right +
                           DEM_RELIGIOUS + DEM_RESIDENCE,
                         data = df)[, -1]

X <- cbind(age_lin = B_age$Xf, inc_lin = B_inc$Xf, X_other)
colnames(X)[1:2] <- c("age_lin", "inc_lin")
```

    ## Fixed-effects design matrix dim: 41976 7

## Random Effects

Country-level random intercepts are encoded with one indicator column
per country, as `mgcv` does for `s(country, bs = "re")`: no sum-to-zero
constraint is imposed, and the zero-mean Gaussian prior with estimated
variance is what identifies the 66 effects against the ordinal
thresholds. The random-effects design stacks the two penalized blocks
and the country indicators, and `Z_group` assigns each column to its own
variance component: one for the age smooth, one for the log-income
smooth, and one for the countries.

``` r
country   <- factor(df$COUNTRY_CODE)
Z_country <- model.matrix(~ country - 1)
colnames(Z_country) <- levels(country)

Z       <- cbind(B_age$Xr, B_inc$Xr, Z_country)
Z_group <- c(rep("age",     ncol(B_age$Xr)),
             rep("income",  ncol(B_inc$Xr)),
             rep("country", ncol(Z_country)))
colnames(Z) <- c(paste0("age_z", seq_len(ncol(B_age$Xr))),
                 paste0("inc_z", seq_len(ncol(B_inc$Xr))),
                 levels(country))
```

    ## Random-effects design matrix dim: 41976 92

    ## Z_group
    ##     age country  income 
    ##      13      66      13

## Model Fitting

The fixed effects receive a Gaussian prior with known covariance
$10\, I$ (`Q0` is its inverse). The standard deviation of each of the
three random-effect blocks receives a half-Cauchy prior with scale
`s_sigma = 1`, represented as $\sigma^2_g \mid a_g \sim IG(1/2, 1/a_g)$,
$a_g \sim IG(1/2, 1/s_\sigma^2)$. We fit the model with
`algorithm = "PMF_mixed"`. Fitting options are collected by
`viord.control()`: the tolerances are relative to the ELBO, so `1e-7` is
a tighter criterion than the default here and stabilises the smaller
variance components, and the iteration caps are raised accordingly so
that both loops converge rather than stopping at the cap.

``` r
p     <- ncol(X)
prior <- list(mu0 = rep(0, p), Q0 = diag(1 / 10, p), s_sigma = 1)

fit <- viord(Y       = df$Yt,
             X       = X,
             Z       = Z,
             Z_group = Z_group,
             prior   = prior,
             algorithm = "PMF_mixed",
             control   = viord.control(tol_inner   = 1e-7,
                                       tol_outer   = 1e-7,
                                       maxit_inner = 200,
                                       maxit_outer = 500))

summary(fit)
```

    ## 
    ## Summary of VI Ordinal Model
    ## Inference algorithm: PMF_mixed 
    ## 
    ## Posterior estimates:
    ##                 Estimate Std. Error
    ## age_lin         -0.0224   0.0350   
    ## inc_lin          0.0283   0.1651   
    ## DEM_GENDER_male  0.0048   0.0105   
    ## DEM_EDU          0.1965   0.0076   
    ## DEM_POL_right   -0.1784   0.0050   
    ## DEM_RELIGIOUS    0.0218   0.0043   
    ## DEM_RESIDENCE    0.1637   0.0124   
    ## 
    ## Threshold parameters (cutpoints):
    ##          Estimate
    ## alpha[1] -1.8118 
    ## alpha[2] -1.1130 
    ## alpha[3] -0.1270 
    ## alpha[4]  0.9822 
    ## 
    ## Random-effect variance posterior:
    ##                   a       b       Mean    E[1/sigma2]
    ## sigma_u2[age]      7.0000  0.0897  0.0149 78.0459    
    ## sigma_u2[income]   7.0000 33.5119  5.5853  0.2089    
    ## sigma_u2[country] 33.5000  2.2710  0.0699 14.7513    
    ## 
    ## Converged in 2 iterations. Approx. log marginal likelihood: -54419

The `summary` output reports the posterior mean and standard deviation
of the fixed effects, the estimated ordinal thresholds, and the
approximate posterior of the three variance components. Education,
political orientation (right-leaning), and urban residence have the
largest and most precisely estimated linear effects. The estimated
country-level variance (posterior mean ≈ 0.07) is modest but
non-negligible.

## Smooth Effects of Age and Log-Income

Each smooth combines its linear fixed effect with its spline random
effects. We recover the posterior mean and pointwise 95% credible
intervals from the joint approximate posterior of all coefficients
(`fit$est$m_joint`, `fit$est$S_joint`). Each curve is centred at its
average over the observed data, since the level of the smooth is not
separately identified from the thresholds.

``` r
m_all <- fit$est$m_joint
S_all <- fit$est$S_joint

# Posterior mean and pointwise sd of a centred smooth, given the positions
# of its coefficients in the joint vector
smooth_post <- function(W_grid, W_obs, idx) {
  C  <- sweep(W_grid, 2, colMeans(W_obs))
  f  <- drop(C %*% m_all[idx])
  se <- sqrt(rowSums((C %*% S_all[idx, idx]) * C))
  list(f = f, se = se)
}

idx_age <- c(which(colnames(X) == "age_lin"),
             p + which(Z_group == "age"))
idx_inc <- c(which(colnames(X) == "inc_lin"),
             p + which(Z_group == "income"))

sm_age <- smooth_post(cbind(Bn_age$Xf, Bn_age$Xr),
                      cbind(B_age$Xf, B_age$Xr), idx_age)
sm_inc <- smooth_post(cbind(Bn_inc$Xf, Bn_inc$Xr),
                      cbind(B_inc$Xf, B_inc$Xr), idx_inc)

z95 <- qnorm(0.975)

plot_smooth <- function(grid, sm, xlab, main, rug_x) {
  lwr <- sm$f - z95 * sm$se
  upr <- sm$f + z95 * sm$se
  plot(grid, sm$f, type = "n", ylim = range(lwr, upr),
       xlab = xlab, ylab = "Linear predictor", main = main)
  polygon(c(grid, rev(grid)), c(upr, rev(lwr)),
          col = adjustcolor("steelblue", alpha.f = 0.25), border = NA)
  lines(grid, sm$f, col = "steelblue", lwd = 2)
  lines(grid, lwr,  col = "steelblue", lwd = 1, lty = 2)
  lines(grid, upr,  col = "steelblue", lwd = 1, lty = 2)
  rug(rug_x, col = "dodgerblue")
}

par(mfrow = c(1, 2))
plot_smooth(age_grid,    sm_age, xlab = "Age",
            main = "Smooth effect of age",        rug_x = age)
plot_smooth(income_grid, sm_inc, xlab = "Log-income (USD)",
            main = "Smooth effect of log-income", rug_x = income)
par(mfrow = c(1, 1))
```

<img src="smooth_effects.png" alt="" width="2400" />

The effect of age is close to linear: support decreases mildly with age,
with wide credible intervals above 80 where data are sparse. The
estimated variance of that block is small (0.015), so the smooth is
heavily penalized towards the linear fit. The effect of log-income is
markedly non-monotone, with several local features and a rise in the
upper part of the distribution; its variance component is much larger
(5.59), so the data support a genuinely wiggly shape. Fitting the same
model with `mgcv::gam` by REML gives the same picture, with 2.1
effective degrees of freedom for age against 10.5 for log-income.

## Country Random Effects

We extract the posterior means of the 66 country effects with `ranef()`
and map them onto a world map.

``` r
library(ggplot2)
library(maps)

u_all <- ranef(fit)$country

country_effects <- data.frame(COUNTRY_CODE = names(u_all),
                              effect       = as.numeric(u_all))

world     <- subset(map_data("world"), region != "Antarctica")
world$iso <- iso.alpha(world$region, n = 3)
world     <- merge(world, country_effects,
                   by.x = "iso", by.y = "COUNTRY_CODE", all.x = TRUE)
world     <- world[order(world$order), ]

ggplot(world, aes(long, lat, group = group, fill = effect)) +
  geom_polygon(color = "gray80", linewidth = 0.15) +
  coord_quickmap() +
  scale_fill_gradient2(low = "#2166ac", mid = "lightyellow", high = "#d6604d",
                       midpoint = 0, na.value = "gray92",
                       name = "Country\neffect") +
  theme_minimal() +
  theme(panel.grid    = element_blank(),
        axis.text     = element_blank(),
        axis.ticks    = element_blank(),
        axis.title    = element_blank(),
        panel.background = element_rect(fill = "#ddeeff", color = NA)) +
  labs(title = "Country random effects on climate policy support")
```

<img src="country_map.png" alt="" width="3200" />

Countries shown in red have higher-than-average climate policy support
after conditioning on individual-level demographics; blue countries have
lower-than-average support. The range of country effects is
approximately (-0.71, 0.57), indicating meaningful cross-national
variation that would be ignored by a fixed-effects-only model.

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-main" class="csl-entry">

Aliverti, Emanuele. 2025. “Approximate Bayesian Inference for Cumulative
Probit Regression Models.” *arXiv Preprint arXiv:2511.06967*.

</div>

<div id="ref-tisp" class="csl-entry">

Mede, Niels G., and Viktoria Cologna. 2023. *The TISP Dataset*. OSF.
<https://doi.org/10.17605/OSF.IO/5C3QD>.

</div>

<div id="ref-wand:2008" class="csl-entry">

Wand, M. P., and J. T. Ormerod. 2008. “On Semiparametric Regression with
O’Sullivan Penalized Splines.” *Australian & New Zealand Journal of
Statistics* 50 (2): 179–98.

</div>

<div id="ref-wand:2011" class="csl-entry">

Wand, M. P., J. T. Ormerod, S. A. Padoan, and R. Frühwirth. 2011. “Mean
Field Variational Bayes for Elaborate Distributions.” *Bayesian
Analysis* 6 (4): 847–900.

</div>

</div>
