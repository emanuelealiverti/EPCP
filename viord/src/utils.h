#ifndef VIORD_UTILS_H
#define VIORD_UTILS_H

#include <RcppArmadillo.h>

// include truncated normal utilities
#include "truncnorm.h"

using namespace Rcpp;
using namespace arma;

// -----------------------------------------------------------------------------
// Truncated normal utilities
// -----------------------------------------------------------------------------

// Asymptotic expansion (Abramowitz & Stegun 26.2.13)
double zeta0(double& x);

// Expectation of truncated normal
double zeta1(double ca, double cb);

// Function combining first and second moments of truncated normal
double z1p2_z2(double ca, double cb);

// -----------------------------------------------------------------------------
// ELBO functions
// -----------------------------------------------------------------------------

// Mean-field ELBO
double elbo_mf(arma::vec mu_b, arma::mat Q0, arma::mat mu0, double lp);

// Partially mean-field ELBO
double elbo_pmf(arma::mat X, arma::mat V, arma::mat XV,
                arma::vec xiZ, arma::vec sigma2Z, arma::vec meanZ,
                arma::vec Xmu0, arma::vec mu0, double lp);

// -----------------------------------------------------------------------------
// Gaussian log normalizing constant
// -----------------------------------------------------------------------------

// Log of Gaussian normalizing constant (used for EP marginal likelihood)
double lPsi(const arma::vec& r, const arma::mat& Q);

#endif // VIORD_UTILS_H
