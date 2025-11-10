#ifndef VIORD_TRUNCNORM_H
#define VIORD_TRUNCNORM_H

//#include <R.h>
//#include <Rmath.h>
//#include <Rinternals.h>


#include <RcppArmadillo.h>

// -----------------------------------------------------------------------------
// Truncated normal distribution utilities
// -----------------------------------------------------------------------------

// Expected value (mean) of left-truncated normal
double e_lefttruncnorm(double a, double mean, double sd);

// Expected value (mean) of doubly truncated normal
double e_truncnorm(double a, double b, double mean, double sd);

// Expected value (mean) of right-truncated normal
double e_righttruncnorm(double b, double mean, double sd);

// Variance of left-truncated normal
double v_lefttruncnorm(double a, double mean, double sd);

// Variance of right-truncated normal
double v_righttruncnorm(double b, double mean, double sd);

// Variance of doubly truncated normal
double v_truncnorm(double a, double b, double mean, double sd);


#endif // VIORD_TRUNCNORM_H
