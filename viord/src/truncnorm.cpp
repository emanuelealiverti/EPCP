/*
 * truncnorm.c - Implementation of truncated normal distribution
 *
 * Authors:
 *  Heike Trautmann  <trautmann@statistik.uni-dortmund.de>
 *  Detlef Steuer    <detlef.steuer@hsu-hamburg.de>
 *  Olaf Mersmann    <olafm@statistik.uni-dortmund.de>
 *
 *  Adapted by Emanuele Aliverti on Oct 2024
 */

//#include <R.h>
//#include <Rmath.h>
//#include <Rinternals.h>
#include <RcppArmadillo.h>

#include "truncnorm.h"

/*
 * These routines calculate the expected value and variance of the
 * left, right and doubly truncated normal distribution. The only
 * tricky bit is the calculation of the variance of the doubly
 * truncated normal distribution. We use a decompostion of the
 * variance of a mixture of distributions to here for numerical
 * reasons. For details see:
 *
 *   Foulley JL. A completion simulator for the two-sided truncated
 *   normal distribution. Genetics, selection, evolution 2000
 *   Nov-Dec;32(6): p. 631-635.
 */
double e_lefttruncnorm(double a, double mean, double sd) {
  const double alpha = (a - mean) / sd;
  const double phi_a = R::dnorm(alpha, 0.0, 1.0, TRUE);
  const double Phi_a = R::pnorm(alpha, 0.0, 1.0, FALSE, TRUE);
  double res = mean + sd * exp(phi_a - Phi_a);
  return res;
}

double e_truncnorm(double a, double b, double mean, double sd) {
  /* Special case numerically instable case when (a, b) is far away from the
   * center of mass. */
  if (b < mean - 6.0 * sd || a > mean + 6.0 * sd)
    return (a + b) / 2.0;

  double delta_phi = 0.0, delta_Phi = 0.0;

  const double alpha = (a - mean) / sd;
  const double beta = (b - mean) / sd;
  const double phi_a = R::dnorm(alpha, 0.0, 1.0, TRUE);
  const double Phi_a = R::pnorm(alpha, 0.0, 1.0, TRUE, TRUE);
  const double phi_b = R::dnorm(beta, 0.0, 1.0, TRUE);
  const double Phi_b = R::pnorm(beta, 0.0, 1.0, TRUE, TRUE);

  if (phi_b < phi_a) {
    delta_phi = R::logspace_sub(phi_a, phi_b);
  } else {
    sd = -sd;
    delta_phi = R::logspace_sub(phi_b, phi_a);
  }

  if (Phi_b > Phi_a) {
    sd = -sd;
    delta_Phi = R::logspace_sub(Phi_b, Phi_a);
  } else {
    delta_Phi = R::logspace_sub(Phi_a, Phi_b);
  }
  return mean + sd * -exp(delta_phi - delta_Phi);
}

double e_righttruncnorm(double b, double mean, double sd) {
  const double beta = (b - mean) / sd;
  const double phi_b = R::dnorm(beta, 0.0, 1.0, TRUE);
  const double Phi_b = R::pnorm(beta, 0.0, 1.0, TRUE, TRUE);
  return mean + sd * -exp(phi_b - Phi_b);
}

double v_lefttruncnorm(double a, double mean, double sd) {
  const double alpha = (a - mean) / sd;
  const double phi_a = R::dnorm(alpha, 0.0, 1.0, FALSE);
  const double Phi_a = R::pnorm(alpha, 0.0, 1.0, TRUE, FALSE);
  const double lambda = phi_a / (1.0 - Phi_a);
  return (sd * sd * (1.0 - lambda * (lambda - alpha)));
}

 double v_righttruncnorm(double b, double mean, double sd) {
  return (v_lefttruncnorm(-b, -mean, sd));
}

 double v_truncnorm(double a, double b, double mean, double sd) {
  /* Special case numerically instable cases. These arise when (a, b) is far
   * away from mean +/- 6*sd */
  if (b < mean - 6.0 * sd || a > mean + 6.0 * sd)
    return 1.0 / 12 * (b - a) * (b - a);

  const double v = sd * sd;
  const double pi1 = R::pnorm(a, mean, sd, TRUE, FALSE);
  const double pi2 =
      R::pnorm(b, mean, sd, TRUE, FALSE) - R::pnorm(a, mean, sd, TRUE, FALSE);
  const double pi3 = R::pnorm(b, mean, sd, FALSE, FALSE); /* 1 - F(b) */
  const double e1 = e_righttruncnorm(a, mean, sd);
  const double e2 = e_truncnorm(a, b, mean, sd);
  const double e3 = e_lefttruncnorm(b, mean, sd);

  const double v1 = v_righttruncnorm(a, mean, sd);
  const double v3 = v_lefttruncnorm(b, mean, sd);

  const double c1 = pi1 * (v1 + (e1 - mean) * (e1 - mean));
  const double c3 = pi3 * (v3 + (e3 - mean) * (e3 - mean));
  //const double cd = pi2 - (e2 - mean) * (e2 - mean);
  return (v - c1 - c3) / pi2 - (e2 - mean) * (e2 - mean);
}

