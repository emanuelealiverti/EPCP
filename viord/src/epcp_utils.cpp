/*
 * Approximate Bayesian inference for ordinal probit regression
 * Author: Emanuele Aliverti, aliverti@stat.unipd.it
 * This file contains serveral utilities for computing the proposed algorithm 
 *
 * - utilities for computing zeta1 and zeta2 functions are based on R package truncnorm (https://github.com/olafmersmann/truncnorm/blob/master/src/truncnorm.c)
 */




#include "utils.h"
#include "truncnorm.h"
/*
 * zeta1 and zeta1^2 + zeta2 correspond to mean (changed of sign) and 1-var of truncated normal distribution
 * Both functions are be made numerically more robust when x is very large
 * relying on asymtotic expansions (see zeta0)
 */

// utilities for comuting truncated normal moments
double zeta0(double& x) {
	// asymptotic expansion (26.2.13) of Abramowitz and Stegun (1964) for very large values of x
	// Taken from sn pacakge by Azzalini (https://cran.r-project.org/web/packages/sn/index.html) for mills ratio
	double x2 = x * x;
	double res = -x / (1 - 1 / (x2 + 2) + 1 / ((x2 + 2) * (x2 + 4)) 
                        - 5 / ((x2 + 2) * (x2 + 4) * (x2 + 6))
                        + 9 / ((x2 + 2) * (x2 + 4) * (x2 + 6) * (x2 + 8)) 
                        - 129 / ((x2 + 2) * (x2 + 4) * (x2 + 6) * (x2 + 8) * (x2 + 10)));
	return -1.0*res;
}

double zeta1(double ca, double cb) {
	double ret;
	if (R_FINITE(ca) && R_FINITE(cb)) {
		ret = e_truncnorm(ca, cb, 0.0, 1.0);
	} else if (R_NegInf == ca && R_FINITE(cb)) {
		if(cb < -50) ret = zeta0(cb);
		ret = e_righttruncnorm(cb, 0.0, 1.0);
	} else if (R_FINITE(ca) && R_PosInf == cb) {
		if(ca > 50) ret = zeta0(ca);
		ret = e_lefttruncnorm(ca, 0.0, 1.0);
	} else if (R_NegInf == ca && R_PosInf == cb) {
		ret = 0.0;
	} else {
		ret = NA_REAL;
	}

	return -1.0*ret;
}


double z1p2_z2(double ca, double cb) {
	//zeta_1 ^2 + zeta_2
	double ret;
	if (R_FINITE(ca) && R_FINITE(cb)) {
		ret = v_truncnorm(ca, cb, 0.0, 1.0);
	} else if (R_NegInf == ca && R_FINITE(cb)) {
		double z = zeta1(ca,cb);
		ret = 1.0 - z*(z+cb);
		if(std::isnan(ret)) ret = 0.0;
	} else if (R_FINITE(ca) && R_PosInf == cb) {
		double z = zeta1(ca,cb);
		ret = 1.0 - z*(z+ca);
		if(std::isinf(ret)) ret = 0.0;
	} else if (R_NegInf == ca && R_PosInf == cb) {
		ret = 1.0; // unbounded
	} else { ret = NA_REAL;
	}
	return 1.0-ret;
}


// ELBOs

double elbo_mf(arma::vec mu_b, arma::mat Q0, arma::mat mu0, double lp){
	// consider only quantites that depend on parameters
	double out = 0.0;
	out = lp - 0.5 * as_scalar(mu_b.t() * Q0 * mu_b) + as_scalar(mu_b.t()*Q0*mu0);
	return out;
}

double elbo_pmf(arma::mat X, arma::mat V, arma::mat XV,  arma::vec xiZ, arma::vec sigma2Z, arma::vec meanZ, arma::vec Xmu0, arma::vec mu0, double lp){
	double out = 0.0;
	arma::rowvec zX(X.n_cols);
	arma::vec zz(X.n_rows);
	zX = meanZ.t() * X; 
	zz = meanZ % meanZ;
	out += -0.5 * accu(zz) + 0.5 * as_scalar(zX * V * zX.t());
	out += as_scalar(zX * mu0) - as_scalar(zX * XV.t() * Xmu0);

	// E_p[q(z)]
	out += 0.5 * accu(zz / sigma2Z) - accu(meanZ % xiZ / sigma2Z) + 0.5 * accu(xiZ % xiZ / sigma2Z);
	out += lp;
	return out;
}


double log2pi = log(2*M_PI);

// Gaussian normalizing constant, used for computing the marginal likelihood for EP
double lPsi(const arma::vec& r, const arma::mat& Q) {
	// Implement lPsi function
	double ldet, res;
	int p = Q.n_cols;
	log_det_sympd(ldet, Q);
	res = 0.5 * arma::as_scalar(r.t() * arma::inv(Q) * r) - 0.5 * ldet + 0.5*p*log2pi;
	return res;
}

// Gaussian normalizing constant, in terms of covariance matrix (avoids a further inversion)
double lPsi_cov(const arma::vec& r, const arma::mat& S) {
	// Implement lPsi function
	double ldet, res;
	int p = S.n_cols;
	log_det_sympd(ldet, S);
	res = 0.5 * arma::as_scalar(r.t() * S * r) + 0.5 * ldet + 0.5*p*log2pi;
	return res;
}


