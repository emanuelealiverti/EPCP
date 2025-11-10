/*
 * VB - PMF for ordinal probit regression
 * Author:
 * Emanuele Aliverti
 *
 */

#include <RcppArmadillo.h>
#include "truncnorm.h"
#include "utils.h"  // include the utilities

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
Rcpp::List pmf_ordinal(
		const arma::vec& Y,  // Response variable
		const arma::mat& X, // Covariates
		const arma::vec& alpha, // (tresholds)
		const arma::vec& mu0, // prior mean
		const arma::mat& S0, // prior variance
		const arma::mat& Q0, // prior precision
		const int maxit = 100, // max number of iterations
		const double tresh = 1e-6, // tolerance
		const bool verbose=false, // print info
		const bool full_out=false // return elbo, z, etc
		)
{
	int n = X.n_rows;
	int p = X.n_cols;

	// Convert tresholds into limits of integration, based on y
	arma::vec u(n);
	arma::vec l(n);

	// i-th element used within parameter updates
	double Ui;
	double Li;

	arma::ivec y_id = conv_to<ivec>::from(Y); // index from 1:K

	for(int i = 0; i < n; i++) {
		// Index from 0
		u(i) = alpha(y_id(i));
		l(i) = alpha(y_id(i)-1);
	}

	// Quantities to monitor convergence
	arma::vec elbo_seq(maxit + 1);
	int it = 0;
	bool conv = false;
	double lp = 0.0;
	
	// X\mu_0
	arma::vec Xmu0(n);
	Xmu0 = X * mu0;



	arma::mat V(p,p);
	arma::inv(V, (X.t() * X + Q0), arma::inv_opts::allow_approx);
	arma::mat XV(n,p);
	XV = X*V;
	
	arma::vec sigma2Z(n);
	//arma::vec H_diag(n);
	sigma2Z = 1.0 / (1 - sum((XV)%X,1));

	arma::vec sigmaZ(n);
	sigmaZ = arma::sqrt(sigma2Z);

	arma::vec xiZ(n, fill::zeros);
	arma::vec meanZ(n, fill::zeros);
	arma::mat D(n,p, fill::zeros);
	arma::vec sigma_xv_prior(n);
	arma::colvec di_Xmu(p); // used to include prior contribution in d
	// initialize with the prior contribution Xmi'Xmi mu0 and then gets updated recursively
	arma::mat XtX(p,p);
	XtX = X.t() * X;
	for(int i = 0; i < n; i++){
		sigma_xv_prior(i) =  sigma2Z(i) * as_scalar(XV.row(i) * (XtX - (X.row(i).t() * X.row(i))) * mu0); 
	}
//	arma::vec varZ(n, fill::zeros); //for elbo only

	int im1;
	while (!conv && it < maxit) {
		it++;
		Rcpp::checkUserInterrupt();

		lp = 0.0;
		for(int i = 0; i < n; i++){
			// need the index (i-1) that should be n when i = 0
			im1 = (i == 0) ? (n-1) : (i-1);
			D.row(i) = D.row(im1) - X.row(i) * meanZ(i) + X.row(im1) * meanZ(im1);
			xiZ(i) = Xmu0(i) - sigma_xv_prior(i) + sigma2Z(i) * arma::as_scalar(XV.row(i) * D.row(i).t());

			Li = (l(i) - xiZ(i))/sigmaZ(i);
			Ui = (u(i) - xiZ(i))/sigmaZ(i);
			meanZ(i) = xiZ(i) - sigmaZ(i) * zeta1(Li,Ui);

			lp += trunc_log(R::pnorm(Ui, 0.0, 1.0, TRUE, FALSE) - R::pnorm(Li, 0.0, 1.0, TRUE, FALSE));
		}

		//elbo_seq(it) = elbo(In_m_H, xiZ, sigma2Z, meanZ,  Xmu0, lp);
		elbo_seq(it) = elbo_pmf(X, V, XV, xiZ, sigma2Z, meanZ, Xmu0, mu0, lp);
		conv = ( abs(elbo_seq(it) - elbo_seq(it-1)) < tresh);
	}

	arma::vec sdZ(n);
	for(int i = 0; i < n; i++){

		im1 = (i == 0) ? (n-1) : (i-1);
		D.row(i) = D.row(im1) - X.row(i) * meanZ(i) + X.row(im1) * meanZ(im1);
		xiZ(i) = Xmu0(i) - sigma_xv_prior(i) + sigma2Z(i) * arma::as_scalar(XV.row(i) * D.row(i).t());
		Li = (l(i) - xiZ(i))/sigmaZ(i);
		Ui = (u(i) - xiZ(i))/sigmaZ(i);
		meanZ(i) = xiZ(i) - sigmaZ(i) * zeta1(Li,Ui);
		sdZ(i) = sigmaZ(i) * std::sqrt((1-z1p2_z2(Li,Ui)));
	}

	arma::vec meanBeta(p);
	arma::mat varBeta(p,p);

	meanBeta = XV.t() * meanZ + (V * Q0 * mu0);
	//arma::mat XVw(n,p);
	//XVw = X*V;
	XV.each_col() %= sdZ;
	varBeta = V + XV.t() * XV;

	// Output
	Rcpp::List out;
	if(full_out) {
		// marginal moments for \beta
		out["m"] = meanBeta;
		out["S"] = varBeta;
		// Quantities for truncated normal simulations
		out["sigmaZ"] = sigmaZ;
		out["xiZ"] = xiZ;
		out["elbo_seq"] = elbo_seq.subvec(1,it);
		out["elbo"] = elbo_seq(it);
		out["conv"] = conv;
		out["it"] = it;

	} else {
		out["m"] = meanBeta;
		out["S"] = varBeta;
	}

	return  out;
}
