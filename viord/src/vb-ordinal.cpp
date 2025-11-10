/*
 * VB for ordinal probit regression
 * Author: Emanuele Aliverti
 *
 */

#include <RcppArmadillo.h>
#include "truncnorm.h"
#include "utils.h"  // include the utilities
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// MAIN VB ROUTINE
// [[Rcpp::export]]
Rcpp::List vb_ordinal(
		const arma::vec& Y,  // Response variable
		const arma::mat& X, // Covariates
		const arma::vec& alpha, // (tresholds)
		const arma::vec& mu0, // prior mean
		const arma::mat& S0, // prior variance
		const arma::mat& Q0, // prior precision
		const int maxit = 100, // max number of iterations
		const double tresh = 1e-6, // tolerance
		const bool verbose=false, // print info
		const bool full_out=false // what is returned as output
		)
{
	int n = X.n_rows;
	int p = X.n_cols;

	// Convert tresholds into limits of integration, based on y
	arma::vec u(n);
	arma::vec l(n);

	arma::vec Ui(n);
	arma::vec Li(n);

	arma::ivec y_id = conv_to<ivec>::from(Y); // index from 1:K

	for(int i = 0; i < n; i++) {
		// Index from 0
		u(i) = alpha(y_id(i));
		l(i) = alpha(y_id(i)-1);
	}

	arma::vec z(n, fill::zeros); // augmented data
	arma::vec mu_beta(p, fill::zeros); //posterior mean of regression coefficients

	arma::vec elbo_seq(maxit + 1);
	int it = 0;
	bool conv = false;
	double lp;
	//

	// precompute X'X + Q0
	arma::mat V(p,p);
	arma::inv(V, (X.t() * X + Q0), arma::inv_opts::allow_approx);
	arma::mat VXt(p,n);
	VXt = V * X.t();
	arma::vec Vprior(p);
	Vprior = V*Q0*mu0;


	while (!conv && it < maxit) {
		it++;
		Rcpp::checkUserInterrupt();


		// update z
		z = X * mu_beta;
		Li = l - z;
		Ui = u - z;

		lp = 0.0;
		for(int i = 0; i < n; i++){
			z(i) -= zeta1(Li(i), Ui(i));
			lp += trunc_log(R::pnorm(Ui(i), 0.0, 1.0, TRUE, FALSE) - R::pnorm(Li(i), 0.0, 1.0, TRUE, FALSE));
		}
		
		// update  mu_beta
		mu_beta = Vprior + VXt * z;

		elbo_seq(it) = elbo_mf(mu_beta, Q0, mu0, lp);

		conv = ( abs(elbo_seq(it) - elbo_seq(it-1)) < tresh);

	}

	// Output
	Rcpp::List out;
	if(full_out) {
		out["m"] = mu_beta;
		out["S"] = V;
		out["it"] = it;
		out["elbo_seq"] = elbo_seq.subvec(1,it);
		out["elbo"] = elbo_seq(it);
		out["conv"] = conv;
		//elbo

	} else {
		out["m"] = mu_beta;
		out["S"] = V;
		out["it"] = it;
	}

	return  out;
}
