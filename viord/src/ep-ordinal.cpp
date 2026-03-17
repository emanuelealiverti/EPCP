/*
 * EP for ordinal probit regression
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


// MAIN EP ROUTINE
// [[Rcpp::export]]
Rcpp::List ep_ordinal(
		const arma::vec& Y,  // Response variable
		const arma::mat& X, // Covariates
		const arma::vec& alpha, // (tresholds)
		const arma::vec& mu0, // prior mean
		const arma::mat& S0, // prior variance
		const arma::mat& Q0, // prior precision
		const int maxit = 100, // max number of iterations
		const double tresh = 1e-6, // tolerance
		const bool verbose=false, // print information
		const bool full_out=false // what is returned as output
		)
{
	int n = X.n_rows;
	int p = X.n_cols;

	// Convert tresholds into limits of integration, based on y
	arma::vec u(n);
	arma::vec l(n);
	arma::ivec y_id = conv_to<ivec>::from(Y); // index from 1:K

	for(int i = 0; i < n; i++) {
		// Index from 0
		u(i) = alpha(y_id(i));
		l(i) = alpha(y_id(i)-1);
	}



	arma::vec logZ(n, arma::fill::zeros); // marginal likelihood

	// EP parameters
	// initialized from the prior
	arma::mat S_ep = S0;
//	arma::mat Q_ep(p,p);
//	Q_ep = Q0; 

	arma::vec r_ep = S_ep*mu0;
	arma::vec S_ep_xi(p); // S_ep x_i

	//double z0 = lPsi(r_ep,Q_ep); // since r_ep = r_0 and Q_ep = Q0, this is the prior contribution 
	double z0 = lPsi_cov(r_ep, S_ep); // since r_ep = r_0 and Q_ep = Q0, this is the prior contribution 
	
	// Remove 2pi factor since does not depend on pars and simplifies at each step
	z0 -=  0.5*p*log(2*M_PI);
	
	// initialize log-determinat of precision (compute from variance and change sign)
	double logd_Qep;
	log_det_sympd(logd_Qep, S_ep);
	logd_Qep = -1.0*logd_Qep;


	// Cavity parameters
	arma::mat Si(p,p, arma::fill::eye);
	arma::vec r_i(p, arma::fill::zeros); 
	arma::vec mu_i(p);


	// Functionals of the cavity
	arma::vec Si_xi(p);    // R_i x_i
	double xi_Si_xi; // x_i' R_i x_i
	double xi_Si_mi; // x_i' R_i m_i


	// Working quantities
	// updates for S_ep / Q_ep and r_ep
	arma::vec k(n, arma::fill::zeros);
	arma::vec w(n, arma::fill::zeros);
	
	// Working scalars
	double si = 0.0; // 1.0 + xi_Si_xi
	double sisq = 0.0; // sqrt(1.0 + xi_Si_xi)
	double di = 0.0; 
	double z = 0.0; 
	double z_ep = 0.0;
	double z_ep_old = 0.0; // check convergence
	double k_old = 0.0;
	double Ui; // lower limit of integration
	double Vi; // upper limit of integration


	double err = -99;


	int it = 0;
	bool conv = false;
	// sequence of p(y,q)
	arma:: vec z_seq(maxit);

	while (!conv && it < maxit) {
		z_ep_old = z_ep;

		for(int i = 0; i < n; i++) {
			Rcpp::checkUserInterrupt();


			S_ep_xi = S_ep * X.row(i).t();

			di = k(i) /  (1.0 - k(i) * arma::as_scalar(X.row(i) * S_ep_xi) );
			if(std::isinf(di)) di = 0.0;

			// Compute parameters of the cavity distribution
			Si = S_ep + di * (S_ep_xi) * (S_ep_xi).t();
			r_i = r_ep - w(i) * X.row(i).t();
			mu_i = Si * r_i;

			// Functionals of the cavity
			Si_xi = Si * X.row(i).t();
			xi_Si_xi = arma::as_scalar(X.row(i) * Si_xi);
			xi_Si_mi = arma::as_scalar(X.row(i) * Si * r_i);

			// Moments of the hybrid
			si = 1.0 + xi_Si_xi;
			sisq = std::sqrt(si);

			if(si > 0.0) {

				// Truncation values
				Ui = (l(i) - xi_Si_mi) / sisq;
				Vi = (u(i) - xi_Si_mi) / sisq;

				z = z1p2_z2(Ui, Vi);

				k_old = k(i); // store old value for logD
				k(i) = z / (1.0 + xi_Si_xi * (1.0 - z));
//				if(std::isnan(k(i)) || isinf(k(i))) k(i) = 0.0;
				//Rcout << "k(i): " << k(i) << std::endl;


				//double dup = trunc_log(1 + (k(i) - k_old) * arma::as_scalar(X.row(i) * S_ep * X.row(i).t()));

				w(i) = xi_Si_mi * k(i) - zeta1(Ui,Vi) * sisq;

				// update determinant
				double dup = trunc_log(1 + (k(i) - k_old) * arma::as_scalar(X.row(i) * S_ep * X.row(i).t()));
				if(std::isnan(dup)) dup = 0.0;
				logd_Qep += dup; 

				// Update moments
				S_ep = Si - (Si_xi * Si_xi.t()) * (z / si);
				r_ep = r_i + w(i) * X.row(i).t();

				// Update marginal likelihood
				logZ(i) = -trunc_log(R::pnorm(Vi, 0.0, 1.0, TRUE, FALSE) - R::pnorm(Ui, 0.0, 1.0, TRUE, FALSE));
				//Z(i) = -trunc_log(arma::normcdf(Vi) - arma::normcdf(Ui));
				di = 1.0 / (1.0 + k(i) * xi_Si_xi) * (2*w(i) * xi_Si_mi + (w(i) * w(i)) * xi_Si_xi - k(i) * xi_Si_mi * xi_Si_mi);
				// check before adding for large p
				logZ(i) += (0.5 * di - 0.5 *trunc_log(1+k(i)*xi_Si_xi));
			}


		}


		// compute log p_ep(y)
		z_ep = 0.5 * arma::as_scalar(r_ep.t() * S_ep * r_ep) - 0.5 * logd_Qep;
		z_ep -= (z0 + sum(logZ));

		// (alternatively, check convergenge of posterior means)
		err = arma::sum(std::abs(z_ep_old - z_ep));
		conv = (err < tresh);
		z_seq(it) = z_ep;
		it++;

		if( (verbose) & (it % 10 == 0) ) {
			Rcout << "it: " << it << " err: " << err << std::endl;
		}
	}

	arma::vec m_ep(p);
	m_ep = S_ep*r_ep;
	
	// Output
	Rcpp::List out;
	if(full_out) {
		out["S"] = S_ep;
		out["m"] = m_ep;
		out["it"] = it;
		out["logZ_seq"] = z_seq.subvec(0,it-1);
		out["w"] = w.col(0); // return (nx1) matr as vectors
		out["k"] = k.col(0);
		out["Z"] = logZ.col(0);
		out["logd_Qep"] = logd_Qep;
		out["logZ"] = z_ep;

	} else {
		// Compute marginal likelihood combining all contributions
		//Z.replace(datum::inf, 0.0);
		out["S"] = S_ep;
		out["m"] = m_ep;
		out["logZ"] = z_ep;
		out["it"] = it;
	}

	return  out;
}
