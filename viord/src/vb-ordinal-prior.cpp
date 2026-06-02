/*
 * VB for ordinal probit regression with variance-component priors.
 *
 * beta | sigma_b2 ~ N(mu0, sigma_b2 I_p)
 * u_g  | sigma_u2_g ~ N(0, sigma_u2_g I)
 *
 * sigma_b2 ~ IG(a0, b0)
 * sigma_u2_g ~ IG(au0, bu0)
 *
 */

#include <RcppArmadillo.h>
#include "truncnorm.h"
#include "utils.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

namespace {

const double log2pi_prior = std::log(2.0 * M_PI);

double trace_product(const arma::mat& A, const arma::mat& B) {
	return arma::accu(A % B.t());
}

void inv_sympd_allow_approx(arma::mat& out, const arma::mat& in) {
	if(!arma::inv_sympd(out, in)) {
		arma::inv(out, in, arma::inv_opts::allow_approx);
	}
}

double expected_log_ig(const double a, const double b,
                       const double E_log_sigma2,
                       const double E_inv_sigma2) {
	return a * std::log(b) - R::lgammafn(a) -
		(a + 1.0) * E_log_sigma2 - b * E_inv_sigma2;
}

double elbo_mf_ordinal_vc(
		const arma::mat& W,
		const arma::mat& WtW,
		const arma::vec& mean_a,
		const arma::vec& location_a,
		const arma::vec& mu_theta,
		const arma::mat& S_theta,
		const arma::vec& mu0,
		const arma::uvec& u_group,
		const arma::uvec& group_counts,
		const double a_sigma_b2,
		const double b_sigma_b2,
		const arma::vec& a_sigma_u2,
		const arma::vec& b_sigma_u2,
		const double a0,
		const double b0,
		const double au0,
		const double bu0,
		const double lp
		)
{
	const int p = mu0.n_elem;
	const int q = u_group.n_elem;
	const int d = mu_theta.n_elem;
	const int n_groups = group_counts.n_elem;

	const double tau_b = a_sigma_b2 / b_sigma_b2;
	const double E_log_sigma_b2 = std::log(b_sigma_b2) - R::digamma(a_sigma_b2);

	double logdetS = 0.0;
	if(!arma::log_det_sympd(logdetS, S_theta)) {
		double sign = 0.0;
		arma::log_det(logdetS, sign, S_theta);
	}

	arma::vec current_location = W * mu_theta;

	double out = lp;

	// E[log p(y, a | theta)] - E[log q(a)].
	out += arma::dot(mean_a, current_location - location_a);
	out -= 0.5 * (arma::dot(current_location, current_location) +
	              trace_product(WtW, S_theta) -
	              arma::dot(location_a, location_a));

	// Fixed-effect prior contribution.
	arma::vec beta_shift = mu_theta.head(p) - mu0;
	const double E_beta_ss = arma::dot(beta_shift, beta_shift) +
		arma::trace(S_theta.submat(0, 0, p - 1, p - 1));

	out += -0.5 * p * log2pi_prior;
	out += -0.5 * p * E_log_sigma_b2;
	out += -0.5 * tau_b * E_beta_ss;
	out += expected_log_ig(a0, b0, E_log_sigma_b2, tau_b);

	// Random-effect prior contributions.
	for(int g = 0; g < n_groups; g++) {
		const double tau_u = a_sigma_u2(g) / b_sigma_u2(g);
		const double E_log_sigma_u2 = std::log(b_sigma_u2(g)) - R::digamma(a_sigma_u2(g));
		double E_u_ss = 0.0;

		for(int j = 0; j < q; j++) {
			if(u_group(j) == static_cast<arma::uword>(g)) {
				const int idx = p + j;
				E_u_ss += mu_theta(idx) * mu_theta(idx) + S_theta(idx, idx);
			}
		}

		out += -0.5 * group_counts(g) * log2pi_prior;
		out += -0.5 * group_counts(g) * E_log_sigma_u2;
		out += -0.5 * tau_u * E_u_ss;
		out += expected_log_ig(au0, bu0, E_log_sigma_u2, tau_u);
	}

	// -E[log q(theta)].
	out += 0.5 * (d * (1.0 + log2pi_prior) + logdetS);

	// -E[log q(sigma_b2)].
	out -= expected_log_ig(a_sigma_b2, b_sigma_b2, E_log_sigma_b2, tau_b);

	// -E[log q(sigma_u2_g)].
	for(int g = 0; g < n_groups; g++) {
		const double tau_u = a_sigma_u2(g) / b_sigma_u2(g);
		const double E_log_sigma_u2 = std::log(b_sigma_u2(g)) - R::digamma(a_sigma_u2(g));
		out -= expected_log_ig(a_sigma_u2(g), b_sigma_u2(g),
		                       E_log_sigma_u2, tau_u);
	}

	return out;
}

} // namespace

// MAIN VB ROUTINE
// [[Rcpp::export]]
Rcpp::List vb_ordinal_prior(
		const arma::vec& Y,  // Response variable
		const arma::mat& X, // Fixed-effect covariates
		const arma::vec& alpha, // thresholds
		const arma::vec& mu0, // fixed-effect prior mean
		const double a0, // Inv-Gamma prior shape for sigma_b2
		const double b0, // Inv-Gamma prior scale for sigma_b2
		const arma::mat& Z, // random-effect design; ignored when it has one column
		const arma::uvec& Z_group, // zero-based variance group per Z column
		const double au0 = NA_REAL, // Inv-Gamma prior shape for sigma_u2
		const double bu0 = NA_REAL, // Inv-Gamma prior scale for sigma_u2
		const int maxit = 100, // max number of iterations
		const double tresh = 1e-6, // tolerance
		const bool verbose=false, // print info
		const bool full_out=false // what is returned as output
		)
{
	const int n = X.n_rows;
	const int p = X.n_cols;

	if(p < 1) {
		Rcpp::stop("X must have at least one column.");
	}
	if(Y.n_elem != static_cast<unsigned int>(n)) {
		Rcpp::stop("Y and X have incompatible dimensions.");
	}
	if(mu0.n_elem != static_cast<unsigned int>(p)) {
		Rcpp::stop("mu0 must have length equal to ncol(X).");
	}
	if(a0 <= 0.0 || b0 <= 0.0) {
		Rcpp::stop("a0 and b0 must be strictly positive.");
	}
	if(maxit < 1) {
		Rcpp::stop("maxit must be positive.");
	}

	const bool has_random = (Z.n_cols > 1);
	const int q = has_random ? Z.n_cols : 0;
	arma::uvec u_group; // group indicators
	arma::uvec group_counts; // number of columns of Z per group
	if(has_random) {
		u_group = Z_group;
		const int n_groups = arma::max(u_group) + 1;
		group_counts.set_size(n_groups);
		group_counts.zeros();
		for(int j = 0; j < q; j++) {
			group_counts(u_group(j))++;
		}
	} else {
		group_counts.set_size(0);
	}

	// Convert thresholds into limits of integration, based on y.
	arma::vec u_bound(n);
	arma::vec l_bound(n);
	arma::ivec y_id = conv_to<ivec>::from(Y); // index from 1:K

	for(int i = 0; i < n; i++) {
		u_bound(i) = alpha(y_id(i));
		l_bound(i) = alpha(y_id(i)-1);
	}

	arma::mat W;
	if(has_random) {
		W = arma::join_rows(X, Z);
	} else {
		W = X;
	}
	const int d = W.n_cols;
	arma::mat WtW = W.t() * W;

	arma::vec prior_mean(d, fill::zeros);
	prior_mean.head(p) = mu0;

	arma::vec a(n, fill::zeros); // augmented data posterior means
	arma::vec a_location(n, fill::zeros); // location used by q(a)
	arma::vec mu_theta = prior_mean; // posterior mean of fixed and random coefficients
	arma::mat S_theta(d, d, fill::eye); // posterior covariance
	arma::vec Ui(n);
	arma::vec Li(n);

	// Initialization
	double a_sigma_b2 = a0 + 0.5 * p;
	double b_sigma_b2 = b0;
	double tau_b = a0 / b0;

	const int n_groups = group_counts.n_elem;
	arma::vec a_sigma_u2(n_groups, fill::zeros);
	arma::vec b_sigma_u2(n_groups, fill::zeros);
	arma::vec tau_u(n_groups, fill::zeros);
	for(int g = 0; g < n_groups; g++) {
		a_sigma_u2(g) = au0 + 0.5 * group_counts(g); // set in the full conditional
		b_sigma_u2(g) = bu0;
		tau_u(g) = au0 / bu0;
	}

	arma::vec prior_precision_diag(d, fill::zeros);
	arma::vec elbo_seq(maxit, fill::zeros);
	arma::vec b_sigma_b2_seq(maxit, fill::zeros);
	arma::vec tau_b_seq(maxit, fill::zeros);
	arma::mat b_sigma_u2_seq(maxit, n_groups, fill::zeros);
	arma::mat tau_u_seq(maxit, n_groups, fill::zeros);
	int it = 0;
	bool conv = false;
	double lp = 0.0;

	while (!conv && it < maxit) {
		Rcpp::checkUserInterrupt();

		// Update q(a).
		a_location = W * mu_theta;
		a = a_location;
		Li = l_bound - a_location;
		Ui = u_bound - a_location;

		lp = 0.0;
		for(int i = 0; i < n; i++) {
			a(i) -= zeta1(Li(i), Ui(i));
			lp += trunc_log(R::pnorm(Ui(i), 0.0, 1.0, TRUE, FALSE) -
			                R::pnorm(Li(i), 0.0, 1.0, TRUE, FALSE));
		}

		// Update q(theta) using E[1 / sigma_b2] and E[1 / sigma_u2_g].
		prior_precision_diag.head(p).fill(tau_b);
		for(int j = 0; j < q; j++) {
			prior_precision_diag(p + j) = tau_u(u_group(j));
		}

		arma::mat theta_precision = WtW;
		theta_precision.diag() += prior_precision_diag;
		inv_sympd_allow_approx(S_theta, theta_precision);
		mu_theta = S_theta * (W.t() * a + prior_precision_diag % prior_mean);

		// Update q(sigma_b2).
		arma::vec beta_shift = mu_theta.head(p) - mu0;
		b_sigma_b2 = b0 + 0.5 * (
			arma::dot(beta_shift, beta_shift) +
			arma::trace(S_theta.submat(0, 0, p - 1, p - 1)));
		tau_b = a_sigma_b2 / b_sigma_b2;

		// Update q(sigma_u2_g).
		for(int g = 0; g < n_groups; g++) {
			double E_u_ss = 0.0;
			for(int j = 0; j < q; j++) {
				if(u_group(j) == static_cast<arma::uword>(g)) {
					const int idx = p + j;
					E_u_ss += mu_theta(idx) * mu_theta(idx) + S_theta(idx, idx);
				}
			}
			b_sigma_u2(g) = bu0 + 0.5 * E_u_ss;
			tau_u(g) = a_sigma_u2(g) / b_sigma_u2(g);
		}

		elbo_seq(it) = elbo_mf_ordinal_vc(
			W, WtW, a, a_location, mu_theta, S_theta, mu0, u_group,
			group_counts, a_sigma_b2, b_sigma_b2, a_sigma_u2,
			b_sigma_u2, a0, b0, au0, bu0, lp);
		b_sigma_b2_seq(it) = b_sigma_b2;
		tau_b_seq(it) = tau_b;
		if(n_groups > 0) {
			b_sigma_u2_seq.row(it) = b_sigma_u2.t();
			tau_u_seq.row(it) = tau_u.t();
		}

		if(it > 0) {
			conv = (std::abs(elbo_seq(it) - elbo_seq(it-1)) < tresh);
		}

		it++;

		if(verbose && (it % 10 == 0 || conv)) {
			Rcout << "it: " << it
			      << " elbo: " << elbo_seq(it-1)
			      << " E_inv_sigma_b2: " << tau_b
			      << std::endl;
		}
	}

	const double sigma_b2_mean =
		(a_sigma_b2 > 1.0) ? b_sigma_b2 / (a_sigma_b2 - 1.0) : R_PosInf;
	arma::vec sigma_u2_mean(n_groups, fill::zeros);
	for(int g = 0; g < n_groups; g++) {
		sigma_u2_mean(g) =
			(a_sigma_u2(g) > 1.0) ? b_sigma_u2(g) / (a_sigma_u2(g) - 1.0) : R_PosInf;
	}

	arma::vec m_beta = mu_theta.head(p);
	arma::mat S_beta = S_theta.submat(0, 0, p - 1, p - 1);
	arma::vec m_u;
	arma::mat S_u;
	if(q > 0) {
		m_u = mu_theta.subvec(p, d - 1);
		S_u = S_theta.submat(p, p, d - 1, d - 1);
	} else {
		m_u.set_size(0);
		S_u.set_size(0, 0);
	}

	// Output
	Rcpp::List out;
	out["m"] = m_beta;
	out["S"] = S_beta;
	out["m_u"] = m_u;
	out["S_u"] = S_u;
	out["m_joint"] = mu_theta;
	out["S_joint"] = S_theta;
	out["it"] = it;
	out["sigma_b2_a"] = a_sigma_b2;
	out["sigma_b2_b"] = b_sigma_b2;
	out["sigma_b2_mean"] = sigma_b2_mean;
	out["sigma_b2_inv_mean"] = tau_b;
	out["sigma_u2_a"] = a_sigma_u2;
	out["sigma_u2_b"] = b_sigma_u2;
	out["sigma_u2_mean"] = sigma_u2_mean;
	out["sigma_u2_inv_mean"] = tau_u;
	arma::uvec group_index;
	if(n_groups > 0) {
		group_index = arma::regspace<arma::uvec>(0, n_groups - 1);
	} else {
		group_index.set_size(0);
	}
	out["u_group"] = group_index;

	if(full_out) {
		out["a"] = a;
		out["elbo_seq"] = elbo_seq.subvec(0, it-1);
		out["elbo"] = elbo_seq(it-1);
		out["sigma_b2_b_seq"] = b_sigma_b2_seq.subvec(0, it-1);
		out["sigma_b2_inv_mean_seq"] = tau_b_seq.subvec(0, it-1);
		if(n_groups > 0) {
			out["sigma_u2_b_seq"] = b_sigma_u2_seq.rows(0, it-1);
			out["sigma_u2_inv_mean_seq"] = tau_u_seq.rows(0, it-1);
		}
		out["conv"] = conv;
	}

	return  out;
}
