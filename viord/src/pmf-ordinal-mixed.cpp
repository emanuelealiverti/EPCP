/*
 * VB-PMF for ordinal probit mixed models with half-Cauchy priors on the
 * random-effect standard deviations.
 *
 * z_i = x_i' beta + w_i' u + eps_i,  eps_i ~ N(0, 1),  z_i in (l_i, u_i]
 * beta ~ N(mu0, Q0^{-1})                              (known prior)
 * u_k | sigma2_g ~ N(0, sigma2_g),   g = g(k)          (G variance blocks)
 * sigma_g ~ C+(0, s)  <=>  sigma2_g | a_g ~ IG(1/2, 1/a_g),  a_g ~ IG(1/2, 1/s^2)
 *
 * Variational family: q(theta | z) prod_i q(z_i) prod_g q(sigma2_g) q(a_g),
 * with theta = (beta, u). For fixed q(sigma2), the optimal q(theta, z) is the
 * PMF solution under the Gaussian prior N(m0, Q^{-1}), m0 = (mu0, 0) and
 * Q = blockdiag(Q0, diag(tau_{g(k)})), tau_g = E[1/sigma2_g].
 *
 * Each outer iteration:
 *   1. rebuild V = (W'W + Q)^{-1} for the current tau   [update q(theta | z)]
 *   2. one sequential PMF sweep over z_i                [update q(z_i)]
 *   3. evaluate the ELBO and check convergence
 *   4. update q(sigma2_g) = IG((1 + q_g)/2, E[1/a_g] + E||u_g||^2 / 2)
 *      and q(a_g) = IG(1, E[1/sigma2_g] + 1/s^2)
 * The loop exits after step 3, so all returned quantities refer to the same
 * variational state as the reported ELBO.
 */

#include <RcppArmadillo.h>
#include "truncnorm.h"
#include "utils.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

namespace {

// E_q[log IG(sigma2; a, b)] given E[log sigma2] and E[1/sigma2]
double expected_log_ig_hc(const double a, const double b,
                          const double E_log, const double E_inv) {
	return a * std::log(b) - R::lgammafn(a) - (a + 1.0) * E_log - b * E_inv;
}

} // namespace

// [[Rcpp::export]]
Rcpp::List pmf_ordinal_mixed(
		const arma::vec& Y,        // response, coded 1:K
		const arma::mat& X,        // fixed-effect design
		const arma::vec& alpha,    // thresholds (-Inf, ..., Inf)
		const arma::vec& mu0,      // fixed-effect prior mean
		const arma::mat& Q0,       // fixed-effect prior precision
		const arma::mat& Z,        // random-effect design
		const arma::uvec& Z_group, // zero-based variance block per Z column
		const double s_sigma,      // half-Cauchy scale
		const int maxit   = 100,
		const double tresh = 1e-6,
		const int min_iter = 1, // iterations before convergence may be declared
		const std::string conv_crit = "elbo", // "elbo" or "coef"
		const bool verbose = false,
		const bool full_out = false,
		Rcpp::Nullable<Rcpp::List> init = R_NilValue // warm start
		)
{
	if(conv_crit != "elbo" && conv_crit != "coef") {
		Rcpp::stop("conv_crit must be either \"elbo\" or \"coef\".");
	}
	const bool crit_coef = (conv_crit == "coef");

	const int n = X.n_rows;
	const int p = X.n_cols;
	const int q = Z.n_cols;
	const int d = p + q;

	if(Z.n_rows != static_cast<unsigned int>(n)) {
		Rcpp::stop("X and Z must have the same number of rows.");
	}
	if(Z_group.n_elem != static_cast<unsigned int>(q)) {
		Rcpp::stop("Z_group must have length ncol(Z).");
	}
	if(!(s_sigma > 0.0) || !R_FINITE(s_sigma)) {
		Rcpp::stop("s_sigma must be a strictly positive finite number.");
	}

	const int G = arma::max(Z_group) + 1;
	arma::vec q_g(G, fill::zeros);
	for(int k = 0; k < q; k++) q_g(Z_group(k)) += 1.0;

	// Threshold bounds
	arma::vec u(n), l(n);
	arma::ivec y_id = conv_to<ivec>::from(Y);
	for(int i = 0; i < n; i++) {
		u(i) = alpha(y_id(i));
		l(i) = alpha(y_id(i) - 1);
	}

	// Fixed pre-computations
	const arma::mat W   = arma::join_rows(X, Z);
	const arma::mat WtW = W.t() * W;
	const arma::vec Xmu0 = X * mu0;
	// Q m0 = (Q0 mu0, 0) does not depend on tau
	arma::vec r0(d, fill::zeros);
	r0.head(p) = Q0 * mu0;

	double logdetQ0 = 0.0;
	if(!arma::log_det_sympd(logdetQ0, Q0)) {
		Rcpp::stop("Q0 must be symmetric positive definite.");
	}

	// Variance-component factors: shapes are fixed
	const double s2 = s_sigma * s_sigma;
	const arma::vec A_sig = 0.5 * (1.0 + q_g);
	const double A_a = 1.0;
	// Initialise at E[1/sigma2_g] = 1/s^2, E[1/a_g] = 1/(2 s^2)
	arma::vec tau_u(G, fill::value(1.0 / s2));
	arma::vec B_sig = A_sig / tau_u;
	arma::vec B_a   = tau_u + 1.0 / s2;
	arma::vec E_inv_a = A_a / B_a;

	// Working quantities
	arma::mat V(d, d), WV(n, d);
	arma::vec sigma2Z(n), sigmaZ(n), WVr0(n);
	double logdetV = 0.0;

	// PMF sequential state
	arma::vec meanZ(n, fill::zeros);
	arma::vec sdZ(n, fill::zeros);
	arma::vec xiZ(n, fill::zeros);
	arma::vec lpZ(n, fill::zeros); // log(Phi(U_i) - Phi(L_i))
	arma::mat D(n, d, fill::zeros);

	arma::vec elbo_seq(maxit, fill::zeros);
	arma::mat tau_u_seq(maxit, G, fill::zeros);
	int it = 0;
	bool conv = false;
	// state monitored by the "coef" criterion: (m_theta, log tau_u)
	arma::vec state(d + G, fill::zeros), state_old;

	// --- Warm start: resume q(z) and q(sigma2), q(a) from a previous fit ---
	if(init.isNotNull()) {
		Rcpp::List ini(init);
		if(ini.containsElementNamed("meanZ")) {
			arma::vec m_in = Rcpp::as<arma::vec>(ini["meanZ"]);
			if(m_in.n_elem == static_cast<unsigned int>(n)) {
				meanZ = m_in;
				// D_i = sum_{j != i} w_j E[z_j]; the sweep recurses from row n-1
				D.row(n - 1) = (W.t() * meanZ).t() - W.row(n - 1) * meanZ(n - 1);
			}
		}
		if(ini.containsElementNamed("sigma_u2_b")) {
			arma::vec b_in = Rcpp::as<arma::vec>(ini["sigma_u2_b"]);
			if(b_in.n_elem == static_cast<unsigned int>(G) && arma::all(b_in > 0.0)) {
				B_sig = b_in;
				tau_u = A_sig / B_sig;
			}
		}
		if(ini.containsElementNamed("a_u_b")) {
			arma::vec ba_in = Rcpp::as<arma::vec>(ini["a_u_b"]);
			if(ba_in.n_elem == static_cast<unsigned int>(G) && arma::all(ba_in > 0.0)) {
				B_a = ba_in;
				E_inv_a = A_a / B_a;
			}
		}
	}

	while(true) {
		Rcpp::checkUserInterrupt();

		// --- q(theta | z): V = (W'W + Q)^{-1} for the current tau ---
		arma::mat prec = WtW;
		prec.submat(0, 0, p - 1, p - 1) += Q0;
		for(int k = 0; k < q; k++) prec(p + k, p + k) += tau_u(Z_group(k));
		if(!arma::inv_sympd(V, prec)) {
			arma::inv(V, prec, arma::inv_opts::allow_approx);
		}
		double ldp = 0.0, sgn = 0.0;
		arma::log_det(ldp, sgn, prec);
		logdetV = -ldp;

		WV      = W * V;
		sigma2Z = 1.0 / (1.0 - arma::sum(WV % W, 1));
		sigmaZ  = arma::sqrt(sigma2Z);
		WVr0    = WV * r0;

		// --- q(z_i): sequential PMF sweep ---
		// xi_i = sigma2_i * w_i' V (sum_{j != i} w_j E[z_j] + Q m0)
		int im1;
		for(int i = 0; i < n; i++) {
			im1 = (i == 0) ? (n - 1) : (i - 1);
			D.row(i) = D.row(im1) - W.row(i) * meanZ(i) + W.row(im1) * meanZ(im1);
			xiZ(i) = sigma2Z(i) * (arma::dot(WV.row(i), D.row(i)) + WVr0(i));

			const double Li = (l(i) - xiZ(i)) / sigmaZ(i);
			const double Ui = (u(i) - xiZ(i)) / sigmaZ(i);
			meanZ(i) = xiZ(i) - sigmaZ(i) * zeta1(Li, Ui);
			sdZ(i)   = sigmaZ(i) * std::sqrt(1.0 - z1p2_z2(Li, Ui));
			lpZ(i)   = trunc_log(R::pnorm(Ui, 0.0, 1.0, TRUE, FALSE)
			                   - R::pnorm(Li, 0.0, 1.0, TRUE, FALSE));
		}
		const arma::vec varZ = sdZ % sdZ;

		// --- ELBO ---
		// PMF part: E_q(z)[log N(z; X mu0, I + W Q^{-1} W')] + H[q(z)], using
		// Sigma_z^{-1} = I - W V W' and log|Sigma_z| = -log|V| - log|Q|.
		const arma::vec zc = meanZ - Xmu0;
		const arma::vec rz = W.t() * zc;
		double elbo = 0.5 * logdetV + 0.5 * logdetQ0;
		elbo -= 0.5 * (arma::dot(zc, zc) - arma::as_scalar(rz.t() * V * rz)
		               + arma::accu(varZ / sigma2Z));
		elbo += arma::accu(lpZ) + 0.5 * arma::accu(arma::log(sigma2Z))
		      + 0.5 * arma::accu((varZ + arma::square(meanZ - xiZ)) / sigma2Z);
		// Variance components. The 0.5 * q_g * log(tau_g) in log|Q| cancels
		// with E[log p(u_g | sigma2_g)] - log N(u_g; 0, tau_g^{-1} I).
		for(int g = 0; g < G; g++) {
			const double E_log_s = std::log(B_sig(g)) - R::digamma(A_sig(g));
			const double E_inv_s = A_sig(g) / B_sig(g);
			const double E_log_a = std::log(B_a(g)) - R::digamma(A_a);
			elbo += -0.5 * q_g(g) * E_log_s;
			// E log p(sigma2 | a) + E log p(a)
			elbo += -0.5 * E_log_a - R::lgammafn(0.5) - 1.5 * E_log_s - E_inv_a(g) * E_inv_s;
			elbo += -0.5 * std::log(s2) - R::lgammafn(0.5) - 1.5 * E_log_a - E_inv_a(g) / s2;
			// - E log q(sigma2) - E log q(a)
			elbo -= expected_log_ig_hc(A_sig(g), B_sig(g), E_log_s, E_inv_s);
			elbo -= expected_log_ig_hc(A_a, B_a(g), E_log_a, E_inv_a(g));
		}

		elbo_seq(it) = elbo;
		tau_u_seq.row(it) = tau_u.t();

		// Marginal moments of theta: Var(theta) = V + (WV)' diag(varZ) WV
		const arma::vec m_theta = WV.t() * meanZ + V * r0;
		const arma::vec v_theta = V.diag() + (WV % WV).t() * varZ;

		double rel_coef = arma::datum::inf;
		state.head(d) = m_theta;
		state.tail(G) = arma::log(tau_u);
		if(it >= min_iter) {
			if(crit_coef) {
				// monitor the variational parameters themselves: the ELBO can be
				// nearly flat while the variance components still move
				rel_coef = max_rel_change(state, state_old);
				conv = (rel_coef < tresh);
			} else {
				// relative tolerance: the ELBO scales with n
				conv = (std::abs(elbo_seq(it) - elbo_seq(it - 1)) <
				        tresh * (1.0 + std::abs(elbo_seq(it))));
			}
		}
		state_old = state;
		it++;

		if(verbose && (it % 10 == 0 || conv)) {
			Rcout << "it: " << it << "  elbo: " << elbo;
			if(crit_coef) Rcout << "  max rel. change: " << rel_coef;
			Rcout << std::endl;
		}
		if(conv || it >= maxit) break;

		// --- q(sigma2_g) and q(a_g) ---
		arma::vec E_u_ss(G, fill::zeros);
		for(int k = 0; k < q; k++) {
			const double mk = m_theta(p + k);
			E_u_ss(Z_group(k)) += mk * mk + v_theta(p + k);
		}
		B_sig   = E_inv_a + 0.5 * E_u_ss;
		tau_u   = A_sig / B_sig;
		B_a     = tau_u + 1.0 / s2;
		E_inv_a = A_a / B_a;
	}

	// --- Marginal moments of theta = (beta, u) under q ---
	const arma::vec m_theta = WV.t() * meanZ + V * r0;
	arma::mat WV_sd = WV;
	WV_sd.each_col() %= sdZ;
	const arma::mat S_theta = V + WV_sd.t() * WV_sd;

	arma::vec sigma_u2_mean(G);
	for(int g = 0; g < G; g++) {
		sigma_u2_mean(g) = (A_sig(g) > 1.0) ? B_sig(g) / (A_sig(g) - 1.0) : R_PosInf;
	}

	Rcpp::List out;
	out["m"]       = m_theta.head(p);
	out["S"]       = S_theta.submat(0, 0, p - 1, p - 1);
	out["m_u"]     = m_theta.tail(q);
	out["S_u"]     = S_theta.submat(p, p, d - 1, d - 1);
	out["m_joint"] = m_theta;
	out["S_joint"] = S_theta;
	out["it"]      = it;
	out["sigma_u2_a"]        = A_sig;
	out["sigma_u2_b"]        = B_sig;
	out["sigma_u2_mean"]     = sigma_u2_mean;
	out["sigma_u2_inv_mean"] = tau_u; // tau used to build V (and xiZ, sigmaZ)
	out["a_u_a"]             = arma::vec(G, fill::value(A_a));
	out["a_u_b"]             = B_a;
	out["u_group"]           = arma::regspace<arma::uvec>(0, G - 1);
	out["elbo"]              = elbo_seq(it - 1);
	out["conv"]              = conv;

	out["meanZ"] = meanZ; // for warm starting a subsequent fit

	if(full_out) {
		out["sigmaZ"]   = sigmaZ;
		out["xiZ"]      = xiZ;
		out["elbo_seq"] = elbo_seq.subvec(0, it - 1);
		out["sigma_u2_inv_mean_seq"] = tau_u_seq.rows(0, it - 1);
	}

	return out;
}
