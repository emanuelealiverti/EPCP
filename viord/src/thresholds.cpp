/*
 * Newton-Raphson step for the thresholds of the cumulative probit model.
 *
 * The only term of the objective (ELBO, or EP log marginal likelihood) that
 * depends on the thresholds is
 *
 *   l(alpha) = sum_i log[ Phi((alpha_{y_i} - xi_i)/sigma_i)
 *                       - Phi((alpha_{y_i-1} - xi_i)/sigma_i) ],
 *
 * with alpha_0 = -Inf and alpha_K = +Inf, where (xi_i, sigma_i) are the
 * location and scale of q(z_i): under MF the linear predictor and 1, under PMF
 * xiZ and sigmaZ.
 *
 * l is concave in alpha for the probit link, and each observation only involves
 * the two thresholds bracketing its category, so the Hessian is tridiagonal and
 * a Newton step costs O(n + K). Ordering is preserved by backtracking: the step
 * is halved until the thresholds are still increasing and the objective has not
 * decreased, which also guards against the rare non-concave excursion caused by
 * rounding in the tails.
 */

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

namespace {

const double MIN_P = 1e-300;

// Phi(a) - Phi(b) evaluated in whichever tail keeps the difference accurate
inline double norm_diff(const double a, const double b) {
	if(!R_FINITE(b) && b < 0) return R::pnorm(a, 0.0, 1.0, TRUE, FALSE);
	if(!R_FINITE(a) && a > 0) return R::pnorm(b, 0.0, 1.0, FALSE, FALSE);
	if(b > 0.0) {
		// both in the upper tail: Phi(-b) - Phi(-a) loses no precision
		return R::pnorm(b, 0.0, 1.0, FALSE, FALSE) - R::pnorm(a, 0.0, 1.0, FALSE, FALSE);
	}
	return R::pnorm(a, 0.0, 1.0, TRUE, FALSE) - R::pnorm(b, 0.0, 1.0, TRUE, FALSE);
}

inline double dens(const double x) {
	return R_FINITE(x) ? R::dnorm(x, 0.0, 1.0, FALSE) : 0.0;
}

double loglik(const arma::ivec& y, const arma::vec& xi, const arma::vec& sigma,
              const arma::vec& alpha, const int K) {
	double out = 0.0;
	for(arma::uword i = 0; i < y.n_elem; i++) {
		const int j = y(i);
		const double a = (j == K) ? R_PosInf : (alpha(j - 1) - xi(i)) / sigma(i);
		const double b = (j == 1) ? R_NegInf : (alpha(j - 2) - xi(i)) / sigma(i);
		out += std::log(std::max(norm_diff(a, b), MIN_P));
	}
	return out;
}

} // namespace

// [[Rcpp::export]]
Rcpp::List newton_thresholds(
		const arma::ivec& y,     // response, coded 1:K
		const arma::vec& xi,     // location of q(z_i)
		const arma::vec& sigma,  // scale of q(z_i)
		const arma::vec& start,  // starting thresholds, length K - 1
		const int maxit = 50,
		const double tol = 1e-10
		)
{
	const arma::uword n = y.n_elem;
	if(xi.n_elem != n || sigma.n_elem != n) {
		Rcpp::stop("y, xi and sigma must have the same length.");
	}
	const int K = start.n_elem + 1;
	if(K < 2) Rcpp::stop("at least two categories are required.");

	arma::vec alpha = start;
	arma::vec g(K - 1);
	arma::mat H(K - 1, K - 1); // tridiagonal
	double ll = loglik(y, xi, sigma, alpha, K);
	int it = 0;
	bool conv = false;

	while(it < maxit) {
		Rcpp::checkUserInterrupt();
		g.zeros();
		H.zeros();

		for(arma::uword i = 0; i < n; i++) {
			const int j = y(i);
			const double s = sigma(i);
			const double a = (j == K) ? R_PosInf : (alpha(j - 1) - xi(i)) / s;
			const double b = (j == 1) ? R_NegInf : (alpha(j - 2) - xi(i)) / s;
			const double P = std::max(norm_diff(a, b), MIN_P);
			const double fa = dens(a), fb = dens(b);
			const double ra = fa / (s * P), rb = fb / (s * P); // d/d alpha of log P

			// gradient: upper threshold of the category gets +ra, lower one -rb
			if(j < K)  g(j - 1) += ra;
			if(j > 1)  g(j - 2) -= rb;

			// Hessian, tridiagonal: only the two thresholds of the category
			if(j < K) {
				const double aa = R_FINITE(a) ? a : 0.0;
				H(j - 1, j - 1) += (-aa * fa / (s * s * P)) - ra * ra;
			}
			if(j > 1) {
				const double bb = R_FINITE(b) ? b : 0.0;
				H(j - 2, j - 2) += (bb * fb / (s * s * P)) - rb * rb;
			}
			if(j > 1 && j < K) {
				const double cross = fa * fb / (s * s * P * P);
				H(j - 1, j - 2) += cross;
				H(j - 2, j - 1) += cross;
			}
		}

		if(arma::abs(g).max() < tol) { conv = true; break; }

		arma::vec step;
		if(!arma::solve(step, -H, g, arma::solve_opts::no_approx)) {
			// fall back on a gradient step if the Hessian is numerically singular
			step = g / std::max(arma::abs(H.diag()).max(), 1.0);
		}

		// backtracking: keep the thresholds ordered and the objective increasing
		double t = 1.0;
		bool ok = false;
		arma::vec cand;
		for(int ls = 0; ls < 30; ls++) {
			cand = alpha + t * step;
			if(cand.is_finite() && (K == 2 || arma::all(arma::diff(cand) > 0.0))) {
				const double ll_new = loglik(y, xi, sigma, cand, K);
				if(ll_new >= ll) { ll = ll_new; ok = true; break; }
			}
			t *= 0.5;
		}
		if(!ok) { conv = true; break; } // no further progress possible

		const double delta = arma::abs(cand - alpha).max();
		alpha = cand;
		it++;
		if(delta < tol) { conv = true; break; }
	}

	return Rcpp::List::create(
		Rcpp::Named("alpha") = alpha,
		Rcpp::Named("loglik") = ll,
		Rcpp::Named("it") = it,
		Rcpp::Named("conv") = conv);
}
