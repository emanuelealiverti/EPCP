/*
 * VB-PMF for ordinal probit regression with IG prior on variance.
 *
 * beta | sigma_b2 ~ N(mu0, sigma_b2 I_p)
 * sigma_b2 ~ IG(a0, b0)
 *
 * The PMF sequential sweep is repeated each outer iteration; V, XV, sigma2Z
 * and sigma_xv_prior are re-derived from the current tau_b = E[1/sigma_b2]
 * after each sweep.
 */

#include <RcppArmadillo.h>
#include "truncnorm.h"
#include "utils.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List pmf_ordinal_prior(
        const arma::vec& Y,
        const arma::mat& X,
        const arma::vec& alpha,
        const arma::vec& mu0,
        const double a0,
        const double b0,
        const int maxit   = 100,
        const double tresh = 1e-6,
        const bool verbose = false,
        const bool full_out = false
)
{
    const int n = X.n_rows;
    const int p = X.n_cols;

    // Threshold bounds
    arma::vec u(n), l(n);
    arma::ivec y_id = conv_to<ivec>::from(Y);
    for(int i = 0; i < n; i++) {
        u(i) = alpha(y_id(i));
        l(i) = alpha(y_id(i) - 1);
    }

    // IG posterior shape — fixed at a0 + p/2
    const double a_sigma_b2 = a0 + 0.5 * p;
    double b_sigma_b2 = b0;
    double tau_b      = a0 / b0;   // initial E[1/sigma_b2]

    // Fixed pre-computations
    const arma::mat XtX    = X.t() * X;
    const arma::vec Xmu0   = X * mu0;
    const arma::vec XtXmu0 = XtX * mu0;  // p-vector, used for sigma_xv_prior

    // Working quantities updated each outer iteration
    arma::mat V(p, p);
    arma::mat XV(n, p);
    arma::vec sigma2Z(n), sigmaZ(n), sigma_xv_prior(n);

    // PMF sequential state
    arma::vec meanZ(n, fill::zeros);
    arma::vec sdZ(n, fill::zeros);
    arma::vec xiZ(n, fill::zeros);
    arma::mat D(n, p, fill::zeros);

    arma::vec elbo_seq(maxit, fill::zeros);
    double lp = 0.0;
    int it = 0;
    bool conv = false;

    while(!conv && it < maxit) {
        Rcpp::checkUserInterrupt();

        // --- Re-derive V, XV, sigma2Z, sigma_xv_prior for current tau_b ---
        arma::mat prec = XtX;
        prec.diag() += tau_b;
        arma::inv(V, prec, arma::inv_opts::allow_approx);
        XV      = X * V;
        sigma2Z = 1.0 / (1.0 - arma::sum(XV % X, 1));
        sigmaZ  = arma::sqrt(sigma2Z);
        // sigma_xv_prior(i) = sigma2Z(i) * XV_i * (XtX - x_i x_i') * mu0
        // equivalent to:  sigma2Z % (XV * XtXmu0 - Xmu0) + Xmu0
        sigma_xv_prior = sigma2Z % (XV * XtXmu0 - Xmu0) + Xmu0;

        // --- PMF sequential sweep ---
        lp = 0.0;
        int im1;
        for(int i = 0; i < n; i++) {
            im1 = (i == 0) ? (n - 1) : (i - 1);
            D.row(i) = D.row(im1) - X.row(i) * meanZ(i) + X.row(im1) * meanZ(im1);
            xiZ(i)   = Xmu0(i) - sigma_xv_prior(i)
                       + sigma2Z(i) * arma::as_scalar(XV.row(i) * D.row(i).t());

            const double Li = (l(i) - xiZ(i)) / sigmaZ(i);
            const double Ui = (u(i) - xiZ(i)) / sigmaZ(i);
            meanZ(i) = xiZ(i) - sigmaZ(i) * zeta1(Li, Ui);
            sdZ(i)   = sigmaZ(i) * std::sqrt(1.0 - z1p2_z2(Li, Ui));

            lp += trunc_log(R::pnorm(Ui, 0.0, 1.0, TRUE, FALSE)
                          - R::pnorm(Li, 0.0, 1.0, TRUE, FALSE));
        }

        // --- Marginal posterior moments of beta ---
        const arma::vec meanBeta = XV.t() * meanZ + tau_b * (V * mu0);

        // --- Update q(sigma_b2) ---
        // E[||beta - mu0||^2] = ||m - mu0||^2 + tr(V) + sum_i sdZ_i^2 ||XV_i||^2
        const arma::vec beta_shift = meanBeta - mu0;
        const double E_beta_ss =
            arma::dot(beta_shift, beta_shift)
            + arma::trace(V)
            + arma::dot(sdZ % sdZ, arma::sum(XV % XV, 1));

        const double tau_b_old = tau_b;
        b_sigma_b2 = b0 + 0.5 * E_beta_ss;
        tau_b      = a_sigma_b2 / b_sigma_b2;

        // ELBO proxy via elbo_pmf (sufficient for outer convergence monitoring)
        elbo_seq(it) = elbo_pmf(X, V, XV, xiZ, sigma2Z, meanZ, Xmu0, mu0, lp);

        if(it > 0) {
            conv = (std::abs(tau_b - tau_b_old) < tresh);
        }
        it++;

        if(verbose && (it % 10 == 0 || conv)) {
            Rcout << "it: "    << it
                  << "  tau_b: " << tau_b
                  << "  elbo: "  << elbo_seq(it - 1)
                  << std::endl;
        }
    }

    // --- Final marginal posterior of beta ---
    const arma::vec meanBeta = XV.t() * meanZ + tau_b * (V * mu0);
    arma::mat XV_sd = XV;
    XV_sd.each_col() %= sdZ;
    const arma::mat varBeta = V + XV_sd.t() * XV_sd;

    const double sigma_b2_mean =
        (a_sigma_b2 > 1.0) ? b_sigma_b2 / (a_sigma_b2 - 1.0) : R_PosInf;

    Rcpp::List out;
    out["m"]                = meanBeta;
    out["S"]                = varBeta;
    out["it"]               = it;
    out["sigma_b2_a"]       = a_sigma_b2;
    out["sigma_b2_b"]       = b_sigma_b2;
    out["sigma_b2_mean"]    = sigma_b2_mean;
    out["sigma_b2_inv_mean"] = tau_b;

    if(full_out) {
        out["sigmaZ"]   = sigmaZ;
        out["xiZ"]      = xiZ;
        out["elbo_seq"] = elbo_seq.subvec(0, it - 1);
        out["elbo"]     = elbo_seq(it - 1);
        out["conv"]     = conv;
    }

    return out;
}
