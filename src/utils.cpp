  // [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>

static double const log2pi = std::log(2.0 * M_PI);

/* C++ version of the dtrmv BLAS function */
void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
  arma::uword const n = trimat.n_cols;

  for(unsigned j = n; j-- > 0;){
    double tmp(0.);
    for(unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
    x[j] = tmp;
  }
}

// [[Rcpp::export]]
arma::vec dmvnrm_arma_fast(arma::mat const &x,
                           arma::rowvec const &mean,
                           arma::mat const &sigma,
                           bool const logd = false) {
  using arma::uword;
  uword const n = x.n_rows,
    xdim = x.n_cols;
  arma::vec out(n);
  arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
  double const rootisum = arma::sum(log(rooti.diag())),
    constants = -(double)xdim/2.0 * log2pi,
    other_terms = rootisum + constants;

  arma::rowvec z;
  for (uword i = 0; i < n; i++) {
    z = (x.row(i) - mean);
    inplace_tri_mat_mult(z, rooti);
    out(i) = other_terms - 0.5 * arma::dot(z, z);
  }

  if (logd)
    return out;
  return exp(out);
}

// [[Rcpp::export]]
arma::vec cpp_eval_ldnorm(arma::mat x, arma::mat mu, arma::cube sigma, arma::uvec z) {
  int n = x.n_rows;

  arma::vec logp(n);

  for(int i = 0; i < n; i++) {
    arma::mat mu_g = mu.row(z(i)-1);
    arma::mat sigma_g = sigma.slice(z(i)-1);
    logp(i) = dmvnrm_arma_fast(x.row(i), mu_g, sigma_g, true)(0);
  }

  return logp;
}



// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  // Ensure symmetry
  arma::mat sigma_sym = 0.5 * (sigma + sigma.t());
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma_sym);
}

// // [[Rcpp::export]]
// arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
//   int ncols = sigma.n_cols;
//   arma::mat Y = arma::randn(n, ncols);
//
//   // Force symmetry
//   arma::mat sigma_sym = 0.5 * (sigma + sigma.t());
//
//   // Add small ridge if not positive definite
//   arma::mat L;
//   bool success = false;
//   double jitter = 1e-8;
//   int max_tries = 5;
//   int tries = 0;
//
//   while (!success && tries < max_tries) {
//     try {
//       L = arma::chol(sigma_sym);
//       success = true;
//     } catch (...) {
//       sigma_sym += jitter * arma::eye(sigma_sym.n_rows, sigma_sym.n_cols);
//       jitter *= 10;  // increase jitter each time
//       tries++;
//     }
//   }
//
//   if (!success) {
//     Rcpp::stop("Cholesky decomposition failed even after ridge regularization.");
//   }
//
//   return arma::repmat(mu, 1, n).t() + Y * L;
// }


// [[Rcpp::export]]
arma::mat cpp_compute_V(arma::mat X,
                        arma::vec omega,
                        arma::mat inv_cov) {

  arma::mat XtOmegaX = X.each_col() % omega;
  XtOmegaX = XtOmegaX.t() * X;
  XtOmegaX += 1e-7 * arma::eye(X.n_cols, X.n_cols);
  arma::mat V = arma::inv(XtOmegaX + inv_cov);

  return V;
}

// [[Rcpp::export]]
arma::mat cpp_compute_m(arma::mat V,
                        arma::mat X,
                        arma::vec z,
                        arma::vec omega,
                        arma::vec C,
                        arma::vec center,
                        arma::mat inv_cov) {
  arma::mat m = V * (X.t() * ((z - 0.5) + (omega % C)) + inv_cov * center);
  return m;
}


// [[Rcpp::export]]
arma::mat sample_alpha(arma::mat X,
                       arma::vec omega,
                       arma::mat inv_cov,
                       arma::vec z,
                       arma::vec C,
                       arma::vec center) {

  arma::mat V = cpp_compute_V(X, omega, inv_cov);
  arma::mat mu = cpp_compute_m(V, X, z, omega, C, center, inv_cov);

  return mvrnormArma(1, mu, V).t();

}


// [[Rcpp::export]]
double logsumexp_cpp(const arma::vec& x) {
  double max_val = x.max();
  return max_val + std::log(arma::sum(arma::exp(x - max_val)));
}
