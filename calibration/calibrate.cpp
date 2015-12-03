// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;


// linear calibration;

// [[Rcpp::export]]
colvec calib_linear (mat X,
                     colvec d,
                     colvec totals,
                     double eps=1e-06) {
  mat Xd = X;
  Xd.each_col() %=d;
  colvec lambda  = pinv(Xd.t() * X, eps) * (totals - (d.t() * X).t());
  colvec w = (1 + X * lambda) % d;
  return (w);
}

// arcsin calibration;

colvec dasinh (colvec u, double alpha) {
  colvec p = asinh(2 * alpha * u);
  colvec d = (p + sqrt((pow(p,2)) + 4*( pow(alpha,2))))/(2*alpha);
  return (d);
}

// [[Rcpp::export]]
colvec calib_asinh(mat X,
                   colvec d,
                   colvec totals,
                   rowvec bounds,
                   double eps=1e-06,
                   double tol = 1e-06,
                   int maxit = 100,
                   double alpha = 1,
                   bool verbose = false) {
  int n = X.n_rows;
  int k = X.n_cols;
  mat lambda(k,1);
  colvec w = d;
  int iter = 0;
  
  while (iter <= maxit &  
         max(abs(X.t() * w - totals)/totals) >= tol &
         !w.has_nan() & 
         !w.has_inf()) {
         
         w = d % dasinh(X * lambda,alpha);
    
    uvec b1 = find(w/d < bounds(0));
    uvec b2 = find(w/d > bounds(1));
    
    if (!b1.is_empty()){
      w.elem(b1) = bounds(1) * d.elem(b1);
    }
    if (!b2.is_empty()) {
      w.elem(b2) = bounds(2) * d.elem(b2);
    }
    
    mat phi = X.t() * w - totals;
    mat dphi = X;
    dphi.each_col() %=d;
    dphi = dphi.t() * X;
    lambda = lambda - pinv(dphi,eps) * phi;
    
    if (verbose) {
      Rcout << "Iteration" << std::endl << iter << std::endl;
    }
    iter = iter + 1;
    
  }
  return (w);
}

// // block calibration
// SEXP block_calibration(SEXP X, sexp d, sexp)
// 
// arguments:
//' @X - array of X values
//' @d - array of weights
//' @totals -- array of totals
// SEXP arma_cube (SEXP array) {
//   NumericVector vecArray(array);
//   IntegerVector arrayDims = vecArray.attr("dim");
//   cube cubeArray(vecArray.begin(), 
//                  arrayDims[0], arrayDims[1], arrayDims[2], false);
//   return(wrap(cubeArray));
// }



