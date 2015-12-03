// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;

// based on http://markovjumps.blogspot.it/2011/12/r-array-to-rcpparmadillo-cube.html;

// [[Rcpp::export]]
SEXP arma_cube (SEXP array) {
  NumericVector vecArray(array);
  IntegerVector arrayDims = vecArray.attr("dim");
  cube cubeArray(vecArray.begin(), 
                 arrayDims[0], 
                 arrayDims[1], 
                 arrayDims[2], 
                 false);
  return(wrap(cubeArray));
}


// example code

/*** R
set.seed(345)
testArray  <- array(rnorm(18), dim=c(3,3,2))
arma_cube(testArray)
*/
