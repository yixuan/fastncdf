#include "fastncdf.h"
#include "Rcpp.h"
#include <stdexcept>

using Rcpp::NumericVector;

//' Piecewise Linear Interpolation for pnorm
//' @description
//' Performs linear interpolation of \code{\link{pnorm}}. See the README
//' at \url{https://github.com/yixuan/fastncdf} for details.
//'
//' The \code{fastpnorm_preallocated} version is faster if one has
//' vector of the same length already.
//'
//' @param q vector of quantiles.
//' @param p pector vector outputs (probabilities). Needs to have the same
//' length as \code{q}.
//' @param precise_tail logical for whether to use \code{\link{pnorm}} in
//' the tail.
//' @export
//'
//' @examples
//' x <- seq(-6, 6, by = 1e-6)
//' system.time(y <- pnorm(x))
//' system.time(fasty <- fastpnorm(x))
//' all.equal(y, fasty) # tiny error
// [[Rcpp::export(rng = false)]]
NumericVector fastpnorm(NumericVector q, bool const precise_tail = false){
  R_len_t const n = q.size();
  NumericVector out(n);
  fastncdf(&q[0], &out[0], n, precise_tail);
  return out;
}

//' @rdname fastpnorm
//' @export
// [[Rcpp::export(rng = false)]]
void fastpnorm_preallocated(NumericVector q, NumericVector p,
                            bool const precise_tail = false){
  R_len_t const n = q.size();
  if(p.size() != n)
    throw std::invalid_argument("'p' and 'q' are not the same size");
  fastncdf(&q[0], &p[0], n, precise_tail);
}
