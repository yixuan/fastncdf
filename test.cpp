#include <Rcpp.h>

// In order to call sourceCpp() in R
// Don't do this elsewhere!
#include "src/fastncdf.c"

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector fastncdf(NumericVector x)
{
    const int n = x.size();
    NumericVector res(n);
    double* px = x.begin();
    double* pr = res.begin();
    for(int i = 0; i < n; i++)
        pr[i] = fastncdf(px[i]);
    return res;
}
