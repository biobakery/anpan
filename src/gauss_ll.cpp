#include <Rcpp.h>
using namespace Rcpp;

//' @export
// [[Rcpp::export]]
int add(int x, int y, int z) {
  int sum = x + y + z;
  return sum;
}

// [[Rcpp::export]]
double l2p() {
  return log(2*M_PI);
}

// [[Rcpp::export]]
NumericVector llij_gauss(NumericVector m1, NumericVector s1, NumericVector l, NumericVector yj, NumericVector s2) {

  int n = m1.length();

  NumericVector div_sq_s1 = 1 / pow(s1, 2);

  NumericVector div_sq_s2 = 1 / pow(s2, 2);

  NumericVector a = -.5 * (div_sq_s1 + div_sq_s2);

  NumericVector b = (m1 * div_sq_s1) + ((yj - l) * div_sq_s2);

  NumericVector c = 0.5 * (-(pow(m1,2) * div_sq_s1) -
      (pow(l,2) - 2*l*yj + pow(yj,2)) * div_sq_s2 -
      (log(2 * M_PI) + 2*log(s1)) -
      (log(2 * M_PI) + 2*log(s2)));

  // NumericVector res = log(sqrt(M_PI) * exp(c - (pow(b, 2) / (4*a))) / sqrt(-a));
  NumericVector res = log(sqrt(M_PI)) + c - (pow(b, 2) / (4*a)) - (log(-a)/2);

  return res;
}

// [[Rcpp::export]]
List get_s12x22inv(NumericMatrix x) {
  int p = x.ncol();
  int nr = x.nrow();

  List res(p);

  for (int i=0; i<p; ++i) {
    NumericVector v = x(_,i);
    v.attr("dim") = Dimension(1,nr);
    res[i] = v;
  }
  return res;
}
