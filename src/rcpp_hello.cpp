#include <Rcpp.h>
using namespace Rcpp;

//' @title
//' rcpp_hellow()
//' @description
//' Prints basic
//'
//' @param q if we parameterize the function
//'
//' @details
//' \code{rcpp_hello} does basic calling
//'
//' @return
//' Returns a list object
//'
//' @export
// [[Rcpp::export]]
List rcpp_hellow(int q) {
  CharacterVector x = CharacterVector::create("foo", "bar");
  NumericVector y   = NumericVector::create(0.0, 1.0);
  List z            = List::create(x, y);
  return z;
}
