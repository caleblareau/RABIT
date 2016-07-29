#include <Rcpp.h>
using namespace Rcpp;

//' @title
//' runRABIT
//' @description
//' Performs regression analysis
//'
//' @param x Input variable matrix. Each column represents a variable to be selected.
//' REQUIRED.
//' @param y Response matrix. Each column represents a response vector and Rabit will
//' run variable selection on each response separately. REQUIRED.
//' @param b Background factors. Each column represents a confounding factor to be controlled.
//' Rabit will exclude these effects from the linear model. Default = NULL
//' @param c Background factors for each column of Y. Different from b which inputs background
//' factors shared across all columns of Y, this input is specific for each column of Y. Rabit
//' will search for shared column names between Y and input background factors here and control
//' these factors. default = NULL
//' @param f FDR threshold. Range (0,1]. Rabit estimates the statistical significance of each
//' variable in X and only keeps significant ones in later stepwise forward regression. Set
//' to 1 to skip this step. default = 0.05
//' @param t Transform Y to normal distribution ( = TRUE) or not ( = FALSE). Transforming response
//' Y to normal distribution increases the statistical power of linear regression t-test. However,
//' if Y is very different from the normal distribution (e.g. binary outcome), then this
//' methodology should not be used and logistic regression should be used in this example.
//' default = TRUE
//' @param s Select one best variable in X for each category ( = TRUE) or not ( = FALSE).
//'  When several  variables in X come from the same category (e.g. Transcription factors
//' may have multiple ChIP-Seq profiles available), multiple names appended together
//' with a . (e.g. CatA.V1, CatA.V2 ...) Rabit will chose the best one based on the
//' maximum t-value and these variableswill be returned. default = FALSE
//' @param r Run forward stepwise selection ( = TRUE) or not ( = FALSE). If set to FALSE,
//' Rabit will only calculate the statistical significance of each variable without
//' forward stepwise selection. default = TRUE.
//'
//' @details
//' \code{runRABIT} as a function call is a one-stop shop. Carefully consider each parameter before
//' proceeding.
//'
//' @return
//' Who knows but hopefully it's okay.
//'
//' @export
// [[Rcpp::export]]
List runRABIT(
        NumericMatrix x,
        NumericMatrix y,
        NumericMatrix b = NumericMatrix(0),
        NumericMatrix c = NumericMatrix(0),
        float f = 0.05,
        bool t = true,
        bool s = false,
        bool r = true) {
  // CharacterVector xx = CharacterVector::create("foo", "bar");
  // NumericVector yy   = NumericVector::create(0.0, 1.0);
  List z            = List::create(x, y);
  return z;
}
