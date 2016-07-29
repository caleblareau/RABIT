#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppGSL.h>

//------------------------------------------------------------------------------
// start gsl_util
//------------------------------------------------------------------------------

/*** R
dyn.load("lm.so")
dyn.load("lin.so")
dyn.load("gsl_util.so")
*/


//------------------------------------------------------------------------------
// end util.cpp
//------------------------------------------------------------------------------

//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


//------------------------------------------------------------------------------
// start util
//------------------------------------------------------------------------------

#include <gsl/gsl_cdf.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace std;
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_blas_types.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas_types.h>
#include <gsl/gsl_blas.h>
#include <math.h>

#define EPS 1e-10
#define SQRT_EPS(x) (fabs(x)<EPS?0:sqrt(x))

template <class T>
class rank_node
{
public:
	rank_node(T v, size_t i):v(v),i(i){}
	bool operator < (const rank_node &r) const {return v < r.v;}

	T v;
	size_t i;
};

// transform values in src to quantiles in dst. If dst == src, transform in place

// [[Rcpp::export]]
void quantile_transform(NumericVector dst, NumericVector src, int N)
{
	size_t i;
	vector<rank_node<double> > sortvec;

	for(i=0;i<N;i++) sortvec.push_back(rank_node<double>(src[i],i));

	sort(sortvec.begin(), sortvec.end());

	for(i=0;i<N;i++) dst[sortvec[i].i] = (double)(i+1)/(N+1);
}

// transform values in src to quantiles in dst. If dst == src, transform in place

// [[Rcpp::export]]
void normal_transform(NumericVector dst, NumericVector src, int N)
{
	size_t i;
	double aver=0, sd=0, v;

	// no point to do normal transformation
	if(N==1) return;

	for (i=0;i<N;i++)
	{
		v = (double)src[i];
		aver += v;
		sd += v*v;
	}

	aver /= N;
	sd = (sd - N*aver*aver)/(N-1);
	sd = (fabs(sd)<EPS?0:sqrt(sd));

	quantile_transform(dst, src, N);
	for(i=0;i<N;i++) dst[i] = gsl_cdf_gaussian_P(dst[i], sd) + aver;
}

// [[Rcpp::export]]
void Benjamini_Hochberg(NumericVector FDR, NumericVector pvalue, int N)
{
	size_t i;
	double qvalue = 1;
	vector<rank_node<double> > sortvec;

	for(i=0;i<N;i++) sortvec.push_back(rank_node<double>(pvalue[i],i));

	sort(sortvec.begin(), sortvec.end());

	for(i=0;i<N;i++)
	{
		qvalue = min<double>(qvalue, sortvec[N-1-i].v*N/(N-i));
		FDR[sortvec[N-1-i].i] = qvalue;
	}
}

// [[Rcpp::export]]
int triangle_index(int i, int j)
{
	if(i==j){
		cerr << "Error: cannot use i==j in triangle index." << endl;
		exit(1);

	}else if(i<j) swap<int>(i,j);

	return ((i*(i-1))>>1) + j;
}

//------------------------------------------------------------------------------
// end util.cpp
//------------------------------------------------------------------------------

//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


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
//' @import RcppGSL
//' @import RcppArmadillo
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
  NumericVector yy   = NumericVector::create(0.0, 1.0, 3.2, 5.5, 4.4);
  NumericVector y2   =  NumericVector::create(0.0, 1.0, 3.2, 5.5, 4.4);
  NumericVector qq   = NumericVector::create(0.1, 1.1, 3.22, 1.5, 2.4);
  Benjamini_Hochberg(yy, qq, 5);
  List z            = List::create(x, yy, triangle_index(4,5));
  return z;
}
