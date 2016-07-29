#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::depends(RcppArmadillo)]]


//------------------------------------------------------------------------------
// start util.h && util.cpp
//------------------------------------------------------------------------------

#include <gsl/gsl_cdf.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace std;

#define	EPS	1e-10
#define	DISPLAY_BOOL(x) ((x)?("1 (yes)"):("0 (no)"))


// intersection of two sorted vectors. Assume vector is ordered.
template <class T>
void intersect(const vector<T> &s1, const vector<T> &s2, vector<T> &result)
{
	result.clear();

	typename vector<T>::const_iterator
		first1 = s1.begin(), last1 = s1.end(),
		first2 = s2.begin(), last2 = s2.end();

	while (first1 != last1 && first2 != last2)
	{
		if (*first1 < *first2)
			first1++;
		else if (*first2 < *first1)
			first2++;
		else {
			result.push_back(*first1);
			first1++;
			first2++;
		}
	}
}


// print vector elements
template <class T>
void print_vector(const vector<T> &v, const string &title="")
{
	typename vector<T>::const_iterator iter;

	cout << title;
	for(iter=v.begin(); iter!=v.end(); iter++) cout << '\t' << *iter;
	cout << endl;
}


// rank node structure for rank transform function
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
template <class T>
void quantile_transform(double dst[], T src[], const size_t N)
{
	size_t i;
	vector<rank_node<T> > sortvec;

	for(i=0;i<N;i++) sortvec.push_back(rank_node<T>(src[i],i));

	sort(sortvec.begin(), sortvec.end());

	for(i=0;i<N;i++) dst[sortvec[i].i] = (double)(i+1)/(N+1);
}


// transform values in src to quantiles in dst. If dst == src, transform in place
template <class T>
void normal_transform(double dst[], T src[], const size_t N)
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

	for(i=0;i<N;i++) dst[i] = gsl_cdf_gaussian_Pinv(dst[i], sd) + aver;
}

size_t triangle_index(size_t i, size_t j);

void Benjamini_Hochberg(double FDR[], double pvalue[], const size_t N);

void Benjamini_Hochberg(double FDR[], double pvalue[], const size_t N)
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

size_t triangle_index(size_t i, size_t j)
{
	if(i==j){
		cerr << "Error: cannot use i==j in triangle index." << endl;
		exit(1);

	}else if(i<j) swap<size_t>(i,j);

	return ((i*(i-1))>>1) + j;
}

//------------------------------------------------------------------------------
// end util.h & util.cpp
//------------------------------------------------------------------------------

//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

//------------------------------------------------------------------------------
// start gsl_util.h && gsl_util.c
//------------------------------------------------------------------------------

#include <gsl/gsl_blas.h>

// EPS correction for sqrt function on nearly 0 value
#define SQRT_EPS(x) (fabs(x)<EPS?0:sqrt(x))

// return 1 if v variation is 0 or size is smaller than 2
// if standardize, make norm2 = 1
int Z_normalize(gsl_vector *v, const int standardize);


// change the dimension of matrices and vectos
void resize_matrix(gsl_matrix *m, const size_t size1, const size_t size2, double *data);
void resize_vector(gsl_vector *v, const size_t size, double *data);

void print_vector(const gsl_vector *v, const char *title, FILE *fp);
void print_matrix(const gsl_matrix *X, const char *title, FILE *fp);


// check the dimension of matrices and vectos
void check_matrix_dimension(const gsl_matrix *X, const size_t n, const size_t p,
		const char *title, const int equal_stride);

void check_vector_dimension(const gsl_vector *X, const size_t n,
		const char *title, const size_t stride);


//------------------------------------------------------------------------------
// end gsl_util.h BUT NEED TO ADD gsl_util.c
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
  // NumericVector yy   = NumericVector::create(0.0, 1.0);
  List z            = List::create(x, y);
  return z;
}
