#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppGSL.h>

//------------------------------------------------------------------------------
// import C code
//------------------------------------------------------------------------------

/*** R
dyn.load("lm.so")
dyn.load("lin.so")
dyn.load("gsl_util.so")
*/


//------------------------------------------------------------------------------
// end C code
//------------------------------------------------------------------------------

//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

//------------------------------------------------------------------------------
// import Matrix map
//------------------------------------------------------------------------------


#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdlib>
using namespace std;


double string_to_number(const string &s, const string &errmsg)
{
	double r;
	char *endstr;

	r = strtod(s.c_str(), &endstr);

	if (*endstr!='\0' || endstr == s.c_str())
	{
		cerr << "Error: " << errmsg << ". \"" << s << "\" is not a number." << endl;
		exit(1);
	}

	return r;
}

#include <string>
#include <vector>
#include <map>
using namespace std;

double string_to_number(const string &s, const string &errmsg="");

class Matrix_map {
public:
	Matrix_map();
	~Matrix_map();

	map<string, size_t> row_map, col_map;
	vector<string> rownames, colnames;

	size_t Nrow, Ncol;

	vector<double*> mat;

	// read matrix from file, append 1 to extra space
	void read(const string &file, const size_t append_count = 0);

	// print matrix
	void print() const;
};


Matrix_map::Matrix_map()
:Nrow(0), Ncol(0) {}


Matrix_map::~Matrix_map()
{
	for (vector<double*>::iterator iter = mat.begin(); iter!=mat.begin(); iter++) delete[] *iter;
}

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
// end util.cpp
//------------------------------------------------------------------------------

//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

//------------------------------------------------------------------------------
// start Rabit.cpp auxillary functions
//------------------------------------------------------------------------------

extern "C" {
#include "gsl_util.h"
#include "lm.h"
}

#include "util.h"
#include <sys/time.h>
#include <cstring>
#include <ctime>
#include <set>
#include <fstream>
#include <sstream>
using namespace std;


// verbose output
int verbose_flag = 0;


void get_common_matrix_name(
		const vector<Matrix_map*> &matrix_vec,
		vector<string> &common_names, const bool row_flag)
{
	// align the row names of all matrices
	vector<string> A, B;
	vector<Matrix_map*>::const_iterator miter;

	for(miter=matrix_vec.begin(); miter!=matrix_vec.end(); miter++)
	{
		if(row_flag)
			B = (*miter)->rownames;
		else
			B = (*miter)->colnames;

		sort(B.begin(), B.end());

		if(miter == matrix_vec.begin())
			common_names = B;
		else
			intersect(A, B, common_names);

		A = common_names;
	}
}


gsl_matrix *matrix_map_to_gsl(
		const vector<string> &rownames,
		const vector<string> &colnames,
		const Matrix_map *mat,
		const bool transpose)
{
	size_t i,j,
		Nrow = rownames.size(), Ncol = colnames.size(),
		*arrinx = new size_t[Ncol];

	double *arr;
	gsl_matrix *result;

	if(transpose)
		result = gsl_matrix_alloc(Ncol, Nrow);
	else
		result = gsl_matrix_alloc(Nrow, Ncol);

	map<string, size_t>::const_iterator miter;

	for (i=0;i<Ncol;i++)
	{
		miter = mat->col_map.find(colnames[i]);

		if(miter==mat->col_map.end())
		{	// I assume this never happens
			cerr << "Fatal error: missing element \"" << colnames[i] << "\"." << endl;
			exit(1);
		}

		arrinx[i] = miter->second;
	}


	for(i=0; i<Nrow; i++)
	{
		miter = mat->row_map.find(rownames[i]);

		if(miter==mat->row_map.end())
		{	// I assume this never happens
			cerr << "Fatal error: missing element \"" << rownames[i] << "\"." << endl;
			exit(1);
		}

		arr = mat->mat[miter->second];

		for(j=0;j<Ncol;j++)
		{
			if(transpose)
				gsl_matrix_set(result, j, i, arr[arrinx[j]]);
			else
				gsl_matrix_set(result, i, j, arr[arrinx[j]]);
		}
	}

	return result;
}



void write_gsl_matrix(const gsl_matrix *m,
		const vector<string> &rownames, const vector<string> &colnames,
		ofstream &fout, const bool transpose)
{
	double v;
	size_t i,j, Nrow = rownames.size(), Ncol = colnames.size();

	for(i=0;i<Ncol;i++) fout << colnames[i] << (i==Ncol-1?'\n':'\t');

	for(i=0;i<Nrow;i++)
	{
		fout << rownames[i];

		for(j=0;j<Ncol;j++)
		{
			if (transpose)
				v = gsl_matrix_get(m,j,i);
			else
				v = gsl_matrix_get(m,i,j);

			fout << '\t' << v;
		}

		fout << '\n';
	}
}


void normal_transform_gsl_matrix(gsl_matrix *m)
{
	size_t i, j, n = m->size1, p = m->size2;

	double *arr = new double[p];

	for(i=0;i<n;i++)
	{
		gsl_vector_view r = gsl_matrix_row(m,i);

		for(j=0;j<p;j++) arr[j] = gsl_vector_get(&r.vector, j);

		normal_transform(arr, arr, p);

		for(j=0;j<p;j++) gsl_vector_set(&r.vector, j, arr[j]);
	}

	delete[] arr;
}


size_t convert_RSS_to_metric(double metric[], const size_t size, const size_t n, const size_t p, double exit_point)
{
	size_t i, inx_min = string::npos;
	double sigma2, metric_min = DBL_MAX, t;

	if(p+size >= n)
	{
		cerr << "Error: Illegal number of features selected p + size >= n" << endl;
		exit(1);
	}

	sigma2 = metric[size-1]/(n-size-p);

	bool smallsigma2 = (fabs(sigma2) < EPS);

	if(verbose_flag && smallsigma2)
		cerr << "Warning: Sigma2 is too small. Use different form for Mallow's Cp calculation." << endl;

	for(i=0 ; i<size ; i++)
	{
		// Mallow's Cp
		if(smallsigma2){
			t = (metric[i] + 2*(p+i+1)*sigma2)/n;

			if(i==0)
				exit_point = exit_point*t;
			else
				if(t < exit_point) break;

		}else
			t = metric[i]/sigma2 + 2*(p+i+1) - n;

		metric[i] = t;

		if(t < metric_min)
		{
			metric_min = t;
			inx_min = i;
		}
	}

	return inx_min;
}


//------------------------------------------------------------------------------
// end Rabit.cpp auxillary functions
//------------------------------------------------------------------------------

//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

//------------------------------------------------------------------------------
// Main function call and documentation shown below now
//------------------------------------------------------------------------------


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

    int status, calculate;

    gsl_matrix *X, *Y, *B, *C;
    gsl_vector *beta, *sderr, *tt, *FDR, *RSS;

    size_t i, j, inx, count, parseCnt = 4, *index_array, n, p, k, m; //sudo random assignment here.


    // Background confounding factor matrices for Y
    set<string> Cfile_vec;
    set<string>::iterator Cfile_vec_iter;

    double FDR_thres = f, *arr, *metric_array, sigma2, exit_point =0.1;

    // transform Y vectors to normal distribution
    bool normal_transform = t, select_best = s, run_forward = r;

    vector<string>
        rownames,	// common row names across all matrices
        colnames		// common column names across Y and all its background matrices
        ;

    vector<rank_node<double> > FDR_sortvec;
    vector<rank_node<double> >::iterator FDR_sortvec_iter;
    set<string> included;

    // read in matrix maps from file
    NumericMatrix *X_map = new NumericMatrix, *Y_map = new NumericMatrix, *B_map = new NumericMatrix, *C_map;
    size_t append_size = 1 + Cfile_vec.size();

    // construct X background matrix B
    if(b == NumericMatrix(0))
    {
        // no X background, construct all 1 background as Intercept
        map<string, size_t>::iterator map_iter;

        // same row dimension with X
        B_map->Nrow = x->Nrow;
        B_map->row_map = x->row_map;
        B_map->rownames = x->rownames;

        // repeated column
        B_map->Ncol = append_size;
        for(i=0;i<append_size;i++) B_map->colnames.push_back("Intercept");

        for(i=0; i<B_map->Nrow; i++)
        {
            arr = new double[append_size];
            for(j=0;j<append_size;j++) arr[j] = 1;
            B_map->mat.push_back(arr);
        }

    }else{
        B_map->b
    }

}




