#include <Rcpp.h>
#include <omp.h>
using namespace Rcpp;
extern int openmp_threads;
using namespace std;
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

// [[Rcpp::export]]
NumericVector colMin(NumericMatrix X) {
  int ncol = X.ncol();
  Rcpp::NumericVector out(ncol);
  for (int col = 0; col < ncol; col++){
    out[col]=Rcpp::min(X(_, col));
  }
  return wrap(out);
}

// [[Rcpp::export]]
double measr_cpp(NumericMatrix x ) {
  x=transpose(x);
  int nr = x.nrow(), nc1 = x.ncol();


  NumericMatrix dist(nc1, nc1);
  bool symmetric = 1;
  NumericVector mindist(nc1);
#pragma omp parallel for                                                     \
  if (openmp_threads > 1 && (nc1 + 0.0) * (nc1 + 0.0) * (nr + 0.0) > 1000)   \
    num_threads(openmp_threads)                                              \
    shared(dist)
    for (int col2 = 0; col2 < nc1; col2++) {
      NumericVector tmp(nr);
      double accum;
      int col1_max = (symmetric) ? col2 + 1 : nc1;
      for (int col1 = 0; col1 < col1_max; col1++) {
        NumericMatrix::Column vx = x(_, col1);  // column <col1> of first matrix
        NumericMatrix::Column vy = x(_, col2);  // column <col2> of second matrix

        accum = sum((vx - vy) * (vx - vy));
        dist(col1, col2) = sqrt(accum);
        dist(col2, col1) = dist(col1, col2);
        dist(col2, col2) = 1000;
      }
    }

    mindist= colMin(dist);
  double gammabar= mean(mindist);
  double s= sum((mindist-gammabar)*(mindist-gammabar));
  double cvd = (sqrt(s/nc1))/gammabar;
  return cvd;
}

// [[Rcpp::export]]
NumericMatrix measr_ff(NumericMatrix x ) {
  x=transpose(x);
  int nr = x.nrow(), nc1 = x.ncol();


  NumericMatrix dist(nc1, nc1);
  bool symmetric = 1;
  NumericVector mindist(nc1);
#pragma omp parallel for                                                     \
  if (openmp_threads > 1 && (nc1 + 0.0) * (nc1 + 0.0) * (nr + 0.0) > 1000)   \
    num_threads(openmp_threads)                                              \
    shared(dist)
    for (int col2 = 0; col2 < nc1; col2++) {
      NumericVector tmp(nr);
      double accum;
      int col1_max = (symmetric) ? col2 + 1 : nc1;
      for (int col1 = 0; col1 < col1_max; col1++) {
        NumericMatrix::Column vx = x(_, col1);  // column <col1> of first matrix
        NumericMatrix::Column vy = x(_, col2);  // column <col2> of second matrix

        accum = sum((vx - vy) * (vx - vy));
        dist(col1, col2) = sqrt(accum);
        dist(col2, col1) = dist(col1, col2);
        dist(col2, col2) = 1000;
      }
    }


    return dist;
}


#ifdef _OPENMP
#include <omp.h>
#endif

int openmp_threads = 1;

// [[Rcpp::export]]
DataFrame CPP_get_openmp_threads() {
  int num_threads = openmp_threads;
#ifdef _OPENMP
  int max_threads = omp_get_max_threads();
#else
  int max_threads = 0;
#endif
  DataFrame res =
    DataFrame::create(_["available"] = max_threads > 0,
                      _["max"] = max_threads,
                      _["threads"] = num_threads);
  res.attr("row.names") = "OpenMP";
  return res;
}

// [[Rcpp::export]]
void CPP_set_openmp_threads(int n) {
  if (n < 1) stop("internal error -- number of threads must be >= 1");
#ifdef _OPENMP
  int max_threads = omp_get_max_threads();
  if (n > max_threads) n = max_threads;
  openmp_threads = n;
#else
  if (n > 1) Rf_warning("OpenMP support not available");
  openmp_threads = 1;
#endif
}
