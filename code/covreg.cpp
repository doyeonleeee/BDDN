#include <omp.h>
#include <math.h>
#include <iostream>
#include <RcppArmadillo.h>
#include <random>

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
void printDimensions_cub(cube a) {
  Rcout << "cubDimensions: " << a.n_rows << " " << a.n_cols << " " << a.n_slices << std::endl;
}
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
void printDimensions_mat(mat a) {
  Rcout << "matDimensions: " << a.n_rows << " " << a.n_cols << std::endl;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// Function to generate random matrix from a Wishart distribution
arma::mat rwish(int nu, const arma::mat& S0) {
  // printDimensions_mat(S0);
  arma::mat sS0 = arma::chol(S0);
  arma::mat Z = arma::randn(nu, S0.n_rows) * sS0;
  return Z.t() * Z;
}

// Function to generate random matrix from an inverse Wishart distribution
arma::mat riwish(int v, const arma::mat& S) {
  return arma::inv(rwish(v, arma::inv(S)));
}

// Function to calculate the product of three matrices
arma::mat tbf(const arma::mat& Y, arma::mat xg, const arma::mat& A){ 

  arma::mat result = trans(Y) * xg * inv(trans(xg) * xg + A);

  return result;
}

// Function to perform matrix multiplication
arma::mat xb(const arma::cube& xx, const arma::cube& b) {
  arma::mat XB = zeros(xx.n_rows,b.n_rows); 
  int d = b.n_slices;

  if (is_finite(d) && d > 0) {
    for (int i = 0; i < d; ++i) {

      XB += xx.slice(i) * trans(b.slice(i));
    }
  }
  return XB;
}


// Function to calculate the sum of squared differences
arma::mat ts(const arma::mat& Y, const arma::cube& xg, const arma::cube& B) {
  arma::mat d = Y - xb(xg, B);

  return trans(d) * d;
}

// Function to generate random matrix from a multivariate normal distribution
arma::mat rmn(const arma::mat& M, const arma::mat& Srow, const arma::mat& Scol) {
  arma::vec eigval_row;
  arma::mat eigvec_row;
  arma::eig_sym(eigval_row, eigvec_row, Srow);
  
  arma::vec eigval_col;
  arma::mat eigvec_col;
  arma::eig_sym(eigval_col, eigvec_col, Scol);
  
  arma::mat tmp_row = eigvec_row * diagmat(sqrt(eigval_row)) * trans(eigvec_row);
  arma::mat tmp_col = eigvec_col * diagmat(sqrt(eigval_col)) * trans(eigvec_col);
  
  arma::mat Z = arma::randn(Srow.n_rows, Scol.n_rows); 
  
  return tmp_row * Z * tmp_col + M;
}

// Function to update g
arma::vec rg_fc(const arma::mat& Y, const arma::mat& Bi, const arma::cube& Bmi,
                const arma::mat& A, const arma::mat& x, const arma::cube& xgmi, int n) {
  arma::mat bx = x * trans(Bi);
  arma::mat Ai = arma::inv(A);
  vec s2 = 1.0 / ((diagvec(bx * Ai * trans(bx))) + 1.0);
  arma::vec m = diagvec(bx * Ai * trans(Y - xb(xgmi, Bmi))) % s2;
  
  return m + sqrt(s2) % randn(n);
}


// Function to update B
arma::mat rB_fc(const arma::mat& Y, const arma::mat& xgi, const arma::cube& xgmi,
                const arma::cube& Bmi, const arma::mat& A, int q, const arma::mat& V0) {
  
  arma::mat tB = tbf(Y - xb(xgmi, Bmi), xgi, V0);

  arma::mat s = arma::inv(trans(xgi) * xgi + V0);
  return rmn(tB, A, s);
} 

// Function to update A
arma::mat rA_fc(const arma::mat& Y, const arma::cube& xg, const arma::cube& B,
                int n, int p, int r, double nu0, const arma::mat& A0, const arma::mat& V0) {
  
  arma::mat s = ts(Y, xg, B);
  for (int i = 0; i < r; ++i) {
    s += B.slice(i) * V0 * trans(B.slice(i));
  }
  s += A0;
  return riwish(n + r * p + nu0, s);
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
// Main MCMC function
List MCMC(const arma::mat& Y, const arma::mat& x, int r = 1,
                 int niter = 1000, int nthin = 1) {
  
  int n = Y.n_rows;
  int p = Y.n_cols;
  int q = x.n_cols;
  
  // Priors
  double gg = 1.0 / n;
  double nu0 = p + 2;
  arma::mat A0 = cov(Y);
  arma::mat V0 = gg * trans(x) * x;
  
  // Starting values
  arma::cube B = randn(p,q,r);
  arma::mat g = randn(n,r); 
  arma::cube xg(n, q, r);
  arma::cube haha = randn(p,q,r);
  
  for (int i = 0; i < r; ++i) {
    xg.slice(i) = x.each_col() % g.col(i);
  }
  arma::mat A = cov(Y);
  
   
  // Output
  arma::cube b_save(p, q, niter / nthin); 
  arma::cube a_save(p, p, niter / nthin);
   
  for (int ns = 1; ns <= niter; ++ns) {
    if(ns % 100 == 0){
      Rcout << "ntier: " << ns << std::endl;
    }
    for (int i = 0; i < r; ++i) {
      arma::cube B1 = B;
      arma::cube xg1 = xg;
      B1.shed_slice(i);  
      xg1.shed_slice(i);
      
      g.col(i) = rg_fc(Y, B.slice(i), B1, A, x, xg1, n);
      
      for (int j = 0; j < x.n_cols; ++j){
       xg.slice(i).col(j) = x.col(j) % g.col(i);
      }
    }
    
    A = rA_fc(Y, xg, B, n, p, r, nu0, A0, V0); //Psi
    
    for (int i = 0; i < r; ++i) {
      arma::cube B2 = B;
      arma::cube xg2 = xg;
      B2.shed_slice(i);  
      xg2.shed_slice(i);
      
      B.slice(i) = rB_fc(Y, xg.slice(i), xg2, B2, A, q, V0);
    }
    
    
    if (ns % nthin == 0) {
      b_save.slice(ns/nthin-1) = B.slice(0);
      a_save.slice(ns/nthin-1) = A;
    }
     
  }
  
  return List::create(Named("B.psamp") = b_save, 
                      Named("A.psamp") = a_save);
}
