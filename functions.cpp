#include <Eigen/Sparse>
#include "functions.hpp"

std::pair<Eigen::VectorXd, double> grid(double start, double end, int N){
  
  Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(N, start, end);

  double step = (end - start) / (N - 1);

  return {x, step};

}

Eigen::SparseMatrix<double> build_kinetic(int N, double step){
  
  Eigen::SparseMatrix<double> kinetic(N,N); // N by N matrix

  const double main_diag = 2 / (step * step);
  const double off_diag = -1 / (step * step);

  // Insert main and off diagonal elements
  for (int i = 0; i < N; ++i){

    kinetic.insert(i,i) = main_diag;
    
    // Avoid inserting off diags out of bounds
    if (i < N - 1){

    kinetic.insert(i, i+1) = off_diag;
    kinetic.insert(i+1, i) = off_diag;

    }

  }

  return kinetic;

}

Eigen::SparseMatrix<double> build_potential(int N, const Eigen::VectorXd& x, int l){
  
  Eigen::SparseMatrix<double> potential(N,N); // N by N matrix

  for (int i = 0; i < N; ++i){
    double V = -2/x[i] + (l*(l+1)) / (x[i] * x[i]);
    potential.insert(i,i) = V;
  }

  return potential;

}