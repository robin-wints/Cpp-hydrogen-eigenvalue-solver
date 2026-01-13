/* Constructs and solves a sparse discretised Hamiltonian for the hydrogen atom
 to obtain the necessary number of eigenvalues, using Eigen and Spectra. */

#include <iostream>
#include <Eigen/Sparse>
#include <Spectra/SymEigsShiftSolver.h>
#include <Spectra/MatOp/SparseSymShiftSolve.h>
#include "functions.hpp"

int main() {

  // Define variables
  const int n_max = 5; // Largest quantum num. of interest 
  const double r_min = 0.005; // r=0 is a singularity.
  const double r_max = 6 * n_max;
  const int N = 1500; // No. of points on 1D grid
  const int ncv = 6 * n_max; // Convergence speed parameter
  const double sigma = -1.0; // Shift parameter

  // Get 1d grid
  const auto [r, dr] = grid(r_min, r_max, N); 

  // Construct Hamiltonian
  Eigen::SparseMatrix<double> K = build_kinetic(N, dr);
  Eigen::SparseMatrix<double> V = build_potential(N, r, 0);
  Eigen::SparseMatrix<double> H = K + V;  

  // Construct matrix operation object
  Spectra::SparseSymShiftSolve<double> op(H); 

  // Construct eigen solver object
  Spectra::SymEigsShiftSolver<Spectra::SparseSymShiftSolve<double>> 
    eigs(op, n_max, ncv, sigma);

  // Initialize and compute
  eigs.init();
  int nconv = eigs.compute(Spectra::SortRule::LargestMagn);
  std::cout << "Number of eigenvalues wanted: " << n_max << std::endl;
  std::cout << "Number of eigenvalues converged: " << nconv << std::endl;

  // Retrieve and print results
  Eigen::VectorXd evalues = eigs.eigenvalues();
  std::cout << evalues.reverse() << std::endl;

  return 0;
    
}