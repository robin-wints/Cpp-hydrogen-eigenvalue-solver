#include <Eigen/Sparse>

std::pair<Eigen::VectorXd, double> grid(double start, double end, int N);

Eigen::SparseMatrix<double> build_kinetic(int N, double step);

Eigen::SparseMatrix<double> build_potential(int N, const Eigen::VectorXd& x, int l = 0);