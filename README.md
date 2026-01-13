# Discretised Hydrogen Atom Hamiltonian Solver
This project constructs and solves a discretised 1D radial Schrodinger equation for the hydrogen atom to obtain the lowest eigenvalues, using Eigen and Spectra.

The radial position, $r$, is in units of $\frac{4\pi\epsilon_0\hbar^2}{m_ee^2}$ or simply $a_0$. 

The energy eigenvalues, $E$, are in units of $\frac{\hbar^2}{2ma_0^2}$.

## Features
* Builds a 1D linear radial grid.
* Constructs a discretised sparse hamiltonian.
* Uses Spectra's symmetric eigenvalue shift solver.
* Computes and outputs lowest eigenvalues

## The Hamiltonian
The *Hamiltonian* is simply the sum of the kinetic and potential components $H = K +  V$. This is written in a tridiagonal reprasentation so the operators were built as sparse matrices using `Eigen::SparseMatrix<double>` to keep memory usage low.


The *kinetic energy* term of the Hamiltonian is a tri-diagonal matrix, with main and off diagonals given by\
$K_{i,i} = \frac{2}{dr^2}$ and 
$K_{i,i+1} = K_{i,i-1} = \frac{-1}{dr^2}$, where $dr$ is the spacing of the linear radial grid.


The *potential* term is a diagonal matrix, when elements are equal to the potential at the $ith$ position. This potential is
$V(r_i) = \frac{-2}{r_i} + \frac{l(l+1)}{r_i^2}$.

## Spectra Usage
The wrapper class `SparseSymShiftSolve` was used to apply a shift-invert operation. This focused the eigensolver on values near the shift, $\sigma$. The corresponding eigensolver used was `SymEigsShiftSolver`. These deal with sparse representations.

## Building and Running
Requires the header-only C++ libraries **Eigen** and **Spectra**.

To build using cmake, with the [CMakeLists.txt file](CMakeLists.txt):
```
mkdir build && cd build
cmake ..
make
```
This creates an executable called hydrogen.exe. 
```
./hydrogen
```


