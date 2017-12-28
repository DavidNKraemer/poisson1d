#include "TriSolver.h"

// MAKE SURE YOU HAVE LAPACKE & CBLAS INCLUDED

namespace ams562 {

  // A_ = [2 2 2 ... 2 -1 -1 ... -1]
  // where N_ 2's and N_-1 -1's
  void TriSolver::assemble() {

    A_.clear();
    u_.clear();
    u_.resize(N_);
    A_.resize(N_ * (N_ - 1));

    for (int i = 0; i < (int) N_; i++)
      A_[i] = 2.;
    for (int i = 0; i < (int) N_ - 1; i++)
      A_[i + (int) N_] = -1.;
  }

  // call tridiagonal solver for SPD
  // keep in mind boundary conditions
  void TriSolver::solve(const std::vector<double> &b) { 
    configure_measurements(b);

    double *up = u_.data();
    double *Ap_diag = A_.data();
    double *Ap_subdiag = Ap_diag + N_;

    // Solve Au = b when A = (P^-1)LU
    LAPACKE_dptsv(LAPACK_ROW_MAJOR, N_, 1, Ap_diag, Ap_subdiag, up, 1);

    // insert the boundary values
    include_bcs();

  }

} // namespace ams562
