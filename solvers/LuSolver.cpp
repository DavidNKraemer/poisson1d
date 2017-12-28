#include "LuSolver.h"

// MAKE SURE YOU HAVE LAPACKE & CBLAS INCLUDED
#include <lapacke.h>
#include <cblas.h>
#include <iostream>

namespace ams562 {

  // assemble a whole dense matrix of N_ by N_
  void LuSolver::assemble() {

    A_.clear();
    A_.resize(N_ * N_);
    u_.clear();
    u_.resize(N_);

    for (int i = 0; i < (int) N_; i++) {
      for (int j = 0; j < (int) N_; j++) {
        if (i == j) {
          A_[IND(i,j,N_)] = 2.;
        }
        else if (abs(i - j) == 1) {
          A_[IND(i,j,N_)] = -1.;
        }
        else {
          A_[IND(i,j,N_)] = 0.;
        }
      }
    }
  }

  // given a right-hand side vector b, solve the problem
  // remember to apply the boundary conditions, i.e. add
  // lb_ (or lb_/h^2 depending on implementation) to b[0]
  // rb_ (or rb_/h^2) to b[N_-1]
  // then call lapack solver to solve the system
  // you can create a vector<int> as buffer for pivotings
  void LuSolver::solve(const std::vector<double> &b) {

    configure_measurements(b);
    double *Ap = A_.data();
    double *up = u_.data();
    int pivots[N_];

    // Solve Au = b when A = (P^-1)LU
    LAPACKE_dgesv(LAPACK_ROW_MAJOR,  N_, 1, Ap, N_, pivots, up, 1);

    // insert the boundary values
    include_bcs();
  }

} // namespace ams562
