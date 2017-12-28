#include "PpSolver.h"

#include <iostream>

using namespace std;

namespace ams562 {

  // assemble the packed UPPER part of A
  // A = [2 -1 0 0]
  //     [x 2 -1 0]
  //     [x x 2 -1]
  //     [x x x  2]
  // where x means empty.
  void PpSolver::assemble() {
    int entries = (N_ * (N_+1)) / 2;

    A_.clear();
    A_.resize(entries);
    u_.clear();
    u_.resize(N_);

    unsigned int idx = 0;
    for (unsigned int j = 0; j < N_; j++) {
      for (unsigned int i = 0; i < j+1; i++) {
        if (i == j) {
          A_[idx] = 2.;
        }
        else if (j - i == 1) {
          A_[idx] = -1.;
        }
        else {
          A_[idx] = 0.;
        }
        idx++;
      }
    }

  }

  // call packed SPD solver, keep in mind
  // boundary conditions
  void PpSolver::solve(const std::vector<double> &b) {

    configure_measurements(b);
    double *Ap = A_.data();
    double *up = u_.data();

    // unsigned int idx = 0;
    // for (idx = 0; idx < N_ * (N_+1)/2; idx++) 
    //   cout << Ap[idx] << endl;

    // Solve Au = b when A = (P^-1)LU
    lapack_int result = LAPACKE_dppsv(LAPACK_COL_MAJOR, 'U', N_, 1, Ap, up, N_);
    cout << "LAPACKE Result: " << (int) result << endl;

    include_bcs();

  }

} // namespace ams562
