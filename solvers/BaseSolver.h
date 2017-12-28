#ifndef __BASESOLVER_H
#define __BASESOLVER_H

/* Numerical linear algebra libraries BLAS and LAPACK */
#include <cblas.h>
#include <lapacke.h>

#include <vector>

#define IND(i,j,cols) ((i) * (cols) + (j))

namespace ams562 {

  class BaseSolver {
    public:
      BaseSolver(double xl, double xr, unsigned int nx)
        : xl_(xl), xr_(xr), N_(nx - 2) {
          h_ = (xr_ - xl_) / (N_ + 1);
          u_.resize(N_); // you should do the same thing in assemble
        }

      virtual ~BaseSolver() {}

      // pure abstract methods
      virtual void assemble() = 0;
      virtual void solve(const std::vector<double> &b) = 0;

      // implementation needed, assign lhs to lb_, rhs to rb_
      void assign_bcs(double lhs, double rhs) {
        lb_ = lhs;
        rb_ = rhs;
      }

      // implemetation needed, return u_
      std::vector<double> &get_solution() {
        return u_;
      }

    protected:
      double xl_, xr_; // left- and rightmost grid points xl_ = a, xr_ = b
      double h_; // grid step size
      unsigned int N_; // linear system size (n - 1)
      double lb_, rb_; // boundary values lb_ = u(a), rb_ = u(b)
      std::vector<double> u_; // solution vector
      std::vector<double> A_; // system matrix

      void configure_measurements(const std::vector<double> &b) {

        double sig = h_ * h_;

        /* First, copy the data of b into u_ */
        const double *bp = b.data();
        double *up = u_.data();

        cblas_daxpy(N_, sig, bp, 1, up, 1); // u_[i] = sig * b[i], i = 0,...,N_-1

        // boundary conditions
        up[0] += lb_;
        up[N_-1] += rb_;
      }

      void include_bcs() {
        // insert the boundary values
        u_.insert(u_.begin(), lb_);
        u_.push_back(rb_);
      }
  };

} // namespace ams562

#endif
