#ifndef __LUSOLVER_H
#define __LUSOLVER_H

#include "BaseSolver.h"
#include <lapacke.h>
#include <cblas.h>


namespace ams562 {

class LuSolver : public BaseSolver {
public:
  LuSolver(double xl, double xr, unsigned int nx) : BaseSolver(xl, xr, nx) {}
  ~LuSolver() {}

  virtual void assemble() override;
  virtual void solve(const std::vector<double> &b) override;
};

} // namespace ams562

#endif
