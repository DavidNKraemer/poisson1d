#ifndef __PPSOLVER_H
#define __PPSOLVER_H

#include "BaseSolver.h"

namespace ams562 {

class PpSolver : public BaseSolver {
public:
  PpSolver(double xl, double xr, unsigned int nx) : BaseSolver(xl, xr, nx) {}

  virtual void assemble() override;
  virtual void solve(const std::vector<double> &b) override;
};

} // namespace ams562

#endif
