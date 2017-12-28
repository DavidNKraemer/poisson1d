#ifndef __TRISOLVER_H
#define __TRISOLVER_H

#include "BaseSolver.h"

namespace ams562 {

class TriSolver : public BaseSolver {
public:
  TriSolver(double xl, double xr, unsigned int nx) : BaseSolver(xl, xr, nx) {}

  virtual void assemble() override;
  virtual void solve(const std::vector<double> &b) override;
};

} // namespace ams562

#endif
