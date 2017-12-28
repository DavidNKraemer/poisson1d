#include <algorithm>
#include <assert.h>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <string>
#include <tuple>
#include "poisson1d.h"

#define PI 3.14159265358979323846264338 // excessive

using namespace ams562;
using namespace std; // laziness

using UnaryDoubleFn = function<double(double)>; // useful shortcut


/* Function headers */
vector<double> linspace(double lower, double upper, int n);
vector<double> measurements(double lower, double upper, int n, UnaryDoubleFn ufun);
vector<double> map_linspace(double lower, double upper, int n, UnaryDoubleFn ufun);
double f0(double x);
double f1(double x);
double f2(double x);
double f3(double x);
double nf0dd(double x);
double nf1dd(double x);
double nf2dd(double x);
double nf3dd(double x);
void run_solvers(int n, double lower, double upper);


/* Main program */
int main() {

  // number of points, lower bound, upper bound
  run_solvers(202, 0., 1.);
  return 0;
}


/* Function implementations */
vector<double> linspace(double lower, double upper, int n) {
  // standard linspace a la NumPy or Julia's implementation
  assert(n > 1 && "linspace requires at least two points");
  vector<double> x(n);

  auto width = (upper - lower) / (n-1);
  x[0] = lower;

  for (int i = 1; i < n; i++)
    x[i] = x[0] + i * width;

  return x;
}

vector<double> measurements(double lower, double upper, int n, UnaryDoubleFn ufun) {
  // compute the measurements for the problem
  assert(n > 2 && "Insufficiently long measurement vector");
  auto big_b = linspace(lower, upper, n);
  transform(big_b.begin(), big_b.end(), big_b.begin(), ufun);
  vector<double> b(big_b.begin()+1, big_b.end() - 1);

  return b;
}

double f0(double x) {
  // f0(x) = x^2
  return x * x;
}

double nf0dd(double x) {
  // -f0''(x) = -2
  x += 0.; // hack, suppressing an unused parameter warning
  return -2.;
}

double f1(double x) {
  // f1(x) = 2 x^2
  return 2 * x * x;
}

double nf1dd(double x) {
  // -f1''(x) = -4
  x += 0.; // hack, suppressing an unused parameter warning
  return -4.;
}

double f2(double x) {
  // f2(x) = - sin(Pi / 2. * x), 
  return sin(PI / 2. * x);
}

double nf2dd(double x) {
  // -f2''(x) = Pi^2 / 4 * sin(PI / 2. * x)
  return PI * PI / 4. * sin(PI / 2. * x);
}

double f3(double x) {
  return exp(-x) * sin(10 * x);
}

double nf3dd(double x) {
  return exp(-x) * (99 * sin(10*x) + 20 * cos(10 * x));
}

vector<double> map_linspace(double lower, double upper, int n, UnaryDoubleFn f) {
  // this maps f to the uniform grid of points {lower, ..., upper}
  vector<double> x = linspace(lower, upper, n);
  vector<double> y(x.size());
  transform(x.begin(), x.end(), y.begin(), f);
  return y;
}

void run_solvers(int n, double lower, double upper) {
  cout << endl << "Poisson 1D Solver: u'' = -f" << endl;
  cout << n << "-point uniform grid" << endl;
  cout << "Interval: (" << lower << "," << upper << ")" << endl << endl;

  /* 
   * This setup lets me just loop through the different solvers and the
   * different functions, performing the same task repeatedly. It's so that
   * this driver program isn't 6x longer than it currently is.
   */
  vector<tuple<BaseSolver*, string*, string*, double, double>> solvers;
  vector<tuple<UnaryDoubleFn, UnaryDoubleFn, string*, string*>> functions;

  // solvers stores information relevant to the linear system solver method
  // each element of solvers is a tuple containing
  //  - a BaseSolver pointer
  //  - a descriptive string for the particular solver
  //  - an abbreviated string for naming a file
  //  - a double for the 2-norm error
  //  - a double for the infinity-norm error
  solvers.push_back(
      make_tuple(new LuSolver(lower, upper, n), 
        new string("uncompressed matrix"), new string("lu"), 0., 0.));
  solvers.push_back(
      make_tuple(new PpSolver(lower, upper, n), 
        new string("packed dense matrix"), new string("pp"), 0., 0.));
  solvers.push_back(
      make_tuple(new TriSolver(lower, upper, n), 
        new string("tridiagonal matrix"), new string("tri"), 0., 0.));

  // functions stores information relevant to each of the functions being
  // integrated. each element of functions is a tuple containing
  //  - a function representing the exact solution to the ODE
  //  - a function representing -f'' in the equation u'' = -f''
  //  - a descriptive string for the function
  //  - an abbreviated string for naming a file
  functions.push_back(
      make_tuple(f0, nf0dd, new string("x^2"), new string("f0")));
  functions.push_back(
      make_tuple(f1, nf1dd, new string("2x^2"), new string("f1")));
  functions.push_back(
      make_tuple(f2, nf2dd, new string("sin(pi/2 x)"), new string("f2")));
  functions.push_back(
      make_tuple(f3, nf3dd, new string("exp(-x) sin(10x)"), new string("f3")));


  // uniform grid points for integration
  auto x = linspace(lower, upper, n);

  for (auto func : functions) { // for each function
    UnaryDoubleFn f = get<0>(func);
    UnaryDoubleFn nfdd = get<1>(func);

    // measurements points determined by the 2nd derivative
    auto b = measurements(lower, upper, n, nfdd);

    // exact solution using the given function
    vector<double> y_exact = map_linspace(lower, upper, n, f);
    vector<double> err2(n);
    vector<double> errinf(n);


    for (auto solver : solvers) { // for each solver
      cout << "Solving with Dirichlet condition: f(x) = " << *get<2>(func);
      cout << " by means of an " << *get<1>(solver) << "." << endl;

      // this part actually computes the solution to the linear system using
      // the same interface for each solver (hooray for inheritence!)
      auto method = get<0>(solver);
      method->assign_bcs(f(lower), f(upper));
      method->assemble();
      method->solve(b);
      auto y = method->get_solution();

      // not required by the project description, however I wanted to be able
      // to plot these nicely and format them in the README.md
      string filename= "data/" + *get<2>(solver) + "-" + *get<3>(func) + ".txt";
      cout << "Writing to '" << filename << "'." << endl;
      write3v2txt(filename, x, y, y_exact);

      // Compute the difference vector, store it in y_exact
      for (unsigned int k = 0; k < err2.size(); k++) {
        err2[k] = y[k] - y_exact[k];
        errinf[k] = y[k] - y_exact[k];
      }

      // compute the 2norm error
      auto n2_err = get<3>(solver);
      n2_err = cblas_dnrm2(err2.size(), err2.data(), 1);

      // compute the inf norm error
      auto ninf_err = get<4>(solver);
      ninf_err = errinf[cblas_idamax(errinf.size(), errinf.data(), 1)];

      cout << "\t2-norm error: " << n2_err << endl;
      cout << "\tInfinity-norm error: " << ninf_err << endl << endl;
    }
  }


}
