// ~\~ language=C++ filename=src/initial_conditions/normalize_power_spectrum.cc
// ~\~ begin <<adhesion_example.md|src/initial_conditions/normalize_power_spectrum.cc>>[0]
#include "initial_conditions.hh"
#include <iostream>

// ~\~ begin <<appendix.md|gsl-integrate-qagiu>>[0]
#include <gsl/gsl_integration.h>

using Function = std::function<double (double)>;

double integration_helper(double x, void *params)
{
  auto f = reinterpret_cast<Function *>(params);
  return (*f)(x);
}

template <typename F>
double integrate_qagiu(
    F const &func, double lower,
    double epsabs, double epsrel, size_t limit=1024)
{
  double x, abserr;
  Function integrant = func;
  gsl_integration_workspace *workspace =
    gsl_integration_workspace_alloc(limit);

  gsl_function f;
  f.function = &integration_helper;
  f.params = reinterpret_cast<void *>(&integrant);

  gsl_integration_qagiu(
    &f, lower, epsabs, epsrel, limit,
    workspace, &x, &abserr);

  gsl_integration_workspace_free(workspace);
  return x;
}
// ~\~ end
// ~\~ begin <<adhesion_example.md|top-hat-function>>[0]
double W_th(double y) {
  return 3.0 / pow(y, 3) * (sin(y) - y * cos(y));
}
// ~\~ end

PowerSpectrum normalize_power_spectrum(
    BoxParam const &box,
    PowerSpectrum const &P,
    Config const &cosmology)
{
  double k_lower = 2 * M_PI / box.L;
  double epsabs = 1e-6, epsrel = 1e-6;
  double sigma8 = cosmology["sigma8"].as<double>();

  // ~\~ begin <<adhesion_example.md|define-integrand>>[0]
  auto integrand = [&] (double k) {
    return P(k) / (2 * M_PI*M_PI) * pow(W_th(8.0 * k) * k, 2);
  };
  // ~\~ end

  double x = integrate_qagiu(
    integrand, k_lower, epsabs, epsrel);

  double A = sigma8 * sigma8 / x;
  std::clog << "Normalised power spectrum, A = " << A << ".\n";

  return [=] (double k) {
    return A * P(k);
  };
}
// ~\~ end
