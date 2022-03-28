// ~\~ language=C++ filename=src/initial_conditions/eisenstein-hu.cc
// ~\~ begin <<adhesion_example.md|src/initial_conditions/eisenstein-hu.cc>>[0]
#include "initial_conditions.hh"

PowerSpectrum EisensteinHu(Config const &cosmology)
{
  double const
    e         = exp(1),
    Theta_CMB = 2.7255/2.7,
    Omega0    = cosmology["Omega0"].as<double>(),
    h         = cosmology["h"].as<double>(),
    ns        = cosmology["ns"].as<double>();

  return [=] (double k)
  {
    double q  = k * pow(Theta_CMB, 2)/(Omega0 * h),
           L0 = log(2*e + 1.8*q),
           C0 = 14.2 + 731.0/(1 + 62.5*q),
           T0 = L0 / (L0 + C0 * pow(q, 2));

    return pow(k, ns) * pow(T0, 2);
  };
}
// ~\~ end
