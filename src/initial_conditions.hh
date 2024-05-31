// ~/~ begin <<adhesion_example.md#src/initial_conditions.hh>>[init]
#pragma once
#include "boxparam.hh"
#include <vector>
#include <yaml-cpp/yaml.h>
#include <functional>

using PowerSpectrum = std::function<double (double)>;
using Config = YAML::Node;

extern std::vector<double>
generate_white_noise(
    BoxParam const &box,
    unsigned long seed);

extern PowerSpectrum EisensteinHu(
    Config const &cosmology);

extern PowerSpectrum normalize_power_spectrum(
    BoxParam const &box,
    PowerSpectrum const &P,
    Config const &cosmology);

extern void compute_potential(
    BoxParam const &box,
    std::vector<double> &white_noise,
    PowerSpectrum const &P);

extern std::vector<double>
generate_initial_potential(
    BoxParam const &box,
    Config const &config);
// ~/~ end