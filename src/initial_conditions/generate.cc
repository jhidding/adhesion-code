// ~/~ begin <<adhesion_example.md#src/initial_conditions/generate.cc>>[init]
#include "initial_conditions.hh"
#include <iostream>

std::vector<double> generate_initial_potential(
    BoxParam const &box,
    Config const &config)
{
  std::clog << "# Generating white noise with seed:\n"
            << config["initial-conditions"]["random"]["seed"] << "\n";
  auto seed = config["initial-conditions"]["random"]["seed"].as<unsigned long>();
  auto field = generate_white_noise(box, seed);

  std::clog << "Applying power spectrum with cosmology:\n"
            << config["cosmology"] << "\n";
  auto cosmology = config["cosmology"];
  auto power_spectrum = normalize_power_spectrum(
    box, EisensteinHu(cosmology), cosmology);
  compute_potential(box, field, power_spectrum);

  return field;
}
// ~/~ end