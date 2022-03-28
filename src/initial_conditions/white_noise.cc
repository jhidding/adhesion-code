// ~\~ language=C++ filename=src/initial_conditions/white_noise.cc
// ~\~ begin <<adhesion_example.md|src/initial_conditions/white_noise.cc>>[0]
#include "initial_conditions.hh"
#include <random>

std::vector<double>
generate_white_noise(
    BoxParam const &box, unsigned long seed)
{
  auto result = std::vector<double>(box.size());

  std::mt19937 random(seed);
  std::normal_distribution<double> normal;

  for (double &value : result) {
    value = normal(random);
  }

  return result;
}
// ~\~ end
