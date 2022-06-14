// ~\~ language=C++ filename=src/main.cc
// ~\~ begin <<adhesion_example.md|src/main.cc>>[init]
#include <iostream>
#include <argagg/argagg.hpp>
#include <yaml-cpp/yaml.h>

#include "run.hh"

// ~\~ begin <<adhesion_example.md|version>>[init]
#define VERSION "2.0"
// ~\~ end
// ~\~ begin <<adhesion_example.md|main-arguments>>[init]
argagg::parser argparser {{
  { "help",    {"-h", "--help"},
    "Show this help message.", 0 },
  { "version", {"--version"},
    "Show the software version.", 0 },
  { "config",  {"-c", "--config"},
    "Supply configuration file.", 1 },
  { "defaults", {"-d", "--defaults"},
    "Show default configuration.", 0 }
}};
// ~\~ end

const char *default_config = R"YAML(
# ~\~ begin <<adhesion_example.md|default-config>>[init]
# Default configuration

cosmology:
  h:        0.674   # Hubble parameter / 100
  ns:       0.965   # primordial power spectrum index
  Omega0:   1.0     # density in units of critical density
  sigma8:   0.811   # amplitude over 8 Mpc/h

initial-conditions:
  random:
    seed:              8
    resolution:      128       # logical box size
    box-size:         50.0     # physical box size in Mpc/h
    power-spectrum:  Eisenstein & Hu (no baryons)

run:
  time:     [0.2, 0.5, 1.0]

output:
  hdf5:            output/lcdm.h5
  walls:           output/lcdm-{time:02.1f}-walls.obj
  filaments:       output/lcdm-{time:02.1f}-filaments.obj
  threshold:       1.0
# ~\~ end
)YAML";

int main(int argc, char **argv)
{
  // ~\~ begin <<adhesion_example.md|parse-arguments>>[init]
  argagg::parser_results args;
  try {
    args = argparser.parse(argc, argv);
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  // ~\~ end
  // ~\~ begin <<adhesion_example.md|maybe-print-help>>[init]
  if (args["help"]) {
    std::cout << "Adhesion model example code -- (C) 2018 Johan Hidding\n";
    std::cout << argparser;
    return EXIT_SUCCESS;
  }
  // ~\~ end
  // ~\~ begin <<adhesion_example.md|maybe-print-version>>[init]
  if (args["version"]) {
    std::cout << "amec v" << VERSION << "\n";
    return EXIT_SUCCESS;
  }
  // ~\~ end
  // ~\~ begin <<adhesion_example.md|maybe-print-defaults>>[init]
  if (args["defaults"]) {
    std::cout << default_config;
    return EXIT_SUCCESS;
  }
  // ~\~ end
  // ~\~ begin <<adhesion_example.md|load-config>>[init]
  YAML::Node config;
  if (args["config"]) {
    auto config_file = args["config"].as<std::string>();
    std::clog << "Reading `" << config_file << "` for input.\n";
    config = YAML::LoadFile(config_file);
  } else {
    std::clog << "No configuration given, proceeding with defaults.\n";
    config = YAML::Load(default_config);
  }
  // ~\~ end
  // ~\~ begin <<adhesion_example.md|run>>[init]
  run(config);
  // ~\~ end
  return EXIT_SUCCESS;
}
// ~\~ end
