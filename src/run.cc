// ~\~ language=C++ filename=src/run.cc
// ~\~ begin <<adhesion_example.md|src/run.cc>>[init]
#include <iostream>
#include <fstream>
#include <exception>
#include <filesystem>
#include <tuple>
#include <H5Cpp.h>
#include <fmt/format.h>

#include "run.hh"
#include "initial_conditions.hh"
#include "adhesion.hh"
#include "mesh_manipulation.hh"
#include "plane.hh"
#include "sphere.hh"
#include "writers.hh"
#include "write_obj.hh"

namespace fs = std::filesystem;

std::tuple<BoxParam,std::vector<double>> random_ic(YAML::Node const &config)
{
  std::clog << "# Using box with parameters:\n"
            << config["initial-conditions"]["random"] << "\n";
  BoxParam box(
    config["initial-conditions"]["random"]["resolution"].as<int>(),
    config["initial-conditions"]["random"]["box-size"].as<double>());
  auto potential = generate_initial_potential(
    box, config);
  return std::make_tuple(box, potential);
}

std::tuple<BoxParam,std::vector<double>> load_ic(YAML::Node const &config)
{
  auto input_filename = config["initial-conditions"]["file"].as<std::string>();
  H5::H5File input_file(input_filename);

  unsigned N = read_attribute<unsigned>(input_file, "N");
  double L = read_attribute<double>(input_file, "L");

  BoxParam box(N, L);
  auto potential = read_vector<double>(input_file, "potential");
  return std::make_tuple(box, potential);
}

void run(YAML::Node const &config)
{
  auto [box, potential] = (config["initial-conditions"]["file"] ?
      load_ic(config) : random_ic(config));

  // ~\~ begin <<adhesion_example.md|workflow>>[init]

  // ~\~ end
  // ~\~ begin <<adhesion_example.md|workflow>>[1]
  std::string output_filename
    = config["output"]["hdf5"].as<std::string>();
  fs::create_directories(
    fs::path(output_filename).parent_path());
  H5::H5File output_file(output_filename, H5F_ACC_TRUNC);

  write_attribute(output_file, "N", box.N);
  write_attribute(output_file, "L", box.L);
  write_vector_with_shape(
    output_file, "potential", potential, box.shape());
  // ~\~ end
  // ~\~ begin <<adhesion_example.md|workflow>>[2]
  std::vector<std::unique_ptr<Surface<Point>>> mesh_shape;
  Point centre(box.L/2, box.L/2, box.L/2);
  Vector dz(box.L/15, 0.0, 0.0);
  mesh_shape.emplace_back(new Sphere<K>(centre, 0.4 * box.L));
  mesh_shape.emplace_back(new Plane<K>(centre + dz, dz));
  mesh_shape.emplace_back(new Plane<K>(centre - dz, -dz));
  // ~\~ end
  // ~\~ begin <<adhesion_example.md|workflow>>[3]
  double threshold
    = config["output"]["threshold"].as<double>(0.0);
  std::string walls_filename
    = config["output"]["walls"].as<std::string>();
  std::string filaments_filename
    = config["output"]["filaments"].as<std::string>();
  // ~\~ end
  // ~\~ begin <<adhesion_example.md|workflow>>[4]
  auto time = config["run"]["time"]
    .as<std::vector<double>>();

  unsigned iteration = 0;
  for (double t : time) {
    // ~\~ begin <<adhesion_example.md|workflow-adhesion>>[init]
    std::clog << "Computing regular triangulation for t = "
              << t << " ...\n";
    Adhesion adhesion(box, potential, t);
    // ~\~ end

    // ~\~ begin <<adhesion_example.md|workflow-create-hdf5-group>>[init]
    auto h5_group = output_file.createGroup(
      fmt::format("{}", iteration));
    write_attribute(h5_group, "time", t);
    // ~\~ end
    // ~\~ begin <<adhesion_example.md|workflow-write-nodes>>[init]
    auto nodes = adhesion.get_nodes(threshold);
    write_vector(h5_group, "nodes", nodes);
    // ~\~ end
    // ~\~ begin <<adhesion_example.md|workflow-write-obj>>[init]
    {
      auto walls = adhesion.get_walls(threshold);
      std::string filename = fmt::format(
        walls_filename, fmt::arg("time", t));
      std::clog << "writing to " << filename << "\n";

      std::ofstream ff(filename);
      auto selected_faces = select_mesh(walls, mesh_shape);
      write_faces_to_obj(ff, selected_faces);
      write_mesh(h5_group, "faces", walls);
    }
    // ~\~ end
    // ~\~ begin <<adhesion_example.md|workflow-write-obj>>[1]
    {
      auto filaments = adhesion.get_filaments(threshold);
      std::string filename = fmt::format(
        filaments_filename, fmt::arg("time", t));
      std::clog << "writing to " << filename << "\n";

      std::ofstream ff(filename);
      auto selected_edges = select_mesh(filaments, mesh_shape, false);
      write_edges_to_obj(ff, selected_edges);
      write_mesh(h5_group, "edges", filaments);
    }
    // ~\~ end
    ++iteration;
  }
  // ~\~ end
}
// ~\~ end
