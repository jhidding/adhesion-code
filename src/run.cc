// ~\~ language=C++ filename=src/run.cc
// ~\~ begin <<adhesion_example.md|src/run.cc>>[0]
#include <iostream>
#include <fstream>
#include <exception>
#include <filesystem>
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

void run(YAML::Node const &config)
{
  // ~\~ begin <<adhesion_example.md|workflow>>[0]
  std::clog << "# Using box with parameters:\n"
            << config["box"] << "\n";
  BoxParam box(
    config["box"]["N"].as<int>(),
    config["box"]["L"].as<double>());

  auto potential = generate_initial_potential(
    box, config);
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
    // ~\~ begin <<adhesion_example.md|workflow-adhesion>>[0]
    std::clog << "Computing regular triangulation for t = "
              << t << " ...\n";
    Adhesion adhesion(box, potential, t);
    // ~\~ end

    // ~\~ begin <<adhesion_example.md|workflow-create-hdf5-group>>[0]
    auto h5_group = output_file.createGroup(
      fmt::format("{}", iteration));
    write_attribute(h5_group, "time", t);
    // ~\~ end
    // ~\~ begin <<adhesion_example.md|workflow-write-nodes>>[0]
    auto nodes = adhesion.get_nodes(threshold);
    write_vector(h5_group, "nodes", nodes);
    // ~\~ end
    // ~\~ begin <<adhesion_example.md|workflow-write-obj>>[0]
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
