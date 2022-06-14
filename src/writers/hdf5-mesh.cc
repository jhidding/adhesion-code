// ~\~ language=C++ filename=src/writers/hdf5-mesh.cc
// ~\~ begin <<appendix.md|src/writers/hdf5-mesh.cc>>[init]
#include "writers.hh"

void write_mesh(
    FileOrGroup &file,
    std::string const &name,
    Mesh<Point, double> const &mesh)
{
  auto group = file.createGroup(name);

  std::vector<double> vertex_data;
  for (Point const &v : mesh.vertices) {
    vertex_data.push_back(v[0]);
    vertex_data.push_back(v[1]);
    vertex_data.push_back(v[2]);
  }

  std::vector<hsize_t> vertex_shape { mesh.vertices.size(), 3 };
  write_vector_with_shape(group, "vertices", vertex_data, vertex_shape);
  write_vector(group, "info",  mesh.info);
  write_vector(group, "sizes", mesh.sizes);
  write_vector(group, "data",  mesh.data);
}
// ~\~ end
