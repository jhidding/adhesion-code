// ~/~ begin <<appendix.md#src/write_obj.hh>>[init]
#pragma once
#include <iostream>

#include "mesh.hh"

template <typename Point>
void write_to_obj(
    std::ostream &out, std::string const &pre,
    Mesh<Point, double> const &mesh)
{
  for (Point const &p : mesh.vertices) {
    out << "v " << p << " 1.0\n";
  }
  out << "\n";

  for (double a : mesh.info) {
    out << "vt " << a << " 0\n";
  }
  out << "\n";

  unsigned i = 1, j = 0;
  for (unsigned s : mesh.sizes) {
    out << pre;
    for (unsigned k = 0; k < s; ++k) {
      out << " " << mesh.data[j] + 1 << "/" << i;
      ++j;
    }
    out << "\n";
    ++i;
  }
}

template <typename Point>
void write_edges_to_obj(
    std::ostream &out,
    Mesh<Point, double> const &mesh)
{
  write_to_obj(out, "l", mesh);
}

template <typename Point>
void write_faces_to_obj(std::ostream &out, Mesh<Point, double> const &mesh)
{
  write_to_obj(out, "f", mesh);
}
// ~/~ end