// ~\~ language=C++ filename=src/adhesion/get_walls.cc
// ~\~ begin <<adhesion_example.md|src/adhesion/get_walls.cc>>[0]
#include "adhesion.hh"
#include "power_diagram.hh"

Mesh<Point, double> Adhesion::get_walls(
    double threshold) const
{
  return power_diagram_faces(rt, threshold);
}

Mesh<Point, double> Adhesion::get_filaments(
    double threshold) const
{
  return power_diagram_edges(rt, threshold);
}
// ~\~ end
