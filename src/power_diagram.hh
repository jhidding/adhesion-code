// ~/~ begin <<adhesion_example.md#src/power_diagram.hh>>[init]
#include "cgal_base.hh"
#include "mesh.hh"

extern Mesh<Point, double> power_diagram_faces(
  RT const &rt, double threshold);

extern Mesh<Point, double> power_diagram_edges(
  RT const &rt, double threshold);
// ~/~ end