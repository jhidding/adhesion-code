// ~\~ language=C++ filename=src/power_diagram/faces.cc
// ~\~ begin <<adhesion_example.md|src/power_diagram/faces.cc>>[0]
#include "power_diagram.hh"

// ~\~ begin <<adhesion_example.md|pd-is-big-edge>>[0]
inline bool is_big_edge(
    RT const &rt, RT::Edge const &e, double threshold)
{
  double l = rt.segment(e).squared_length();
  return l < threshold;
}
// ~\~ end

Mesh<Point, double> power_diagram_faces(
  RT const &rt,
  double threshold)
{
  Mesh<Point, double> mesh;
  // ~\~ begin <<adhesion_example.md|pd-dual-vertex>>[0]
  std::map<RT::Cell_handle, unsigned> cell_index;

  auto get_dual_vertex = [&rt, &cell_index, &mesh] (
      RT::Cell_handle const &h) -> unsigned
  {
    if (cell_index.count(h) == 0)
    {
      cell_index[h] = mesh.vertices.size();
      mesh.vertices.push_back(rt.dual(h));
    }

    return cell_index[h];
  };
  // ~\~ end
  // ~\~ begin <<adhesion_example.md|pd-walls-loop>>[0]
  for (auto e = rt.finite_edges_begin();
        e != rt.finite_edges_end();
        ++e)
  {
    // ~\~ begin <<adhesion_example.md|pd-edge-check>>[0]
    if (is_big_edge(rt, *e, threshold)) {
      continue;
    }
    // ~\~ end
    // ~\~ begin <<adhesion_example.md|pd-collect-dual>>[0]
    std::vector<unsigned> polygon;
    auto first = rt.incident_cells(*e), c = first;
    bool ok = true;

    do {
      if (rt.is_infinite(++c)) {
          ok = false;
          break;
      }

      polygon.push_back(get_dual_vertex(c));
    } while (c != first);
    // ~\~ end
    // ~\~ begin <<adhesion_example.md|pd-add-to-mesh>>[0]
    if (ok) {
      double l = sqrt(rt.segment(*e).squared_length());
      mesh.push_back(polygon, l);
    }
    // ~\~ end
  }
  // ~\~ end
  return mesh;
}
// ~\~ end
