// ~\~ language=C++ filename=src/power_diagram/edges.cc
// ~\~ begin <<adhesion_example.md|src/power_diagram/edges.cc>>[init]
#include "power_diagram.hh"

// ~\~ begin <<adhesion_example.md|pd-is-big-facet>>[init]
inline bool is_big_facet(
    RT const &rt, RT::Facet const &f, double threshold)
{
  RT::Cell_handle c = f.first;
  unsigned        k = f.second;

  for (unsigned i = 0; i < 4; ++i) {
    if (i == k) {
      continue;
    }

    for (unsigned j = 0; j < i; ++j) {
      if (j == k) {
        continue;
      }

      if (rt.segment(c, i, j).squared_length() < threshold) {
        return false;
      }
    }
  }

  return true;
}
// ~\~ end

Mesh<Point, double> power_diagram_edges(
  RT const &rt,
  double threshold)
{
  Mesh<Point, double> mesh;
  // ~\~ begin <<adhesion_example.md|pd-dual-vertex>>[init]
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

  for (auto f = rt.finite_facets_begin();
       f != rt.finite_facets_end();
       ++f)
  {
    if (!is_big_facet(rt, *f, threshold)) {
      continue;
    }

    double area = sqrt(rt.triangle(*f).squared_area());
    auto mirror = rt.mirror_facet(*f);

    if (rt.is_infinite(f->first) || rt.is_infinite(mirror.first)) {
      continue;
    }

    std::vector<unsigned> polygon;
    polygon.push_back(get_dual_vertex(f->first));
    polygon.push_back(get_dual_vertex(mirror.first));
    mesh.push_back(polygon, area);
  }

  return mesh;
}
// ~\~ end
