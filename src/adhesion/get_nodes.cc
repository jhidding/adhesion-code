// ~/~ begin <<adhesion_example.md#src/adhesion/get_nodes.cc>>[init]
#include "adhesion.hh"

// ~/~ begin <<adhesion_example.md#adhesion-type-from-edge-count>>[init]
inline Adhesion::NodeType type_from_edge_count(int n)
{
  switch (n) {
    case 0: return Adhesion::VOID;
    case 1:
    case 2: return Adhesion::KURTOPARABOLIC;
    case 3:
    case 4: return Adhesion::WALL;
    case 5: return Adhesion::FILAMENT;
    case 6: return Adhesion::CLUSTER;
  }
  return Adhesion::UNDEFINED_NODE_TYPE;
}
// ~/~ end

std::vector<Adhesion::Node>
Adhesion::get_nodes(double threshold) const
{
  std::vector<Adhesion::Node> result;

  for (auto c  = rt.finite_cells_begin();
            c != rt.finite_cells_end();
            ++c) {
    std::array<double, 3>  ps, vs;
    Point    p = rt.dual(c);
    Vector   v = velocity(c);
    for (unsigned i = 0; i < 3; ++i) {
      ps[i] = p[i];
      vs[i] = v[i];
    }
    NodeType n = type_from_edge_count(
                   edge_count(c, threshold));
    double   m = rt.tetrahedron(c).volume();

    result.push_back(Adhesion::Node({ps, vs, m, n}));
  }

  return result;
}
// ~/~ end