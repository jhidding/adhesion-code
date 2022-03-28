// ~\~ language=C++ filename=src/adhesion/velocity.cc
// ~\~ begin <<adhesion_example.md|src/adhesion/velocity.cc>>[0]
#include "adhesion.hh"
#include <limits>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Kernel_d/Hyperplane_d.h>

typedef CGAL::Cartesian_d<double>    LiftedK;
typedef CGAL::Point_d<LiftedK>       LiftedPoint;
typedef CGAL::Hyperplane_d<LiftedK>  HyperPlane;

// ~\~ begin <<adhesion_example.md|velocity-define-infinity>>[0]
constexpr double infinity
  = std::numeric_limits<double>::infinity();
// ~\~ end

// ~\~ begin <<adhesion_example.md|velocity-lifted-point>>[0]
inline LiftedPoint lifted_point(
    double x, double y, double z, double w)
{
  double p[4] = { x, y, z, w };
  return LiftedPoint(4, p, p + 4);
}
// ~\~ end

Vector Adhesion::velocity(RT::Cell_handle c) const {
  // ~\~ begin <<adhesion_example.md|velocity-implementation>>[0]
  LiftedPoint points[4];

  for (unsigned i = 0; i < 4; ++i)
  {
    Weighted_point wp = rt.point(c, i);
    auto p = wp.point();
    auto w = wp.weight();
    points[i] = lifted_point(p.x(), p.y(), p.z(), w);
  }
  // ~\~ end
  // ~\~ begin <<adhesion_example.md|velocity-implementation>>[1]
  auto guide = lifted_point(0, 0, 0, -infinity);
  HyperPlane h(points, points + 4, guide, CGAL::ON_NEGATIVE_SIDE);
  // ~\~ end
  // ~\~ begin <<adhesion_example.md|velocity-implementation>>[2]
  auto normal = h.orthogonal_vector();
  auto v  = normal / (2 * time * normal[3]);

  return Vector(v[0], v[1], v[2]);
  // ~\~ end
}
// ~\~ end
