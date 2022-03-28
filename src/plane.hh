// ~\~ language=C++ filename=src/plane.hh
// ~\~ begin <<appendix.md|src/plane.hh>>[0]
#pragma once
#include "surface.hh"

// ~\~ begin <<appendix.md|sign-function>>[0]
template <typename T>
inline int sign(T a) {
  return (a < 0 ? -1 : (a == 0 ? 0 : 1));
}
// ~\~ end

template <typename K>
class Plane: public Surface<typename K::Point_3>
{
  using Point   = typename K::Point_3;
  using Vector  = typename K::Vector_3;

  Point  centre;
  Vector normal;

public:
  Plane(Point const &centre, Vector const &normal)
    : centre(centre)
    , normal(normal)
  {}

  // ~\~ begin <<appendix.md|plane-oriented-side>>[0]
  int oriented_side(Point const &a) const {
    return sign((a - centre) * normal);
  }
  // ~\~ end
  // ~\~ begin <<appendix.md|plane-intersect>>[0]
  std::optional<Point> intersect(Point const &a, Point const &b) const
  {
    Vector u = centre - a,
           v = b - a;

    if (v * normal == 0)
      return std::nullopt;

    double t = (u * normal) / (v * normal);

    if (t < 0 || t > 1)
      return std::nullopt;

    return a + t * v;
  }
  // ~\~ end
};
// ~\~ end
