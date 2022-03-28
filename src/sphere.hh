// ~\~ language=C++ filename=src/sphere.hh
// ~\~ begin <<appendix.md|src/sphere.hh>>[0]
#pragma once
#include "surface.hh"

template <typename K>
class Sphere: public Surface<typename K::Point_3>
{
  using Point   = typename K::Point_3;
  using Vector  = typename K::Vector_3;
  // ~\~ begin <<appendix.md|sphere-members>>[0]
  Point  origin;
  double radius_squared;
  // ~\~ end

public:
  Sphere(Point const &p, double r):
      origin(p), radius_squared(r*r) {}

  // ~\~ begin <<appendix.md|sphere-oriented-side>>[0]
  int oriented_side(Point const &p) const
  {
    double d = (p - origin).squared_length();

    if (d < radius_squared)
      return -1;

    if (d > radius_squared)
      return +1;

    return 0;
  }
  // ~\~ end
  // ~\~ begin <<appendix.md|sphere-intersect>>[0]
  std::optional<Point> intersect(Point const &a, Point const &b) const
  {
    if (oriented_side(a) * oriented_side(b) >= 0)
      return std::nullopt;

    Vector m = b - a, n = a - origin;
    double m_sqr = m.squared_length(),
            n_sqr = n.squared_length(),
            mn    = m*n;

    double D = mn*mn - (m_sqr * (n_sqr - radius_squared));

    if (D < 0)
      return std::nullopt;               // shouldn't happen

    double sol_m = (- mn - sqrt(D)) / m_sqr,
           sol_p = (- mn + sqrt(D)) / m_sqr;

    if ((sol_m >= 0) and (sol_m <= 1.0))
        return a + m*sol_m;

    if ((sol_p >= 0) and (sol_p <= 1.0))
        return a + m*sol_p;

    return std::nullopt;                 // shouldn't happen
  }
  // ~\~ end
};
// ~\~ end
