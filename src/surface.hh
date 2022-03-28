// ~\~ language=C++ filename=src/surface.hh
// ~\~ begin <<appendix.md|src/surface.hh>>[0]
#pragma once
// ~\~ begin <<appendix.md|include-optional>>[0]
#if __cplusplus == 201402L
#include <experimental/optional>
namespace std {
using namespace std::experimental;
}
#elif __cplusplus == 201703L
#include <optional>
#else
#error "Unrecognised C++ version."
#endif
// ~\~ end

template <typename Point>
class Surface
{
public:
  virtual int oriented_side(Point const &p) const = 0;
  virtual std::optional<Point> intersect(
    Point const &a, Point const &b) const = 0;

  virtual ~Surface() {}
};
// ~\~ end
