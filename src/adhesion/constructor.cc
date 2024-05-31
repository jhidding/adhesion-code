// ~/~ begin <<adhesion_example.md#src/adhesion/constructor.cc>>[init]
#include "adhesion.hh"

Adhesion::Adhesion(
    BoxParam const &box,
    std::vector<double> const &potential,
    double t)
  : time(t)
  // ~/~ begin <<adhesion_example.md#tbb-initialise-lock>>[init]
  #ifdef CGAL_LINKED_WITH_TBB
  , lock(CGAL::Bbox_3(
      0.0, 0.0, 0.0,
      box.L, box.L, box.L), box.N)
  , rt(K(), &lock)
  #endif
  // ~/~ end
{
  vertices.reserve(box.size());
  for (size_t i = 0; i < box.size(); ++i)
  {
    vertices.emplace_back(
      box.point<Point>(i),
      potential[i] * 2 * time);
  }

  rt.insert(
    vertices.begin(),
    vertices.end());
}
// ~/~ end