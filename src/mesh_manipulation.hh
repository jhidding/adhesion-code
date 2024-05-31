// ~/~ begin <<appendix.md#src/mesh_manipulation.hh>>[init]
#pragma once
#include <tuple>
#include "mesh.hh"
#include "surface.hh"

// ~/~ begin <<appendix.md#split-polygon>>[init]
template <typename Point>
using Polygon = std::tuple<
    std::vector<Point> *,
    std::vector<unsigned>>;
// ~/~ end
// ~/~ begin <<appendix.md#split-polygon>>[1]
template <typename Point>
using PolygonPair = std::tuple<
    std::vector<Point> *,
    std::vector<unsigned>,
    std::vector<unsigned>>;
// ~/~ end
// ~/~ begin <<appendix.md#split-polygon>>[2]
template <typename Point>
PolygonPair<Point> split_polygon(
    Polygon<Point> const &polygon,
    Surface<Point> const &surface,
    bool closed = true)
{
  std::vector<Point> *vertices = std::get<0>(polygon);
  std::vector<unsigned> orig = std::get<1>(polygon), r1, r2;

  auto is_below = [&vertices, &surface] (unsigned i) -> bool {
    return (surface.oriented_side((*vertices)[i]) == -1);
  };

  if (closed)
    orig.push_back(orig.front());
  else
    r1.push_back(orig.front());

  auto i = orig.begin();
  auto j = i; ++j;
  bool below = is_below(*i);

  while (j != orig.end())
  {
    if (below != is_below(*j)) // surface crossed
    {
      if (auto q = surface.intersect(
            (*vertices)[*i], (*vertices)[*j])) {
        r1.push_back(vertices->size());
        r2.push_back(vertices->size());
        vertices->push_back(*q);
        std::swap(r1, r2);
      } else {
        return PolygonPair<Point>();
      }

      below = not below;
    } else {
      r1.push_back(*j);
      ++i; ++j;
    }
  }

  if (below)
    return PolygonPair<Point>(vertices, r1, r2);
  else
    return PolygonPair<Point>(vertices, r2, r1);
}
// ~/~ end
// ~/~ begin <<appendix.md#select-mesh>>[init]
template <typename Point, typename Info>
Mesh<Point, Info> select_mesh(
    Mesh<Point, Info> const &mesh,
    Surface<Point> const &surface,
    bool closed=true)
{
  Mesh<Point, Info> result;
  result.vertices = mesh.vertices;

  unsigned i = 0, j = 0;
  for (unsigned s : mesh.sizes)
  {
    auto x = mesh.data.begin() + j;
    std::vector<unsigned> vs(x, x + s);

    auto below = std::get<1>(split_polygon(
      Polygon<Point>(&result.vertices, vs),
      surface, closed));

    if (below.size() > 0) {
      result.push_back(below, mesh.info[i]);
    }

    ++i;
    j += s;
  }

  return clean(result);
}

template <typename Point, typename Info>
Mesh<Point, Info> select_mesh(
    Mesh<Point, Info> const &mesh,
    std::vector<std::unique_ptr<Surface<Point>>> const &surfaces,
    bool closed=true)
{
  Mesh<Point, Info> m = mesh;
  for (auto const &s : surfaces)
    m = select_mesh(m, *s, closed);
  return m;
}
// ~/~ end
// ~/~ begin <<appendix.md#clean-mesh>>[init]
template <typename Point, typename Info>
Mesh<Point, Info> clean(
    Mesh<Point, Info> const &source)
{
  Mesh<Point, Info> result;
  std::map<unsigned, unsigned> vertex_map;

  for (unsigned i : source.data) {
    if (vertex_map.count(i) == 0) {
      vertex_map[i] = result.vertices.size();
      result.vertices.push_back(source.vertices[i]);
    }

    result.data.push_back(vertex_map[i]);
  }

  result.info = source.info;
  result.sizes = source.sizes;
  return result;
}
// ~/~ end
// ~/~ end