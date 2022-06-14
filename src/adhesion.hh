// ~\~ language=C++ filename=src/adhesion.hh
// ~\~ begin <<adhesion_example.md|src/adhesion.hh>>[init]
#pragma once
#include "boxparam.hh"
#include "cgal_base.hh"
#include "mesh.hh"

#include <memory>
#include <array>
#include <algorithm>
#include <cstdint>

class Adhesion
{
  // ~\~ begin <<adhesion_example.md|adhesion-members>>[init]
  double                      time;
  // ~\~ begin <<adhesion_example.md|tbb-lock-member>>[init]
  #ifdef CGAL_LINKED_WITH_TBB
  RT::Lock_data_structure     lock;
  #endif
  // ~\~ end
  RT                          rt;
  std::vector<Weighted_point> vertices;
  // ~\~ end

public:
  // ~\~ begin <<adhesion_example.md|adhesion-node-type>>[init]
  enum NodeType : uint32_t {
    VOID, KURTOPARABOLIC, WALL, FILAMENT, CLUSTER, UNDEFINED_NODE_TYPE
  };
  // ~\~ end
  // ~\~ begin <<adhesion_example.md|adhesion-node-struct>>[init]
  struct Node {
    std::array<double, 3> position;
    std::array<double, 3> velocity;
    double    mass;
    NodeType  node_type;
  };
  // ~\~ end

  Adhesion(
    BoxParam const &box,
    std::vector<double> const &potential,
    double t);

  // ~\~ begin <<adhesion_example.md|adhesion-methods>>[init]
  int edge_count(RT::Cell_handle h, double threshold) const;
  Vector velocity(RT::Cell_handle c) const;

  Mesh<Point, double> get_walls(double threshold) const;
  Mesh<Point, double> get_filaments(double threshold) const;
  std::vector<Node> get_nodes(double threshold) const;
  // ~\~ end
};
// ~\~ end
