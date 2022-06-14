// ~\~ language=C++ filename=src/cgal_base.hh
// ~\~ begin <<adhesion_example.md|src/cgal_base.hh>>[init]
#pragma once
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_3.h>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;

// ~\~ begin <<adhesion_example.md|regular-triangulation-type>>[init]
#ifdef CGAL_LINKED_WITH_TBB
  using TDS = CGAL::Triangulation_data_structure_3<
      CGAL::Regular_triangulation_vertex_base_3<K>,
      CGAL::Regular_triangulation_cell_base_3<K>,
      CGAL::Parallel_tag>;

  using RT = CGAL::Regular_triangulation_3<K, TDS>;
#else
  using RT = CGAL::Regular_triangulation_3<K>;
#endif
// ~\~ end

typedef K::Vector_3           Vector;
typedef RT::Bare_point        Point;
typedef RT::Edge              Edge;
typedef RT::Weighted_point    Weighted_point;
typedef RT::Segment           Segment;
typedef RT::Tetrahedron       Tetrahedron;
// ~\~ end
