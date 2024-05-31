// ~/~ begin <<appendix.md#src/writers/h5_node_type.cc>>[init]
#include "adhesion.hh"
#include "writers.hh"
#include <cstddef>

H5::CompType h5_node_type()
{
  hsize_t dim = 3;

  H5::FloatType ft(H5::PredType::NATIVE_DOUBLE);
  H5::IntType   it(H5::PredType::NATIVE_UINT32);
  H5::ArrayType at(ft, 1, &dim);
  H5::CompType  ct(sizeof(Adhesion::Node));

  ct.insertMember("position",  0, at);
  ct.insertMember("velocity",  offsetof(Adhesion::Node, velocity), at);
  ct.insertMember("mass",      offsetof(Adhesion::Node, mass), ft);
  ct.insertMember("node_type", offsetof(Adhesion::Node, node_type), it);
  return ct;
}
// ~/~ end