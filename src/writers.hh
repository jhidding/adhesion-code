// ~\~ language=C++ filename=src/writers.hh
// ~\~ begin <<appendix.md|src/writers.hh>>[0]
#pragma once
#include "cgal_base.hh"
#include "adhesion.hh"
#include "mesh.hh"
#include <H5Cpp.h>

// ~\~ begin <<appendix.md|hdf5-file-or-group>>[0]
#if H5_VERSION_GE(1, 10, 1)
using FileOrGroup = H5::Group;
#else
using FileOrGroup = H5::CommonFG;
#endif
// ~\~ end

template <typename T>
struct H5TypeFactory {};

template <>
struct H5TypeFactory<double>
{
  static H5::FloatType get()
  {
    return H5::FloatType(H5::PredType::NATIVE_DOUBLE);
  }
};

template <>
struct H5TypeFactory<unsigned>
{
  static H5::IntType get()
  {
    return H5::IntType(H5::PredType::NATIVE_UINT);
  }
};

#include "writers/h5_node_type.ih"

template <typename V, typename S>
void write_vector_with_shape(
    FileOrGroup &group,
    std::string const &name,
    V const &v,
    S const &shape)
{
  std::vector<hsize_t> hshape(shape.begin(), shape.end());
  auto          data_type = H5TypeFactory<typename V::value_type>::get();
  H5::DataSpace data_space(hshape.size(), hshape.data());
  auto          data_set = group.createDataSet(name, data_type, data_space);
  data_set.write(v.data(), data_type);
}

template <typename V>
void write_vector(
    FileOrGroup &group,
    std::string const &name,
    V const &v)
{
  std::vector<hsize_t> shape { v.size() };
  write_vector_with_shape(group, name, v, shape);
}

template <typename Group, typename T>
void write_attribute(
    Group &group,
    std::string const &name,
    T const &value)
{
  auto attr = group.createAttribute(
    name, H5TypeFactory<T>::get(), H5::DataSpace());
  attr.write(H5TypeFactory<T>::get(), &value);
}

extern void write_mesh(
    FileOrGroup &group,
    std::string const &name,
    Mesh<Point, double> const &mesh);
// ~\~ end
