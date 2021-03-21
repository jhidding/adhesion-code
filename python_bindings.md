# Python bindings

``` {.cpp file=src/python/test.cc}
#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include "mesh.hh"
#include "cgal_base.hh"
#include "surface.hh"
#include "sphere.hh"
#include "plane.hh"
#include "mesh_manipulation.hh"

PYBIND11_MAKE_OPAQUE(std::vector<Point>);
PYBIND11_MAKE_OPAQUE(std::vector<unsigned>);
PYBIND11_MAKE_OPAQUE(std::vector<double>);

namespace py = pybind11;

int add(int i, int j) {
    return i + j;
}

Mesh<Point, double> select_mesh_py(
    Mesh<Point, double> const &mesh,
    py::list surfaces,
    bool closed)
{
    Mesh<Point, double> m = mesh;
    for (auto s : surfaces) {
        Surface<Point> *s_ptr = s.cast<Surface<Point> *>();
        m = select_mesh(m, *s_ptr, closed);
    }
    return m;
}

PYBIND11_MODULE(adhesion, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    m.def("add", &add, "A function which adds two numbers");

    py::bind_vector<std::vector<unsigned>>(m, "IndexVector");
    py::bind_vector<std::vector<Point>>(m, "Vertices");
    py::bind_vector<std::vector<double>>(m, "InfoVector");

    py::class_<Point>(m, "Point")
        .def(py::init<double, double, double>())
        .def_property_readonly("x", &Point::x)
        .def_property_readonly("y", &Point::y)
        .def_property_readonly("z", &Point::z);

    py::class_<Vector>(m, "Vector")
        .def(py::init<double, double, double>())
        .def_property_readonly("x", &Vector::x)
        .def_property_readonly("y", &Vector::y)
        .def_property_readonly("z", &Vector::z);

    py::class_<Mesh<Point, double>>(m, "Mesh")
        .def(py::init<>())
        .def_readonly("vertices", &Mesh<Point,double>::vertices)
        .def_readonly("data",     &Mesh<Point,double>::data)
        .def_readonly("sizes",    &Mesh<Point,double>::sizes)
        .def_readonly("info",     &Mesh<Point,double>::info);

    py::class_<Surface<Point>> surface(m, "Surface");
    py::class_<Sphere<K>>(m, "Sphere", surface)
        .def(py::init<Point, double>());
    py::class_<Plane<K>>(m, "Plane", surface)
        .def(py::init<Point, Vector>());

    m.def("select_mesh", &select_mesh_py, "Select part of a mesh", py::arg("mesh"), py::arg("surfaces"), py::arg("closed") = true);
}
```

