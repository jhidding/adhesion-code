# Appendix

This appendix provides convenience interfaces to the FFTW3, GSL integration and HDF5 libraries.

## Numerical integration

The GNU Scientific Library (GSL) is a collection of numerical algorithms written in C. In many cases it is convenient to use a C++ wrapper in stead of the C API directly. Here we define an interface around the `gsl_integration_qagiu` function, which integrates functions from a lower bound to upper infinity. The C API expects a pointer to a function and a `gsl_integration_workspace` to be provided. The function pointed to should accept a `double` and a `void *`. The latter argument can be used to pass around extra arguments of which the function knows the intended type. In our case we reinterpret the void pointer as a generic C++ `std::function` pointer and then call that. This allows any callable C++ function to be passed to `integrate_qagiu`.

Since we only call the integrator once, we do allocation and destruction of the `gsl_integration_workspace` in the body of `integrate_qagiu`. In a situation where `integrate_qagiu` would be called repeatedly, it is better to use a more complicated wrapper.

``` {.cpp #gsl-integrate-qagiu}
#include <gsl/gsl_integration.h>

using Function = std::function<double (double)>;

double integration_helper(double x, void *params)
{
  auto f = reinterpret_cast<Function *>(params);
  return (*f)(x);
}

template <typename F>
double integrate_qagiu(
    F const &func, double lower,
    double epsabs, double epsrel, size_t limit=1024)
{
  double x, abserr;
  Function integrant = func;
  gsl_integration_workspace *workspace =
    gsl_integration_workspace_alloc(limit);

  gsl_function f;
  f.function = &integration_helper;
  f.params = reinterpret_cast<void *>(&integrant);

  gsl_integration_qagiu(
    &f, lower, epsabs, epsrel, limit,
    workspace, &x, &abserr);

  gsl_integration_workspace_free(workspace);
  return x;
}
```

## Manipulating meshes

To visualise a structure it is important to limit the visualisation to a specific region. Otherwise the image is flooded with too many polygons and we lose the aim of visualisation: making structures visible.

In the `mesh_manipulation.hh` header file we define methods to cut a mesh using a spherical surface or a plane.


``` {.cpp file=src/mesh_manipulation.hh}
#pragma once
#include <tuple>
#include "mesh.hh"
#include "surface.hh"

<<split-polygon>>
<<select-mesh>>
<<clean-mesh>>
```

### Surfaces

To select parts of a mesh we need to define a surface that can tell us on what side a point lies, and if we have two points, if and where the segment between those points intersects. This concept of a `Surface` is embodied by the following abstract base class:

``` {.cpp file=src/surface.hh}
#pragma once
<<include-optional>>

template <typename Point>
class Surface
{
public:
  virtual int oriented_side(Point const &p) const = 0;
  virtual std::optional<Point> intersect(
    Point const &a, Point const &b) const = 0;

  virtual ~Surface() {}
};
```

#### `std::optional`

We're using a library feature of C++17, namely `std::optional`. In C++14, this is included in the `std::experimental` namespace. With this macro, we can use `std::optional` in both cases.

``` {.cpp #include-optional}
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
```

### Split a polygon

Given an implementation of such a class, we  can implement a function that will split a polygon in two parts, each on either side of the surface. This makes certain assumptions about the surface and the polygon that will not always hold, but suffice for our purposes of visualisation.

We still have to define the input `Polygon` type,

``` {.cpp #split-polygon}
template <typename Point>
using Polygon = std::tuple<
    std::vector<Point> *,
    std::vector<unsigned>>;
```

The algorithm returns two polygons, one of which lies below the surface, the other above.

``` {.cpp #split-polygon}
template <typename Point>
using PolygonPair = std::tuple<
    std::vector<Point> *,
    std::vector<unsigned>,
    std::vector<unsigned>>;
```

The third argument `closed` tells if the given path is cyclic or not. This is a distinction between having a polygon or a chain of line segments.

``` {.cpp #split-polygon}
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
```

### The sphere

Unfortunately CGAL seems to have no function to compute the intersection of a `Sphere` with a `Segment`. We'll take the opportunity to do a bit computational geometry our selves.

``` {.cpp file=src/sphere.hh}
#pragma once
#include "surface.hh"

template <typename K>
class Sphere: public Surface<typename K::Point_3>
{
  using Point   = typename K::Point_3;
  using Vector  = typename K::Vector_3;
  <<sphere-members>>

public:
  Sphere(Point const &p, double r):
      origin(p), radius_squared(r*r) {}

  <<sphere-oriented-side>>
  <<sphere-intersect>>
};
```

The `Sphere` data type will have two members:

``` {.cpp #sphere-members}
Point  origin;
double radius_squared;
```

Now, any point will have an orientation with respect to this sphere,

$${\rm orient}(\vec{o}, r, \vec{p}) = \begin{cases}
-1 & \quad {\rm if}~(\vec{p} - \vec{o})^2 < r^2\\
0 & \quad {\rm if}~(\vec{p} - \vec{o})^2 = r^2\\
1 & \quad {\rm if}~(\vec{p} - \vec{o})^2 > r^2\\
\end{cases}$$

``` {.cpp #sphere-oriented-side}
int oriented_side(Point const &p) const
{
  double d = (p - origin).squared_length();

  if (d < radius_squared)
    return -1;

  if (d > radius_squared)
    return +1;

  return 0;
}
```

Next, given a directed segment from $\vec{a}$ to $\vec{b}$, we'd like to know if and where the segment will intersect the sphere first. This is done by solving a quadratic equation in terms of the vector going from $\vec{a}$ to $\vec{b}$, $\vec{m} = \vec{b} - \vec{a}$, and the vector from the origin $\vec{o}$ to $\vec{a}$, $\vec{n} = \vec{a} - \vec{o}$.

The points on the segment are given by the function
$$s(t) = \vec{a} + t\vec{m}.$$
And the sphere is defined by the equation
$$(\vec{x} - \vec{o})^2 = r^2.$$
Then the resulting equation is,
$$s^2 = (\vec{a} + t\vec{m} - \vec{o})^2 = n^2 + 2t \vec{n} \vec{m} + m^2 t^2,$$
which is solved by,
$$t = \vec{m}\cdot\vec{n} \pm \sqrt{m^2 (n^2 - r^2)}.$$

If the discriminant $D = m^2 (n^2 - r^2)$ is negative, there is no intersection. We reflect this by returning a `std::optional<Point>`. Also, we only return a point if the found solution is for $0 \le t \le 1$.

``` {.cpp #sphere-intersect}
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
```

### Plane

Similar to the definition of a sphere, we'll define a plane.

``` {.cpp file=src/plane.hh}
#pragma once
#include "surface.hh"

<<sign-function>>

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

  <<plane-oriented-side>>
  <<plane-intersect>>
};
```

#### `sign` function

This generic function returns -1, 0 or 1 according to the sign of the given number.

``` {.cpp #sign-function}
template <typename T>
inline int sign(T a) {
  return (a < 0 ? -1 : (a == 0 ? 0 : 1));
}
```

#### Plane Oriented side

The oriented side of a point ${\bf a}$ to a plane with normal ${\bf n}$ through point ${\bf c}$ is given by

$${\rm sign}(({\bf a} - {\bf c}) \cdot {\bf n}).$$

A point on the positive side lies in the direction of the normal vector.

``` {.cpp #plane-oriented-side}
int oriented_side(Point const &a) const {
  return sign((a - centre) * normal);
}
```

#### Plane Intersection

Given two points ${\bf a}$ and ${\bf b}$, we have a vector ${\bf v} = {\bf b} - {\bf a}, and a vector ${\bf u} = {\bf c} - {\bf a}$ to the point in the plane. Then we can define a time $t$ to parametrise the line segment from ${\bf a}$ at $t = 0$ to ${\bf b}$ at $t = 1$. The intersection with the plane is at

$$t = \frac{{\bf u} \cdot {\bf n}}{{\bf v} \cdot {\bf n}}.$$

``` {.cpp #plane-intersect}
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
```

### Clean mesh

Once we have selected a part of the mesh by cutting polygons, the mesh contains a lot of vertices that have no associated polygon. These vertices need to be removed. This also means that the vertex indices contained in the logical polygon data need to be remapped. The `clean` function takes a `Mesh`, only copies the vertices that are referenced in the polygons and remaps the polygons, creating a new `Mesh` instance.

``` {.cpp #clean-mesh}
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
```

### Select a mesh

``` {.cpp #select-mesh}
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
```

## Writing to disk

For reading and writing HDF5, there is a C++ interface available. This interface is still too low level for our purposes.

### Interface

``` {.cpp file=src/writers.hh}
#pragma once
#include "cgal_base.hh"
#include "adhesion.hh"
#include "mesh.hh"
#include <H5Cpp.h>

<<hdf5-file-or-group>>

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

template <typename T>
std::vector<T> read_vector(
   FileOrGroup &group,
   std::string const &name)
{
  auto data_set = group.openDataSet(name);
  auto data_space = data_set.getSpace();
  auto size = data_space.getSimpleExtentNpoints();
  std::clog << "Reading " << size << " elements from " << name << std::endl;
  std::vector<T> data(size);
  data_set.read(data.data(), H5TypeFactory<T>::get());
  return data;
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

template <typename T, typename Group>
T read_attribute(
    Group &group,
    std::string const &name)
{
  T value;
  auto attr = group.openAttribute(name);
  attr.read(H5TypeFactory<T>::get(), &value);
  return value;
}

extern void write_mesh(
    FileOrGroup &group,
    std::string const &name,
    Mesh<Point, double> const &mesh);
```

There has been a very recent change in the API of HDF5. So we define a proxy type that is common to `H5::Group` and `H5::H5File`.

``` {.cpp #hdf5-file-or-group}
#if H5_VERSION_GE(1, 10, 1)
using FileOrGroup = H5::Group;
#else
using FileOrGroup = H5::CommonFG;
#endif
```

### Saving to HDF5

``` {.cpp file=src/writers/h5_node_type.ih}
extern H5::CompType h5_node_type();

template <>
struct H5TypeFactory<Adhesion::Node>
{
  static H5::CompType get()
  {
    return h5_node_type();
  }
};
```

``` {.cpp file=src/writers/h5_node_type.cc}
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
```

### Writing meshes to HDF5

``` {.cpp file=src/writers/hdf5-mesh.cc}
#include "writers.hh"

void write_mesh(
    FileOrGroup &file,
    std::string const &name,
    Mesh<Point, double> const &mesh)
{
  auto group = file.createGroup(name);

  std::vector<double> vertex_data;
  for (Point const &v : mesh.vertices) {
    vertex_data.push_back(v[0]);
    vertex_data.push_back(v[1]);
    vertex_data.push_back(v[2]);
  }

  std::vector<hsize_t> vertex_shape { mesh.vertices.size(), 3 };
  write_vector_with_shape(group, "vertices", vertex_data, vertex_shape);
  write_vector(group, "info",  mesh.info);
  write_vector(group, "sizes", mesh.sizes);
  write_vector(group, "data",  mesh.data);
}
```

### Writing an OBJ file {#obj-file-format}

``` {.cpp file=src/write_obj.hh}
#pragma once
#include <iostream>

#include "mesh.hh"

template <typename Point>
void write_to_obj(
    std::ostream &out, std::string const &pre,
    Mesh<Point, double> const &mesh)
{
  for (Point const &p : mesh.vertices) {
    out << "v " << p << " 1.0\n";
  }
  out << "\n";

  for (double a : mesh.info) {
    out << "vt " << a << " 0\n";
  }
  out << "\n";

  unsigned i = 1, j = 0;
  for (unsigned s : mesh.sizes) {
    out << pre;
    for (unsigned k = 0; k < s; ++k) {
      out << " " << mesh.data[j] + 1 << "/" << i;
      ++j;
    }
    out << "\n";
    ++i;
  }
}

template <typename Point>
void write_edges_to_obj(
    std::ostream &out,
    Mesh<Point, double> const &mesh)
{
  write_to_obj(out, "l", mesh);
}

template <typename Point>
void write_faces_to_obj(std::ostream &out, Mesh<Point, double> const &mesh)
{
  write_to_obj(out, "f", mesh);
}
```

## Fourier interface

``` {.cpp file=src/fft.hh}
#pragma once
#include <fftw3.h>
#include <memory>
#include <complex>
#include <vector>

#include "boxparam.hh"
#include <iostream>

template <typename T>
class FFTW_allocator: public std::allocator<T>
{
public:
  typedef T           value_type;
  typedef T *         pointer;
  typedef T &         reference;
  typedef T const *   const_pointer;
  typedef T const &   const_reference;
  typedef size_t      size_type;
  typedef ptrdiff_t   difference_type;

  pointer allocate(
      size_t n,
      std::allocator<void>::const_pointer hint = 0)
  {
    return reinterpret_cast<T *>(fftw_malloc(n * sizeof(T)));
  }

  void deallocate(
      pointer p,
      size_t n)
  {
    fftw_free(p);
  }
};

class RFFT3
{
public:
  using c64 = std::complex<double>;
  std::vector<c64, FFTW_allocator<c64>>
    fourier_space;
  std::vector<double, FFTW_allocator<double>>
    real_space;

private:
  BoxParam    box;
  fftw_plan   d_plan_fwd, d_plan_bwd;

public:
  RFFT3(BoxParam const &box_):
    fourier_space(box_.rfft_size()),
    real_space(box_.size()),
    box(box_)
  {
    int N = static_cast<int>(box.N);

    d_plan_fwd = fftw_plan_dft_r2c_3d(N, N, N,
      reinterpret_cast<double *>(real_space.data()),
      reinterpret_cast<fftw_complex *>(fourier_space.data()),
      FFTW_ESTIMATE);

    d_plan_bwd = fftw_plan_dft_c2r_3d(N, N, N,
      reinterpret_cast<fftw_complex *>(fourier_space.data()),
      reinterpret_cast<double *>(real_space.data()),
      FFTW_ESTIMATE);
  }

  void forward_transform()
  {
      fftw_execute(d_plan_fwd);
  }

  void backward_transform()
  {
      fftw_execute(d_plan_bwd);
      for (double &z : real_space) z /= box.size();
  }
};
```

\pagebreak
