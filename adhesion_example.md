---
title: Computing the adhesion model using C++ and CGAL
author: Johan Hidding
date: September 1, 2018
---

# Abstract

We present a (relatively) small example of using the CGAL library to run the adhesion model. A basic CGAL program can seem intimidating to start with. The CGAL manual provides basic examples that offer a good starting point. In our case we are interested in creating regular triangulations. We adapted the example from the manual.

## Version

``` {.cpp #version}
#define VERSION "0.1"
```

# Introduction

This document is aimed to be self-containing. This means that all the code to build a working adhesion model is included. A significant fraction of the code in here is dealing with generating initial conditions for the actual model.

We've tried to limit the involvement of too much boilerplate code by using existing libraries where possible.

If you are only interested only in the CGAL parts, it should be safe to skip the section on generating cosmological initial conditions.

## Literate programming

This example is written in a style of *literate programming*. This document contains a complete and functioning example of working with CGAL to compute the adhesion model. For didactic reasons we don't always give the listing of an entire source file in one go. In stead, we use a system of references known as *noweb*.

Inside source fragments you may encounter a line like,

``` {.cpp}
#include <cmath>

std::pair<float, float> solve_quadratic(
    float a, float b, float c)
{
  <<abc-formula>>
}
```

which is elsewhere specified as:

``` {.cpp #abc-formula}
float d = sqrt(b*b - 4*a*c);
```

A definition can be appended with more code as follows:

``` {.cpp #abc-formula append=true}
return std::make_pair(
    (-b - d)/(2*a),
    (-b + d)/(2*a));
```

## Prerequisites

- C++ 14 compiler
- CGAL 4.12 - The Computational Geometry Algorithm Library
- XTensor 0.17 - A template library for handling array data in C++
- FFTW3 3.3 - The Fastest Fourier Transform in the West
- XTensor-FFTW - Tying FFTW3 to an XTensor interface
- hdf5-cpp - [HDF5](https://support.hdfgroup.org/HDF5/doc/cpplus_RM/index.html) is used to store large blobs of binary data and meta data.
- yaml-cpp - [YAML-cpp](https://github.com/jbeder/yaml-cpp) is a YAML parser for C++. We use it to parse configuration files.
- argagg - [ArgAgg](https://github.com/vietjtnguyen/argagg) stands for Argument Aggregator and is a C++ command-line argument parser.
- fmt - [fmt](http://fmtlib.net/latest/index.html) is a string formatting library that has a similar interface as Python's.

All of these packages are available in the Debian GNU/Linux package repositories.

## CGAL Geometry kernels

CGAL comes with a set of *geometry kernels*. Each kernel bundles basic type definitions like `Point`, `Vector`, `Circle`, etc. and geometric operations on those types. Depending on the requirements of the programmer, we can choose different implementations of these concepts. These implementations vary in representation of real numbers, vector quantities, and how geometric operations are computed on them. Some abstract applications require an exact representation of numbers while other more cosmological applications can afford to be more liberal with regards to exactness.

The algorithms that actually do the advanced geometric computations, like the regular triangulations we use, are implemented in generic terms using C++ template techniques. We need to supply those algorithms the correct geometry kernel for our application. This is why all CGAL programs start with a list of template type definitions.

In our case, what we need is a double precision floating point representation of numbers, while retaining logical consistency in geometric predicates, of which the co-linearity test is the most obvious example. This is provided by the `Exact_predicates_inexact_constructions_kernel` kernel.

We collect those type definitions in a separate header file:

``` {.cpp file=src/cgal_base.hh}
#pragma once
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
// #include <CGAL/Periodic_3_regular_triangulation_traits_3.h>
// #include <CGAL/Periodic_3_regular_triangulation_3.h>

#include <CGAL/Regular_triangulation_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

// typedef CGAL::Periodic_3_regular_triangulation_traits_3<K>  Gt;
// typedef CGAL::Periodic_3_regular_triangulation_3<Gt>        RT;
typedef CGAL::Regular_triangulation_3<K>               RT;

typedef K::Vector_3           Vector;
typedef K::Point_3            Point;
// typedef RT::Bare_point        Point;
typedef RT::Edge              Edge;
// typedef RT::Iso_cuboid        Iso_cuboid;
typedef RT::Weighted_point    Weighted_point;
typedef RT::Segment           Segment;
typedef RT::Tetrahedron       Tetrahedron;
```

Since we'll be using bare (weightless) points, weighted points, and vectors we defined aliases for those types. Note that CGAL is particular about the difference between points and vectors. Points are locations without absolute properties, whereas vectors describe how to get from one point to the other. Internally they may have the same numerical representation, but this may not strictly be the case for all geometry kernels.

# Initial conditions

The initial conditions are randomly generated on a grid. We suppose a platonic ideal Gaussian random field that underlies our realisation. This is a function that is only defined in probabalistic terms. In cosmology it comes natural that these probabalities do not depend on location. For example, in the case of completely uncorrelated Gaussian white noise, we can ask: what is the probability that this function attains a certain value,

$$P(f(x) = y) = \frac{1}{\sqrt{2\pi \sigma^2}} \exp \left(-\frac{(y - \mu)^2}{2 \sigma^2}\right).$$

We're looking at quantities, like the density perturbation, that have mean $\mu = 0$. When we generate white noise, we're sampling a realisation of such a function $f$ at a limited set of points. This should be considered in contrast with seeing a realisation as as integral quantities to a grid cell. Any integral of a white noise over a finite area results exactly in the mean value.

To get an instance of a physically meaningful field, with non-zero integrals, requires that the values of the function $f$ are positively correlated at the small scale. Taking any two positions $x_1$ and $x_2$, their correlation is

$$\xi(x_1, x_2) = \langle f(x_1) f(x_2) \rangle.$$

Often we write the correlation function $\xi(r)$ because our fields are isotropic and homogeneous. We can now ask the next question: what is the probability that the function $f$ at position $\vec{x}$ attains the value $f(\vec{x}) = y_1$ and at position $\vec{x} + \vec{r}$ attains the value $f(\vec{x} + \vec{r}) = y_2$,

$$P(f(\vec{x}) = y_1, f(\vec{x} + \vec{r}) = y_2) = \frac{1}{\sqrt{2\pi |\Sigma(r)|}} \exp \left(-\frac{1}{2} \begin{pmatrix}y_1\\y_2\end{pmatrix}^T \Sigma^{-1}(r) \begin{pmatrix}y_1\\y_2\end{pmatrix}\right).$$

Here, $\Sigma(r)$ is the corellation matrix,

$$\Sigma(r) = \begin{pmatrix}\sigma^2 & \xi(r)\\\xi(r) & \sigma^2\end{pmatrix},$$

``` {.cpp file=src/initial_conditions.hh}
#pragma once
#include "boxparam.hh"

#include <memory>
#include <xtensor/xtensor.hpp>
#include <yaml-cpp/yaml.h>

extern std::unique_ptr<xt::xtensor<double, 3>>
generate_white_noise(
    BoxParam const &box,
    unsigned long seed);

<<power-spectra>>

extern void compute_potential(
    BoxParam const &box,
    xt::xtensor<double, 3> &white_noise,
    PowerSpectrum const &P);
```

## The simulation box

Next, we need to define the box that we will use. We collect the required parameters, box size in pixels and the physical length, in a structure called `BoxParam`.

``` {.cpp file=src/boxparam.hh}
#pragma once
#include <cstdlib>
#include <cmath>
#include <array>

<<increment-index>>

template <typename T>
inline T sqr(T x) { return x*x; }

struct BoxParam {
  unsigned N;         // number of grid points in a row
  size_t   size;      // number of points in the entire box
  double   L;         // physical size of the box
  double   res;       // resolution of the box
  double   res_sqr;   // square of the resolution

  BoxParam(unsigned N_, double L_):
    N(N_),
    size(N*N*N),
    L(L_),
    res(L/N),
    res_sqr(res*res)
  {}

  std::array<size_t, 3> shape() const
  {
    return std::array<size_t, 3>{ N, N, N };
  }

  <<fourier-properties>>

  <<boxparam-methods>>
};
```

We add a method to compute the $n$-th point on a grid.

``` {.cpp #boxparam-methods}
template <typename Point>
Point point(size_t i) const
{
  int x = i % N;
  int y = (i / N) % N;
  int z = i / (N * N);

  return Point(x * res, y * res, z * res);
}
```

Note that we set the $x$-coordinate to be the fastest changing coordinate in the flattened array. This is known as *row-major* ordering, which the same as how indexing into C/C++ and Python/NumPy arrays works. Later on, we will be adding more methods to the `BoxParam` structure.

### Fourier properties

These functions we'll need when we compute Fourier transforms. The real FFT algorithm saves precious memory by using only half the space of the complex FFT. With the exception of the Nyquist frequencies that makes $N/2 + 1$ for the last axis.

``` {.cpp #fourier-properties append=true}
std::array<size_t, 3> rfft_shape() const {
  return std::array<size_t, 3>{ N, N, N/2 + 1 };
}

size_t rfft_size() const {
  return N * N * (N / 2 + 1);
}
```

Next up, we need to compute the wave number of the Fourier mode represented at a certain index. With a physical box-size of $L$ and a logical size of $N$, we use the following convention

$$k_i = i \frac{2 \pi}{L},~{\rm for}~i \in [0, 1, \dots, N/2, -N/2 - 1, \dots, -1].$$

``` {.cpp #fourier-properties append=true}
double wave_number(int i) const {
  return (int(i) > int(N)/2
          ? int(i) - int(N)
          : int(i)) * (2*M_PI)/L;
}

double k_abs(std::array<size_t, 3> const &loc) const {
  double x = 0.0;
  for (size_t i : loc)
    x += sqr(wave_number(i));
  return sqrt(x);
}
```

### Iterating multi-dimensional arrays

We'll be indexing multi-dimensional arrays. To prevent having to write nested for-loops, we use the `increment_index` helper function.

``` {.cpp #increment-index}
template <unsigned R>
inline unsigned increment_index(
    std::array<size_t, R> const &shape,
    std::array<size_t, R> &index)
{
  for (unsigned i = 0; i < R; ++i) {
    unsigned k = R - i - 1;
    if (++index[k] < shape[k])
      return k;

    index[k] = 0;
  }

  return R;
}
```

## White noise

The `white_noise` function fills a newly created array with random values, following a normal distribution with $\sigma = 1$.

``` {.cpp file=src/white_noise.cc}
#include "initial_conditions.hh"
#include <random>

std::unique_ptr<xt::xtensor<double, 3>>
generate_white_noise(
    BoxParam const &box, unsigned long seed)
{
  auto result = std::make_unique<xt::xtensor<double, 3>>(box.shape());

  std::mt19937 random(seed);
  std::normal_distribution<double> normal;

  for (double &value : *result) {
    value = normal(random);
  }

  return std::move(result);
}
```

## Power spectrum

A power spectrum is a function taking in a value of $k$ in units of $h {\rm Mpc}^{-1}$ giving an amplitude.

``` {.cpp #power-spectra}
using PowerSpectrum = std::function<double (double)>;
using Config = YAML::Node;

extern PowerSpectrum EisensteinHu(
    Config const &cosmology);

extern PowerSpectrum normalize_power_spectrum(
    BoxParam const &box,
    PowerSpectrum const &P,
    Config const &cosmology);
```

### Eisenstein-Hu power spectrum

The power spectrum for CDM is given by an almost scale-free spectrum modified by a *transfer function* $T_0$ which embodies post-inflation physics.

$$P(k) = A k^{n_s} T_0^2(k)$$

To compute $T_0$ from the Boltzmann equation there exists a code called CMBfast, but most often people use a fitting function. Eisenstein & Hu (1997) give the following fitting formula for the CDM transfer function:

$$\begin{aligned}
T_0(q) &=& \frac{L_0}{L_0 + C_0 q^2}\\
L_0(q) &=& \log(2 e + 1.8 q)\\
C_0(q) &=& 14.2 + \frac{731}{1 + 62.5 q},
\end{aligned}$$

where $q$ is the wave number re-scaled to find the *knee* in the CDM power spectrum.

$$q = \frac{k}{h {\rm Mpc}^{-1}} \Theta_{2.7}^2 / \Gamma,$$

where $\Theta_{2.7}$ is the temperature of the CMB divided by 2.7 and $\Gamma = \Omega_0 h$.

``` {.cpp file=src/eisenstein-hu.cc}
#include "initial_conditions.hh"

PowerSpectrum EisensteinHu(Config const &cosmology)
{
  double const
    e         = exp(1),
    Theta_CMB = 2.7255/2.7,
    Omega0    = cosmology["Omega0"].as<double>(),
    h         = cosmology["h"].as<double>(),
    ns        = cosmology["ns"].as<double>(),
    A         = 1122670;

  return [=] (double k)
  {
    double q  = k * pow(Theta_CMB, 2)/(Omega0 * h),
           L0 = log(2*e + 1.8*q),
           C0 = 14.2 + 731.0/(1 + 62.5*q),
           T0 = L0 / (L0 + C0 * pow(q, 2));

    return A * pow(k, ns) * pow(T0, 2);
  };
}
```

### Normalisation

## Applying the power spectrum

We now apply the desired power spectrum to the previously generated white noise. This is done by transforming the white noise to the Fourier domain, multiplying it by the square root of the power spectrum, and then transforming back again.

``` {.cpp file=src/apply_power_spectrum.cc}
#include "initial_conditions.hh"
#include "fft.hh"

void compute_potential(
    BoxParam const &box,
    xt::xtensor<double, 3> &field,
    PowerSpectrum const &P)
{
  RFFT3 rfft(box);
  auto f_shape = box.rfft_shape();

  std::copy(field.begin(), field.end(), rfft.real_space.begin());
  rfft.forward_transform();

  std::array<size_t, 3> loc = {0, 0, 1};
  for (size_t i = 1; i < box.rfft_size(); ++i) {
    double k = box.k_abs(loc);
    rfft.fourier_space[i] *= sqrt(P(k)) / (k * k);
    increment_index<3>(f_shape, loc);
  }

  rfft.backward_transform();
  std::copy(rfft.real_space.begin(), rfft.real_space.end(), field.begin());
}
```

# The Adhesion model

We're solving the inviscid Burgers equation,
$$\partial_t \vec{v} + (\vec{v} \cdot \vec{\nabla}) \vec{v} = \nu \nabla^2 \vec{v},$$
in the limit of $\nu \to 0$. @hopf1950 gave the solution to this equation. 
This solution is given by maximising the function
$$G(\vec{q}, \vec{x}, t) = \Phi_0(\vec{q}) - \frac{(\vec{x} - \vec{q})^2}{2t},$$
to obtain the Eulerian velocity potential
$$\Phi(\vec{x}, t) = \max_q G(\vec{q}, \vec{x}, t).$$
This solution can be computed through the power diagram, given a set of points $S \subset \mathbb{R}^n$, and a weight $w_u$ associated with each point $\vec{u} \in S$, the *power cell* is defined as,
$$V_u = \left\{\vec{x} \in \mathbb{R}^n \big| (\vec{x} - \vec{u})^2 + w_u \le (\vec{x} - \vec{v})^2 + w_v\ \forall \vec{v} \in S\right\}.$$
Setting the weights to,
$$w(\vec{q}) = 2 t \Phi_0(\vec{q}),$$
the adhesion model uses regular triangulations to compute structures directly from the initial conditions.

``` {.cpp file=src/adhesion.hh}
#pragma once
#include "boxparam.hh"
#include "cgal_base.hh"
#include "mesh.hh"
#include <memory>

class Adhesion
{
  RT                          rt;
  std::vector<Weighted_point> vertices;
  double                      time;

public:
  <<adhesion-constructor>>

  int edge_count(RT::Cell_handle h, double threshold) const;
  Vector velocity(RT::Cell_handle c) const;

  Mesh<Point, double> get_walls(double threshold) const;
};
```

## Computing the triangulation

We give each point in the grid a weight proportional to the velocity potential,
$$w_i = 2 t \Phi(q_i).$$
Then we insert these weighted points into the triangulation.

``` {.cpp #adhesion-constructor}
template <typename Array>
Adhesion(BoxParam const &box, Array &&potential, double t)
  : // rt(Iso_cuboid(0, 0, 0, box.L, box.L, box.L))
   time(t)
{
  for (size_t i = 0; i < box.size; ++i)
  {
    vertices.emplace_back(
      box.point<Point>(i),
      potential[i] * 2 * time);
  }

  rt.insert(
    vertices.begin(),
    vertices.end());
}
```

## Filtering for structures

When we want to select filaments or clusters we need to count how many edges of a certain cell in the regular triangulation exceeds a given threshold. The function `Adhesion::edge_count` takes a cell handle and a threshold and returns the number of long edges. This count determines if the cell is part of a void, wall, filament or node.

| edge count | structure |
|-----------:|-----------|
|          0 | void      |
|          1 | kurtoparabolic point |
|          2 | wall      |
|          3 | filament  |
|          4 | node      |

``` {.cpp file=src/adhesion_edge_count.cc}
#include "adhesion.hh"

int Adhesion::edge_count(RT::Cell_handle h, double threshold) const
{
  int count = 0;
  for (unsigned i = 1; i < 4; ++i) {
    for (unsigned j = 0; j < i; ++j) {
      // auto segment = rt.periodic_segment(h, i, j);
      // double l = rt.construct_segment(segment)
      //              .squared_length();

      auto segment = rt.segment(h, i, j);
      double l = segment.squared_length();

      if (l > threshold) {
        ++count;
      }
    }
  }
  return count;
}
```

## Velocity

To compute the velocity of a particle (a node in the power diagram), we need to compute the gradient of the velocity potential over the corresponding cell in the regular triangulation. We can use CGAL here to do the hard work for us. The d-dimensional geometry kernel lets us compute the hyperplane associated with the four vertices of the cell in the triangulation.

``` {.cpp file=src/adhesion_velocity.cc}
#include "adhesion.hh"
#include <limits>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Kernel_d/Hyperplane_d.h>

typedef CGAL::Cartesian_d<double>    LiftedK;
typedef CGAL::Point_d<LiftedK>       LiftedPoint;
typedef CGAL::Hyperplane_d<LiftedK>  HyperPlane;

<<velocity-define-infinity>>

<<velocity-lifted-point>>

Vector Adhesion::velocity(RT::Cell_handle c) const {
  <<velocity-implementation>>
}
```

We'll need a point at infinity to give the hyperplane an orientation.

``` {.cpp #velocity-define-infinity}
constexpr double infinity
  = std::numeric_limits<double>::infinity();
```

CGAL's d-dimensional kernel needs to be told how many coordinates there are in a point, so we write a little wrapper.

``` {.cpp #velocity-lifted-point}
inline LiftedPoint lifted_point(
    double x, double y, double z, double w)
{
  double p[4] = { x, y, z, w };
  return LiftedPoint(4, p, p + 4);
}
```

We first need to convert the cell handle to its four lifted points.

``` {.cpp #velocity-implementation}
LiftedPoint points[4];

for (unsigned i = 0; i < 4; ++i)
{
  Weighted_point wp = rt.point(c, i);
  auto p = wp.point();
  auto w = wp.weight();
  points[i] = lifted_point(p.x(), p.y(), p.z(), w);
}
```

Then we create the hyperplane associated with these points, taking care to have the orientation such that the normal is pointing in possitive `w` direction. This is done by having the guide point $(0, 0, 0, -\infty)$ on the negative side of the hyperplane.

``` {.cpp #velocity-implementation}
auto guide = lifted_point(0, 0, 0, -infinity);
HyperPlane h(points, points + 4, guide, CGAL::ON_NEGATIVE_SIDE);
```

Given the normal $\vec{n}$, the velocity vector is given by
$$v_i = \frac{n_i}{2 t n_w},$$
where $i$ indexes the $x$, $y$ and $z$ components.

``` {.cpp #velocity-implementation}
auto normal = h.orthogonal_vector();
auto v  = normal / (2 * time * normal[3]);

return Vector(v[0], v[1], v[2]);
```

## Getting the power diagram


``` {.cpp file=src/adhesion_get_walls.cc}
#include "adhesion.hh"
#include "power_diagram.hh"

Mesh<Point, double> Adhesion::get_walls(
    double threshold) const
{
  return power_diagram_faces(rt, threshold);
}
```

``` {.cpp file=src/power_diagram.hh}
#include "cgal_base.hh"
#include "mesh.hh"

extern Mesh<Point, double> power_diagram_faces(
  RT const &rt, double threshold);
```

``` {.cpp file=src/power_diagram.cc}
#include "power_diagram.hh"

Mesh<Point, double> power_diagram_faces(
  RT const &rt,
  double threshold)
{
  Mesh<Point, double> mesh;

  <<power-diagram-dual-vertex>>

  for (auto e = rt.finite_edges_begin();
       e != rt.finite_edges_end();
       ++e)
  {
    std::vector<unsigned> vs;

    double l = rt.segment(*e).squared_length();
    if (l < threshold) continue;

    auto first = rt.incident_cells(*e), c = first;
    bool ok = true;
    do {
      if (rt.is_infinite(++c)) {
          ok = false;
          break;
      }

      vs.push_back(get_dual_vertex(c));
    } while (c != first);

    if (ok) {
      mesh.polygons.emplace_back(vs, sqrt(l));
    }
  }

  return mesh;
}
```

### Dual vertex

Every cell in the regular triangulation is associated with a vertex in the power diagram. We write a small helper function that obtains this dual vertex and caches it in a map.

``` {.cpp #power-diagram-dual-vertex}
std::map<RT::Cell_handle, unsigned> cell_index;

auto get_dual_vertex = [&rt, &cell_index, &mesh] (
    RT::Cell_handle const &h) -> unsigned
{
  if (cell_index.count(h) == 0)
  {
    cell_index[h] = mesh.vertices.size();
    mesh.vertices.push_back(rt.dual(h));
  }

  return cell_index[h];
};
```

# The main program

## Configuration

We read the configuration from a YAML file. Let's take the latest values from the Planck collaboration.

``` {.yaml #default-config file=examples/lcdm.yaml}
box:
  N:      128      # logical box size
  L:      32.0     # physical box size

cosmology:
  power-spectrum: Eisenstein & Hu
  h:        0.674   # Hubble parameter / 100
  ns:       0.965   # primordial power spectrum index
  Omega0:   1.0     # density in units of critical density
  sigma8:   0.811   # amplitude over 8 Mpc/h

run:
  seed:     8
  time:     [1.0, 2.0, 3.0]

output:
  hdf5:            output/lcdm.h5
  walls:           output/lcdm-{time:02.1f}-walls.obj
  wall-threshold:  10.0
```

## Run function

``` {.cpp file=src/run.hh}
#pragma once
#include <yaml-cpp/yaml.h>
#include <fmt/format.h>

extern void run(YAML::Node const &config);
```

``` {.cpp file=src/run.cc}
#include <iostream>
#include <exception>
#include <H5Cpp.h>

#include "run.hh"
#include "initial_conditions.hh"
#include "adhesion.hh"
#include "sphere.hh"
#include "write_selection_to_obj.hh"

void run(YAML::Node const &config)
{
  <<workflow>>
}
```

### Create box

``` {.cpp #workflow}
std::cerr << "Using box with parameters:\n"
          << config["box"] << "\n";
BoxParam box(
  config["box"]["N"].as<int>(),
  config["box"]["L"].as<double>());
```

### Generate initial conditions

``` {.cpp #workflow}
std::cerr << "Generating initial conditions:\n"
          << config["cosmology"] << "\n";
auto seed = config["run"]["seed"].as<unsigned long>();
auto field = generate_white_noise(box, seed);
compute_potential(box, *field, EisensteinHu(config["cosmology"]));
```

### Write initial conditions to file

``` {.cpp #workflow}
std::string output_filename = config["output"]["hdf5"].as<std::string>();
H5::H5File file(output_filename, H5F_ACC_TRUNC);
std::array<hsize_t, 3> shape = { box.N, box.N, box.N };
H5::DataSpace dataspace(3, shape.data());
H5::FloatType datatype(H5::PredType::NATIVE_DOUBLE);
H5::DataSet dataset = file.createDataSet("potential", datatype, dataspace);
dataset.write(field->data(), H5::PredType::NATIVE_DOUBLE);
```

### Compute adhesion model

``` {.cpp #workflow}
std::vector<double> time;
auto time_cfg = config["run"]["time"];
switch (time_cfg.Type()) {
  case YAML::NodeType::Scalar:
    time.push_back(time_cfg.as<double>());
    break;
  case YAML::NodeType::Sequence:
    time = time_cfg.as<std::vector<double>>();
    break;
  default: throw std::runtime_error("No time given.");
}

double threshold = config["output"]["wall-threshold"].as<double>(0.0);
std::string walls_filename = config["output"]["walls"].as<std::string>();
Sphere<K> sphere(Point(box.L/2, box.L/2, box.L/2), 0.4 * box.L);

for (double t : time) {
  std::cerr << "Computing regular triangulation ...\n";
  std::cerr << "time: " << t << " \n";
  Adhesion adhesion(box, *field, t);

  <<run-write-obj>>
}
```

### Write a sphere to OBJ

``` {.cpp #run-write-obj}
auto walls = adhesion.get_walls(threshold);
std::cerr << "Walls: " << walls.vertices.size() << " vertices and " << walls.polygons.size() << " polygons.\n";
std::string filename = fmt::format(walls_filename, fmt::arg("time", t));
std::cerr << "writing to " << filename << "\n";
write_selection_to_obj(filename, walls, sphere);
```

## Main function

The main program has not many arguments. It reads configuration from standard input in the YAML format. Any binary data will be stored in an auxiliary HDF5 file.

``` {.cpp #main-arguments}
argagg::parser argparser {{
  { "help",    {"-h", "--help"},
    "Show this help message.", 0 },
  { "version", {"--version"},
    "Show the software version.", 0 },
  { "config",  {"-c", "--config"},
	  "Supply configuration file.", 1 }
}};
```

We include the default configuration as a fallback.

``` {.cpp file=src/main.cc}
#include <iostream>
#include <argagg/argagg.hpp>
#include <yaml-cpp/yaml.h>

#include "run.hh"

<<version>>

<<main-arguments>>

const char *default_config = R"YAML(
<<default-config>>
)YAML";

int main(int argc, char **argv)
{
  argagg::parser_results args;
  try {
    args = argparser.parse(argc, argv);
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  if (args["help"]) {
    std::cerr << "Adhesion model example code -- (C) 2018 Johan Hidding\n";
    std::cerr << argparser;
    return EXIT_SUCCESS;
  }

  if (args["version"]) {
    std::cerr << "amec v" << VERSION << "\n";
    return EXIT_SUCCESS;
  }

  YAML::Node config;
  if (args["config"]) {
    auto config_file = args["config"].as<std::string>();
    std::cerr << "Reading `" << config_file << "` for input.\n";
    config = YAML::LoadFile(config_file);
  } else {
    std::cerr << "No configuration given, proceeding with defaults.\n";
    config = YAML::Load(default_config);
  }

  run(config);

  return EXIT_SUCCESS;
}
```

# Appendix

## Keeping a mesh in memory

``` {.cpp file=src/mesh.hh}
#pragma once

#include <tuple>
#include <memory>
#include <vector>

template <typename T, typename Info>
class Decorated: public T
{
    Info _info;

public:
    using T::T;

    Decorated(T const &t, Info const &info_):
        T(t), _info(info_) {}

    Decorated(Info const &info_):
        _info(info_) {}

    Info const &info() const
    {
        return _info;
    }
};

template <typename Point>
using Polygon = std::tuple<
                  std::vector<Point> *,
                  std::vector<unsigned>>;

template <typename Point>
using PolygonPair = std::tuple<
                      std::vector<Point> *,
                      std::vector<unsigned>,
                      std::vector<unsigned>>;

template <typename Point, typename Info>
struct Mesh
{
  using PolygonData = Decorated<std::vector<unsigned>, Info>;

  std::vector<Point> vertices;
  std::vector<PolygonData> polygons;
};
```

## Cutting polygons

``` {.cpp file=src/split_polygon.hh}
#pragma once

#include "mesh.hh"
#include <algorithm>

template <typename Point, typename Surface>
PolygonPair<Point> split_polygon(
    Polygon<Point> const &polygon,
    Surface const &surface,
    bool closed = true)
{
  std::vector<Point> *vertices = std::get<0>(polygon);
  std::vector<unsigned> orig = std::get<1>(polygon), r1, r2;

  auto is_below = [&vertices, &surface] (unsigned i) -> bool {
    return (surface.oriented_side((*vertices)[i]) == -1);
  };

  if (closed)
    orig.push_back(orig.front());

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

### Cutting a spherical region

``` {.cpp file=src/sphere.hh}
#pragma once
#include <optional>

template <typename K>
class Sphere
{
  using Point   = typename K::Point_3;
  using Vector  = typename K::Vector_3;

  Point origin;
  double radius_squared;

public:
  Sphere(Point const &p, double r):
      origin(p), radius_squared(r*r) {}

  int oriented_side(Point const &p) const
  {
    double d = (p - origin).squared_length();

    if (d < radius_squared)
      return -1;

    if (d > radius_squared)
      return +1;

    return 0;
  }

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
};
```

### Writing an OBJ file

``` {.cpp file=src/write_obj.hh}
#pragma once
#include <iostream>

#include "mesh.hh"

template <typename Point>
void write_obj(std::ostream &out, Mesh<Point, double> const &mesh)
{
  for (Point const &p : mesh.vertices) {
    out << "v " << p << " 1.0\n";
  }
  out << "\n";

  std::vector<double> val;
  double min = 1e6, max = 0.0;
  for (auto const &p : mesh.polygons) {
    double a = p.info();
    if (a < min) min = a;
    if (a > max) max = a;
  }

  unsigned i = 0;
  for (auto const &p : mesh.polygons) {
    double a = p.info();
    // out << "vt " << (a - min)/(max - min) << " 0\n";
    out << "vt " << a << " 0\n";
    ++i;
  }
  out << "\n";

  i = 1;
  for (auto const &p : mesh.polygons) {
    out << "f";
    for (unsigned j : p)
        out << " " << j+1 << "/" << i;
    out << "\n";
    ++i;
  }
}
```

### Clean mesh

``` {.cpp file=src/clean_mesh.hh}
#pragma once
#include "mesh.hh"

template <typename Point, typename Info>
Mesh<Point, Info> clean(
    Mesh<Point, Info> const &source)
{
  Mesh<Point, Info> result;
  std::map<unsigned, unsigned> vertex_map;

  for (auto const &v : source.polygons)
  {
    result.polygons.emplace_back(v.info());

    for (unsigned i : v) {
      if (vertex_map.count(i) == 0) {
        vertex_map[i] = result.vertices.size();
        result.vertices.push_back(source.vertices[i]);
      }

      result.polygons.back().push_back(vertex_map[i]);
    }
  }

  return result;
}
```

### Write selection

``` {.cpp file=src/write_selection_to_obj.hh}
#pragma once
#include <iostream>
#include <fstream>

#include "write_obj.hh"
#include "mesh.hh"
#include "clean_mesh.hh"
#include "split_polygon.hh"

template <typename Point, typename Surface, typename Info>
Mesh<Point, Info> select_mesh(
    Mesh<Point, Info> const &mesh,
    Surface const &surface)
{
  Mesh<Point, Info> result;
  result.vertices = mesh.vertices;

  for (auto const &v : mesh.polygons)
  {
      Info info  = v.info();
      std::vector<unsigned> polygon_data(v);

      auto below = std::get<1>(split_polygon(
        Polygon<Point>(&result.vertices, polygon_data),
        surface));

      if (below.size() > 0) {
        result.polygons.emplace_back(below, info);
      }
  }
  
  return clean(result);
}

template <typename Surface>
void write_selection_to_obj(
    std::string const &filename,
    Mesh<Point, double> const &mesh,
    Surface const &selector)
{
  std::cerr << "writing output to: " << filename << "\n";
  std::ofstream fo(filename);
  write_obj(fo, select_mesh(mesh, selector));
  fo.close();
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
    if (hint != 0)
        fftw_free(hint);

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
    real_space(box_.size),
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
      for (double &z : real_space) z /= box.size;
  }
};
```

## Testing

``` {.cpp file=tests/main.cc}
#include <gtest/gtest.h>

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
```

``` {.cpp file=tests/initial_conditions.cc}
#include <gtest/gtest.h>
#include "initial_conditions.hh"
#include "xtensor/xmath.hpp"

TEST(InitialConditions, BoxParam) {
  BoxParam box(128, 100.0);
  EXPECT_EQ(box.N, 128);
  EXPECT_EQ(box.size, 2097152);
  EXPECT_FLOAT_EQ(box.L, 100.0);
  EXPECT_FLOAT_EQ(box.res, 0.78125);
}

TEST(InitialConditions, GeneratingNoise) {
  BoxParam box(128, 100.0);
  auto x = generate_white_noise(box, 0);
  ASSERT_TRUE(x);
  EXPECT_EQ(x->size(), box.size);
  double total = xt::mean(*x)[0];
  EXPECT_NEAR(total, 0.0, 1e-2);
}

TEST(InitialConditions, Fourier) {
  BoxParam box(128, 100.0);
  std::array<size_t, 3> loc = {0, 0, 0};
  double k1 = box.k_abs(loc);
  EXPECT_FLOAT_EQ(k1, 0.0);
  increment_index<3>(box.shape(), loc);
  EXPECT_EQ(loc[2], 1);
  double k2 = box.k_abs(loc);
  EXPECT_FLOAT_EQ(k2, 2*M_PI/100.0);
}
```
