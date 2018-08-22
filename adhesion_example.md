# Computing the adhesion model using C++ and CGAL

We present a (relatively) small example of using the CGAL library to run the adhesion model. A basic CGAL program can seem intimidating to start with. The CGAL manual provides basic examples that offer a good starting point. In our case we are interested in creating regular triangulations. We adapted the example from the manual.

## Introduction

This document is aimed to be self-containing. This means that all the code to build a working adhesion model is included. A significant fraction of the code in here is dealing with generating initial conditions for the actual model.

We've tried to limit the involvement of too much boilerplate code by using existing libraries where possible.

If you are only interested only in the CGAL parts, it should be safe to skip the section on generating cosmological initial conditions.

### Literate programming

This example is written in a style of *literate programming*. This document contains a complete and functioning example of working with CGAL to compute the adhesion model. For didactic reasons we don't always give the listing of an entire source file. In stead, we use a system of references known as *noweb*.

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

### Prerequisites

- C++ 14 compiler
- CGAL 4.12 - The Computational Geometry Algorithm Library
- XTensor 0.17 - A template library for handling array data in C++
- FFTW3 3.3 - The Fastest Fourier Transform in the West
- XTensor-FFTW - Tying FFTW3 to an XTensor interface
- hdf5-cpp - [HDF5](https://support.hdfgroup.org/HDF5/doc/cpplus_RM/index.html) is used to store large blobs of binary data and meta data.
- yaml-cpp - [YAML-cpp](https://github.com/jbeder/yaml-cpp) is a YAML parser for C++. We use it to parse configuration files.

## CGAL Geometry kernels

CGAL comes with a set of *geometry kernels*. Each kernel bundles basic type definitions like `Point`, `Vector`, `Circle`, etc. and geometric operations on those types. Depending on the requirements of the programmer, we can choose different implementations of these concepts. These implementations vary in representation of real numbers, vector quantities, and how geometric operations are computed on them. Some abstract applications require an exact representation of numbers while other more cosmological applications can afford to be more liberal with regards to exactness.

The algorithms that actually do the advanced geometric computations, like the regular triangulations we use, are implemented in generic terms using C++ template techniques. We need to supply those algorithms the correct geometry kernel for our application. This is why all CGAL programs start with a list of template type definitions.

In our case, what we need is a double precision floating point representation of numbers, while retaining logical consistency in geometric predicates, of which the co-linearity test is the most obvious example. This is provided by the `Exact_predicates_inexact_constructions_kernel` kernel.

We collect those type definitions in a separate header file:

``` {.cpp file=src/cgal_base.hh}
#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Regular_triangulation_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Regular_triangulation_euclidean_traits_3<K>   Traits;
typedef Traits::RT                                          Weight;
typedef Traits::Bare_point                                  Point;
typedef Traits::Vector_3                                    Vector;
typedef Traits::Weighted_point                              Weighted_point;
typedef CGAL::Regular_triangulation_3<Traits>               Rt;
```

Since we'll be using bare (unweighted) points, weighted points, and vectors we defined aliases for those types. Note that CGAL is particular about the difference between points and vectors. Points are locations without absolute properties, whereas vectors describe how to get from one point to the other. Internally they may have the same numerical representation, but this may not strictly be the case for all geometry kernels.

## Initial conditions

``` {.cpp file=src/initial_conditions.hh}
#pragma once
#include "boxparam.hh"

#include <memory>
#include <xtensor/xtensor.hpp>
#include <yaml-cpp/yaml.h>

extern std::unique_ptr<xt::xtensor<double, 3>> white_noise(
    BoxParam const &box,
    unsigned long seed);

<<power-spectra>>

extern void apply_power_spectrum(
    BoxParam const &box,
    xt::xtensor<double, 3> &field,
    PowerSpectrum const &P);
```

### The simulation box

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

  std::array<size_t, 3> rfft_shape() const
  {
    return std::array<size_t, 3>{ N, N, N/2 + 1 };
  }

  size_t rfft_size() const
  {
    return N * N * (N / 2 + 1);
  }

  double wave_number(int i) const
  {
    return (int(i) > int(N)/2 ? int(i) - int(N) : int(i)) * (2*M_PI)/L;
  }

  double k_abs(std::array<size_t, 3> const &loc) const
  {
    double x = 0.0;
    for (size_t i : loc)
      x += sqr(wave_number(i));
    return sqrt(x);
  }

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
  int z = i / (N*N);

  return Point(x * res, y * res, z * res);
}
```

Note that we set the $x$-coordinate to be the fastest changing coordinate in the flattened array. This is known as *row-major* ordering, which the same as how indexing into C/C++ and Python/NumPy arrays works. Later on, we will be adding more methods to the `BoxParam` structure.

#### Iterating multi-dimensional arrays

We'll be indexing multi-dimensional arrays. To prevent having to write nested for-loops, we use the `increment_index` helper function.

``` {.cpp #increment-index}
template <unsigned R>
inline unsigned increment_index(
    std::array<size_t, R> const &shape,
    std::array<size_t, R> &index)
{
  for (unsigned i = 0; i < R; ++i) {
    if (++index[i] < shape[i])
      return i;

    if (i == R - 1)
      return R;

    index[i] = 0;
  }
}
```

### White noise

``` {.cpp file=src/white_noise.cc}
#include "initial_conditions.hh"
#include <random>

std::unique_ptr<xt::xtensor<double, 3>> white_noise(
    BoxParam const &box, unsigned long seed)
{
  auto result = std::make_unique<xt::xtensor<double, 3>>(box.shape());

  std::mt19937 random(seed);
  std::normal_distribution<double> normal;

  for (double &value : *result)
  {
    value = normal(random);
  }

  return std::move(result);
}
```

### Power spectrum

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

#### Eisenstein-Hu power spectrum

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
  });
}
```

#### Normalisation

### Applying the power spectrum

This takes some Fourier wizardry.

``` {.cpp file=src/apply_power_spectrum.cc}
#include "initial_conditions.hh"
#include <xtensor-fftw/basic.hpp>

void apply_power_spectrum(
    BoxParam const &box,
    xt::xtensor<double, 3> &field,
    PowerSpectrum const &P)
{
  auto f_k = xt::fftw::rfft(field);
  auto f_shape = box.rfft_shape();

  std::array<size_t, 3> loc(0);
  for (size_t i = 0; i < box.rfft_size(); ++i) {
    f_k[i] *= sqrt(P(box.k_abs(loc)));
    increment_index(f_shape, loc);
  }

  field = xt::fftw::irfft(f_k) / box.size();
}
```

## Lloyd Iteration

## Detecting Structures

### Computing the triangulation

### Filtering for structures

## Input/Output

