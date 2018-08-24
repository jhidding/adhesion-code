# Computing the adhesion model using C++ and CGAL

We present a (relatively) small example of using the CGAL library to run the adhesion model. A basic CGAL program can seem intimidating to start with. The CGAL manual provides basic examples that offer a good starting point. In our case we are interested in creating regular triangulations. We adapted the example from the manual.

### Version

``` {.cpp #version}
#define VERSION "0.1"
```

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
- argagg - [ArgAgg](https://github.com/vietjtnguyen/argagg) stands for Argument Aggregator and is a C++ command-line argument parser.

All of these packages are available in the Debian GNU/Linux package repositories.

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

Since we'll be using bare (weightless) points, weighted points, and vectors we defined aliases for those types. Note that CGAL is particular about the difference between points and vectors. Points are locations without absolute properties, whereas vectors describe how to get from one point to the other. Internally they may have the same numerical representation, but this may not strictly be the case for all geometry kernels.

## Initial conditions

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

#### Fourier properties

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

#### Iterating multi-dimensional arrays

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

### White noise

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
  };
}
```

#### Normalisation

### Applying the power spectrum

This takes some Fourier wizardry.

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

## Lloyd Iteration

## Detecting Structures

### Computing the triangulation

### Filtering for structures

## The main program

### Configuration

We read the configuration from a YAML file. Let's take the latest values from the Planck collaboration.

``` {.yaml #default-config file=examples/lcdm150.yaml}
box:
  N:      128       # logical box size
  L:      150.0     # physical box size

cosmology:
  power-spectrum: Eisenstein & Hu
  h:        0.674   # Hubble parameter / 100
  ns:       0.965   # primordial power spectrum index
  Omega0:   1.0     # density in units of critical density
  sigma8:   0.811   # amplitude over 8 h⁻¹ Mpc

run:
  seed:     0
  time:     1.0
  output:   lcdm150.h5
```

### Run function

``` {.cpp file=src/run.hh}
#pragma once
#include <yaml-cpp/yaml.h>

extern void run(YAML::Node const &config);
```

``` {.cpp file=src/run.cc}
#include <iostream>
#include <H5Cpp.h>

#include "run.hh"
#include "initial_conditions.hh"

void run(YAML::Node const &config)
{
  <<workflow>>
}
```

#### Create box

``` {.cpp #workflow}
std::cerr << "Using box with parameters:\n"
          << config["box"] << "\n";
BoxParam box(
  config["box"]["N"].as<int>(),
  config["box"]["L"].as<double>());
```

#### Generate initial conditions

``` {.cpp #workflow}
std::cerr << "Generating initial conditions:\n"
          << config["cosmology"] << "\n";
auto seed = config["run"]["seed"].as<unsigned long>();
auto field = generate_white_noise(box, seed);
compute_potential(box, *field, EisensteinHu(config["cosmology"]));
```

#### Write initial conditions to file

``` {.cpp #workflow}
std::string output_filename = config["run"]["output"].as<std::string>();
H5::H5File file(output_filename, H5F_ACC_TRUNC);
std::array<hsize_t, 3> shape = { box.N, box.N, box.N };
H5::DataSpace dataspace(3, shape.data());
H5::FloatType datatype(H5::PredType::NATIVE_DOUBLE);
H5::DataSet dataset = file.createDataSet("potential", datatype, dataspace);
dataset.write(field->data(), H5::PredType::NATIVE_DOUBLE);
```

### Main function

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
