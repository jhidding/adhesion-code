[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1477536.svg)](https://doi.org/10.5281/zenodo.1477536)

# CGAL Adhesion model example code

We present a (relatively) small example of using the CGAL library ([http://cgal.org](http://cgal.org)) to run the adhesion model. This literate C++ code generates an initial potential field and computes the *regular triangulation* to that potential, which is a weighted generalisation of the *Delaunay triangulation*. The output is a selection of its dual, the power diagram or weighted Voronoi tessellation, written in a form that is ready for analysis and visualisation.

![The output structure](figures/web-evolution.png)

## Documentation

Go to [https://jhidding.github.io/adhesion-code](https://jhidding.github.io/adhesion-code).

## Literate programming

This example is written in a style of *literate programming* [@Knuth1984]. This document contains a complete and functioning example of working with CGAL to compute the adhesion model. For didactic reasons we don't always give the listing of an entire source file in one go. In stead, we use a system of references known as *noweb* [@Ramsey1994].

Inside source fragments you may encounter a line with `<<...>>` marks like,

*file: «examples/hello_world.cc»=*
```cpp
#include <cstdlib>
#include <iostream>

<<example-main-function>>
```

which is then elsewhere specified. Order doesn't matter,

*«hello-world»=*
```cpp
std::cout << "Hello, World!" << std::endl;
```

So we can reference the `<<hello-world>>` code block later on.

*«example-main-function»=*
```cpp
int main(int argc, char **argv) {
  <<hello-world>>
}
```

A definition can be appended with more code as follows (in this case, order does matter!):

*«hello-world»+*
```cpp
return EXIT_SUCCESS;
```

These blocks of code can be *tangled* into source files. The source code presented in this report combine into a fully working example of the adhesion model!

## Build instructions

This code has been tested on Debian testing, Ubuntu 18.04, CentOS (with Conda), and MacOS X. There is a `Dockerfile` in the repository that builds the application as well as the `pdf` of the report successfully.

### Prerequisites

For building the software:

| Package  | version | description |
|----------|---------|-------------|
| C++ compiler | C++17 standard | Tested with GCC 8. |
| GNU Make | - | - |
| CGAL     | ≥4.12   | [The Computational Geometry Algorithm Library](http://cgal.org) |
| FFTW3    | ≥3.3    | [The Fastest Fourier Transform in the West](http://www.fftw.org/) |
| hdf5-cpp | ≥1.8.13 | [HDF5](https://support.hdfgroup.org/HDF5/doc/cpplus_RM/index.html) is used to store large blobs of binary data and meta data. |
| yaml-cpp | ≥0.5    | [YAML-cpp](https://github.com/jbeder/yaml-cpp) is a YAML parser for C++. We use it to parse configuration files. |
| argagg   | ≥0.4.6  | [ArgAgg](https://github.com/vietjtnguyen/argagg) stands for Argument Aggregator and is a C++ command-line argument parser. |
| fmt      | ≥4.1    | [fmt](http://fmtlib.net/latest/index.html) is a string formatting library that has a similar interface as Python's. |
| Pandoc   | ≥2.2.3  | [Pandoc](http://pandoc.org/) is a universal document converter. To build this example from the markdown, you need version 2.2.3 or higher and the `pandoc-citeproc` extension. |

These requirements are available in Debian testing or Ubuntu 18.04, with the exception of Pandoc for which the versions may still be too old. More recent binary releases are available on [http://pandoc.org/](http://pandoc.org) for all major distributions of Linux, MacOS and Windows.

### LaTeX

To create the PDF version of the report, `xelatex` and a good many fonts are needed.

### Running on Debian/Ubuntu

Get the following packages

```shell
sudo apt install \
    argagg-dev cmake g++ libcgal-dev libfftw3-dev \
    libfmt-dev libgsl-dev libhdf5-dev libyaml-cpp-dev \
    lmodern make pkg-config python3-pip rsync texlive \
    texlive-fonts-extra texlive-latex-extra \
    texlive-latex-recommended texlive-xetex wget
```

Install a recent version of `pandoc`,

```shell
wget https://github.com/jgm/pandoc/releases/download/2.4/pandoc-2.4-1-amd64.deb
dpkg -i ./pandoc-2.4-1-amd64.deb
```

Install the `pandoc-fignos` and `pandoc-eqnos` plugins,

```shell
pip install --user pandoc-fignos pandoc-eqnos
```

Go to the root folder of this package and inspect the `Makefile`. The `~/.local/include` directory is already configured as include path. Build the executable by running `make`. Build the PDF by running `make report`.

### Running on another GNU/Linux using Conda

- install Pandoc, CGAL and HDF5 using Conda

    conda install -c conda-forge pandoc cgal hdf5

- download YAML-cpp and install with `cmake`

```
git clone https://github.com/jbeder/yaml-cpp.git
cd yaml-cpp && mkdir build && cd build
ccmake ..
# change install path to your liking
make install
```

- download `fmtlib` and install with `cmake`

    git clone https://github.com/fmtlib/fmt.git

TIP: if you would like to keep your `~/.local` directory clean, or use some form of package management, take a look at [`xstow`](http://xstow.sourceforge.net/).

### Running on Mac

Several dependencies can be installed using Homebrew.

    brew install pandoc cgal yaml-cpp hdf5 fmt librsvg gsl pkg-config

Some others need to be installed through `pip`:

    pip install pandoc-eqnos pandoc-fignos

ArgAgg is available as a single header file which is included in this repository.

You may want to change the `compiler` and `linker` variables in the `Makefile` into `clang++`. Also, comment out the `cgal_cflags` variable as `clang` does not support `-frounding-math`. You may need to remove `-lboost_thread` from `cgal_libs` as well.

To build the code run `make`, to build the report run `make report`.

### Getting coverage

Compile with `--coverage`, then run code, then:

```shell
lcov --capture --directory . --output coverage.info --no-external
genhtml coverage.info --output-directory out
```

## License

This package is distributed under the Apache v2 license. See the `LICENSE` file for more information.
