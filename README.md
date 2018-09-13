# CGAL Adhesion model example code

We present a (relatively) small example of using the CGAL library to run the adhesion model. The adhesion model is a model of cosmic structure formation. It computes the shape of emerging structures in the cosmic web straight from an initial potential.

![The output structure](figures/orbit.gif)

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

Some components of this package are bleeding edge. There is a `Dockerfile` in the repository that builds the application as well as the `pdf` of the report successfully. To build this application on your native system can be a bit challenging.

### Prerequisites

| Package  | version | description |
|----------|---------|-------------|
| C++ compiler | C++17 standard | Tested with GCC 8. |
| GNU Make | - | - |
| CGAL     | ≥4.12   | [The Computational Geometry Algorithm Library](http://cgal.org) |
| XTensor  | ≥0.17   | [XTensor](http://quantstack.net/xtensor) is a template library for handling array data in C++. |
| FFTW3    | ≥3.3    | [The Fastest Fourier Transform in the West](http://www.fftw.org/) |
| hdf5-cpp | ≥1.8.13 | [HDF5](https://support.hdfgroup.org/HDF5/doc/cpplus_RM/index.html) is used to store large blobs of binary data and meta data. |
| yaml-cpp | ≥0.5    | [YAML-cpp](https://github.com/jbeder/yaml-cpp) is a YAML parser for C++. We use it to parse configuration files. |
| argagg   | ≥0.4.6  | [ArgAgg](https://github.com/vietjtnguyen/argagg) stands for Argument Aggregator and is a C++ command-line argument parser. |
| fmt      | ≥4.1    | [fmt](http://fmtlib.net/latest/index.html) is a string formatting library that has a similar interface as Python's. |
| Pandoc   | ≥2.2.3  | [Pandoc](http://pandoc.org/) is a universal document converter. To build this example from the markdown, you need version 2.2.3 or higher and the `pandoc-citeproc` extension. |

These requirements are available in Debian testing, except for XTensor, which you need to install manually from [the XTensor github page](https://github.com/quantstack/xtensor).

### LaTeX

To create the PDF version of the report, `xelatex` and a good many fonts are needed.

### Running on Debian (testing)

(as root) Get the following packages

```shell
sudo apt install \
    argagg-dev cmake g++ libcgal-dev libfftw3-dev \
    libfmt-dev libhdf5-dev libyaml-cpp-dev lmodern make \
    pkg-config rsync texlive texlive-fonts-extra \
    texlive-latex-extra texlive-latex-recommended \
    texlive-xetex wget git
```

Install a very recent version of `pandoc`

```shell
wget https://github.com/jgm/pandoc/releases/download/2.2.3.2/pandoc-2.2.3.2-1-amd64.deb
sudo dpkg -i ./pandoc-2.2.3.2-1-amd64.deb
```

Get the development version of XTensor, install as user:

```shell
git clone https://github.com/QuantStack/xtl.git
git clone https://github.com/QuantStack/xtensor.git

cd xtl
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=~/.local
make install

cd ../../xtensor
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=~/.local
make install
```

Go to the root folder of this package and inspect the `Makefile`. The `~/.local/include` directory is already configured as include path. Build the executable by running `make`. Build the PDF by running `make report`.

### Running on Mac

A beta-tester has been found. More's to follow.

### Getting coverage

compile with `--coverage`, then run code, then:

```shell
lcov --capture --directory . --output coverage.info --no-external
genhtml coverage.info --output-directory out
```

## License

This package is distributed under the Apache v2 license. See the `LICENSE` file for more information.