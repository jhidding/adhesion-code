// ~\~ language=C++ filename=src/boxparam.hh
// ~\~ begin <<adhesion_example.md|src/boxparam.hh>>[0]
#pragma once
#include <cstdlib>
#include <cmath>
#include <array>

// ~\~ begin <<adhesion_example.md|helper-functions>>[0]
template <typename T>
inline T sqr(T x) { return x*x; }
// ~\~ end
// ~\~ begin <<adhesion_example.md|increment-index>>[0]
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
// ~\~ end

class BoxParam {
public:
  static constexpr unsigned R = 3; // dimension of the box
  unsigned N;         // number of grid points in a row
  double   L;         // physical size of the box

public:
  BoxParam(unsigned N_, double L_)
    : N(N_)
    , L(L_)
  {}

  // ~\~ begin <<adhesion_example.md|fourier-properties>>[0]
  std::array<size_t, 3> rfft_shape() const {
    return std::array<size_t, 3>{ N, N, N/2 + 1 };
  }

  size_t rfft_size() const {
    return N * N * (N / 2 + 1);
  }
  // ~\~ end
  // ~\~ begin <<adhesion_example.md|fourier-properties>>[1]
  double wave_number(int i) const {
    return ( int(i) > int(N)/2
           ? int(i) - int(N)
           : int(i) ) * (2*M_PI)/L;
  }
  // ~\~ end
  // ~\~ begin <<adhesion_example.md|fourier-properties>>[2]
  double k_abs(std::array<size_t, 3> const &loc) const {
    double x = 0.0;
    for (size_t i : loc)
      x += sqr(wave_number(i));
    return sqrt(x);
  }
  // ~\~ end
  // ~\~ begin <<adhesion_example.md|boxparam-methods>>[0]
  std::array<size_t, 3> shape() const {
    return std::array<size_t, 3>({ N, N, N });
  }

  size_t size() const {
    return N * N * N;
  }
  // ~\~ end
  // ~\~ begin <<adhesion_example.md|boxparam-methods>>[1]
  std::array<size_t, 3> iloc(size_t i) const {
    size_t x = i % N,
           y = (i / N) % N,
           z = i / (N * N);

    return {z, y, x};
  }
  // ~\~ end
  // ~\~ begin <<adhesion_example.md|boxparam-methods>>[2]
  template <typename Point>
  Point point(size_t i) const {
    auto p = iloc(i);
    return Point(p[2] * L/N, p[1] * L/N, p[0] * L/N);
  }
  // ~\~ end
};
// ~\~ end
