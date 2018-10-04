## Convex conjugate method

``` {.cpp file=src/convex_conjugate.hh}
#pragma once
#include "boxparam.hh"
#include <vector>
#include <array>

struct CCResult {
  std::vector<double> potential;
  std::vector<std::array<double>> map;
};

extern CCResult convex_conjugate(
  BoxParam const &box,
  std::vector<double> const &potential,
  double t);
```

``` {.cpp file=src/convex_conjugate/integer_cc.cc}
std::vector<size_t> integer_convex_conjugate(
  BoxParam const &box,
  std::vector<double> const &potential,
  double t)
{
  <<icc-distance-function>>
  <<icc-initialise-result>>
  <<icc-main-loop>>
}
```

The distance function is defined as

$$d({\bf x}, {\bf q}) = ({\bf x} - {\bf q})^2 + 2 t \Phi_0({\rm q}).$${#eq:distance-function}

``` {.cpp #icc-distance-function}
auto distance_function = [&] (size_t x, size_t q) {
  return sqr(box.point(x) - box.point(q)) + 2 * t * potential[q];
};
```

In the initial state, each point maps onto itself,

$${\bf q}(t = 0) = {\bf x}$$.

We also define a buffer, to which each iteration will write. At the end of each iteration the buffer and the result are swapped.

``` {.cpp #icc-initialise-result}
std::vector<size_t> result(box.size()),
                    buffer(box.size());
for (size_t x = 0; x < box.size(); ++x)
  result[x] = x;
```

We loop until the result converges

``` {.cpp #icc-main-loop}
bool converged;
do {
  converged = true;
  <<icc-iterate>>
} while (not converged);
```

Each iteration we compare distance functions with neighbours in a *stencil operation*.

``` {.cpp #icc-iterate}
for (size_t x = 0; x < box.size(); ++x) {
  size_t q = result[x];
  double d = distance_function(x, q);

  std::array<size_t, 3> index = box.iloc(x);
  for (unsigned k = 0; k < 3; ++k) {
    std::array<size_t, 3> stencil = index;
    for (int delta = -1; delta <= 1; delta += 2) {
      stencil[k] = modulo(index[k] + delta, box.N);

      size_t xn = box.iloc(stencil);
      size_t qn = result[xn];
      if (qn != q) {
        double dn = distance_function(x, qn);
        if (dn < d) {
          d = dn;
          q = qn;
          converged = false;
        }
      }
    }
  }
}
```