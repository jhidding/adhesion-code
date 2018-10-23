## Cubic splines

``` {.cpp src/spline.hh}
namespace spline {
  template <unsigned R>
  class Spline {};

  template <>
  class Spline<K, 3> {
  public:
    constexpr unsigned R = 3;
    using Point = std::array<double, 3>;
    using Vector = std::array<double, 3>;
    using Index = std::array<size_t, 3>;

  private:
    static double const A[4096];
    BoxParam const &box;
    std::vector<double> const &data;

  public:
    Spline(BoxParam box_, std::vector<double> const &data_)
      : box(box_)
      , data(data_)
    {}

    Array<double> calc_coeff(Index const &S) const;
    double operator()(Point const &x) const;
    double f(Point const &x) const;
    Vector df(Point const &x) const;
    std::pair<double, Vector> fdf(Point const &x) const;
  };
}
```

``` {.cpp file=src/spline/tricubic.cc}
Array<double> Spline<3>::calc_coeff(Spline<3>::Index const &S) const
{
  std::vector<double> b(64);
  auto const &dx = box->dx;

  auto g = [&] (iVector<R> const &X)
  {
    return data[box->idx(X)];
  };

  for (unsigned i = 0; i < 8; ++i)
  {
    iVector<R>  X = S + box->block[i];

    b[i]    = g(X);

    b[i + 8]  = (g(X + dx[0]) - g(X - dx[0])) / 2;

    b[i + 16]  = (g(X + dx[1]) - g(X - dx[1])) / 2;

    b[i + 24]  = (g(X + dx[2]) - g(X - dx[2])) / 2;

    b[i + 32]  = ((g(X + dx[1] + dx[0]) - g(X - dx[1] + dx[0])) -
                       (g(X + dx[1] - dx[0]) - g(X - dx[1] - dx[0]))) / 4;

    b[i + 40]  = ((g(X + dx[2] + dx[1]) - g(X - dx[2] + dx[1])) -
                       (g(X + dx[2] - dx[1]) - g(X - dx[2] - dx[1]))) / 4;

    b[i + 48]  = ((g(X + dx[0] + dx[2]) - g(X - dx[0] + dx[2])) -
                       (g(X + dx[0] - dx[2]) - g(X - dx[0] - dx[2]))) / 4;
    
    b[i + 56]  = (g(X + dx[0] + dx[1] + dx[2]) - g(X + dx[0] + dx[1] - dx[2])
        -  g(X + dx[0] - dx[1] + dx[2]) + g(X + dx[0] - dx[1] - dx[2])
        -  g(X - dx[0] + dx[1] + dx[2]) + g(X - dx[0] + dx[1] - dx[2])
        +  g(X - dx[0] - dx[1] + dx[2]) - g(X - dx[0] - dx[1] - dx[2])) / 8;
  }

  Array<double> alpha(64, 0.0);
  MdRange<2> Z(iVector<2>(64));

  for (unsigned i = 0; i < 4096; ++i)
  {
    iVector<2> z = Z[i];
    alpha[z[1]] += b[z[0]] * A[i];
  }

  return alpha;
}

double Spline<3>::operator()(dVector<3> const &x) const
{
  return f(x);
}

double Spline<3>::f(dVector<3> const &x_) const
{
  dVector<R> x = x_ / box->res();
  iVector<R> origin = floor_cast(x);
  dVector<R> P = x - dVector<R>(origin);

  auto c = calc_coeff(origin);
  MdRange<R> loop(iVector<R>(4));

  // P(x, y, z) = Sum_i Sum_j Sum_k [c_ijk * x^i y^j z^k]
  double s = 0;
  for (unsigned i = 0; i < 64; ++i)
  {
    iVector<R> e = loop[i];
    s += c[i] * pow(P[0], e[0]) * pow(P[1], e[1]) * pow(P[2], e[2]);
  }

  return s;
}

dVector<3> Spline<3>::df(dVector<3> const &x_) const
{
  dVector<R> x = x_ / box->res();
  iVector<R> origin = floor_cast(x);
  dVector<R> P = x - dVector<R>(origin);

  auto c = calc_coeff(origin);
  MdRange<R> loop(iVector<R>(4));

  dVector<R> v(0);

  for (unsigned n = 0; n < 64; ++n)
  {
    iVector<R> e = loop[n];

    for (unsigned k = 0; k < R; ++k)
    {
      unsigned i = (k + 1) % R;
      unsigned j = (k + 2) % R;

      if (e[k] != 0)
      {
        v[k] += c[n] * e[k] * pow(P[k], e[k] - 1) 
              * pow(P[i], e[i]) * pow(P[j], e[j]);
      }
    }
  }

  return v;
}

std::pair<double, dVector<3>> Spline<3>::fdf(dVector<3> const &x_) const
{
  dVector<R> x = x_ / box->res();
  iVector<R> origin = floor_cast(x);
  dVector<R> P = x - dVector<R>(origin);

  auto c = calc_coeff(origin);
  MdRange<R> loop(iVector<R>(4));

  double s = 0;
  dVector<R> v(0);
  for (unsigned n = 0; n < 64; ++n)
  {
    iVector<R> e = loop[n];
    s += c[n] * pow(P[0], e[0]) * pow(P[1], e[1]) * pow(P[2], e[2]);

    for (unsigned k = 0; k < R; ++k)
    {
      unsigned i = (k + 1) % R;
      unsigned j = (k + 2) % R;

      if (e[k] != 0)
      {
        v[k] += c[n] * e[k] * pow(P[k], e[k] - 1) 
              * pow(P[i], e[i]) * pow(P[j], e[j]);
      }
    }
  }

  return std::make_pair(s, v);
}

```

``` {.cpp src/spline/tricubic_matrix.cc}
#include "spline.hh"

double const Spline<3>::A[4096] = {
  1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
 -3,   3,   0,   0,   0,   0,   0,   0,  -2,  -1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  2,  -2,   0,   0,   0,   0,   0,   0,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -3,   3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
 -2,  -1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   2,  -2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
 -3,   0,   3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -2,   0,  -1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,  -3,   0,   3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
 -2,   0,  -1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  9,  -9,  -9,   9,   0,   0,   0,   0,   6,   3,  -6,  -3,   0,   0,   0,   0,   6,  -6,   3,  -3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  4,   2,   2,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
 -6,   6,   6,  -6,   0,   0,   0,   0,  -3,  -3,   3,   3,   0,   0,   0,   0,  -4,   4,  -2,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
 -2,  -2,  -1,  -1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  2,   0,  -2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   2,   0,  -2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  1,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
 -6,   6,   6,  -6,   0,   0,   0,   0,  -4,  -2,   4,   2,   0,   0,   0,   0,  -3,   3,  -3,   3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
 -2,  -1,  -2,  -1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  4,  -4,  -4,   4,   0,   0,   0,   0,   2,   2,  -2,  -2,   0,   0,   0,   0,   2,  -2,   2,  -2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  1,   1,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -3,   3,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -2,  -1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   2,  -2,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,  -3,   3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -2,  -1,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   2,  -2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   1,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -3,   0,   3,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,  -2,   0,  -1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -3,   0,   3,   0,   0,   0,   0,   0,  -2,   0,  -1,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   9,  -9,  -9,   9,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   6,  -6,   3,  -3,   0,   0,   0,   0,   6,   3,  -6,  -3,   0,   0,   0,   0,   4,   2,   2,   1,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -6,   6,   6,  -6,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,  -4,   4,  -2,   2,   0,   0,   0,   0,  -3,  -3,   3,   3,   0,   0,   0,   0,  -2,  -2,  -1,  -1,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   2,   0,  -2,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   2,   0,  -2,   0,   0,   0,   0,   0,   1,   0,   1,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -6,   6,   6,  -6,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,  -3,   3,  -3,   3,   0,   0,   0,   0,  -4,  -2,   4,   2,   0,   0,   0,   0,  -2,  -1,  -2,  -1,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   4,  -4,  -4,   4,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   2,  -2,   2,  -2,   0,   0,   0,   0,   2,   2,  -2,  -2,   0,   0,   0,   0,   1,   1,   1,   1,   0,   0,   0,   0,
 -3,   0,   0,   0,   3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -2,   0,   0,   0,  -1,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,  -3,   0,   0,   0,   3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -2,   0,   0,   0,  -1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  9,  -9,   0,   0,  -9,   9,   0,   0,   6,   3,   0,   0,  -6,  -3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   6,  -6,   0,   0,   3,  -3,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   4,   2,   0,   0,   2,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
 -6,   6,   0,   0,   6,  -6,   0,   0,  -3,  -3,   0,   0,   3,   3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -4,   4,   0,   0,  -2,   2,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -2,  -2,   0,   0,  -1,  -1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -3,   0,   0,   0,   3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,  -2,   0,   0,   0,  -1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
 -3,   0,   0,   0,   3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -2,   0,   0,   0,  -1,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   9,  -9,   0,   0,  -9,   9,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  6,   3,   0,   0,  -6,  -3,   0,   0,   6,  -6,   0,   0,   3,  -3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   4,   2,   0,   0,   2,   1,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -6,   6,   0,   0,   6,  -6,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
 -3,  -3,   0,   0,   3,   3,   0,   0,  -4,   4,   0,   0,  -2,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -2,  -2,   0,   0,  -1,  -1,   0,   0,
  9,   0,  -9,   0,  -9,   0,   9,   0,   0,   0,   0,   0,   0,   0,   0,   0,   6,   0,   3,   0,  -6,   0,  -3,   0,   6,   0,  -6,   0,   3,   0,  -3,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   4,   0,   2,   0,   2,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   9,   0,  -9,   0,  -9,   0,   9,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  6,   0,   3,   0,  -6,   0,  -3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   6,   0,  -6,   0,   3,   0,  -3,   0,   4,   0,   2,   0,   2,   0,   1,   0,
-27,  27,  27, -27,  27, -27, -27,  27, -18,  -9,  18,   9,  18,   9, -18,  -9, -18,  18,  -9,   9,  18, -18,   9,  -9, -18,  18,  18, -18,  -9,   9,   9,  -9,
-12,  -6,  -6,  -3,  12,   6,   6,   3, -12,  12,  -6,   6,  -6,   6,  -3,   3, -12,  -6,  12,   6,  -6,  -3,   6,   3,  -8,  -4,  -4,  -2,  -4,  -2,  -2,  -1,
 18, -18, -18,  18, -18,  18,  18, -18,   9,   9,  -9,  -9,  -9,  -9,   9,   9,  12, -12,   6,  -6, -12,  12,  -6,   6,  12, -12, -12,  12,   6,  -6,  -6,   6,
  6,   6,   3,   3,  -6,  -6,  -3,  -3,   8,  -8,   4,  -4,   4,  -4,   2,  -2,   6,   6,  -6,  -6,   3,   3,  -3,  -3,   4,   4,   2,   2,   2,   2,   1,   1,
 -6,   0,   6,   0,   6,   0,  -6,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -3,   0,  -3,   0,   3,   0,   3,   0,  -4,   0,   4,   0,  -2,   0,   2,   0,
  0,   0,   0,   0,   0,   0,   0,   0,  -2,   0,  -2,   0,  -1,   0,  -1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,  -6,   0,   6,   0,   6,   0,  -6,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
 -3,   0,  -3,   0,   3,   0,   3,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -4,   0,   4,   0,  -2,   0,   2,   0,  -2,   0,  -2,   0,  -1,   0,  -1,   0,
 18, -18, -18,  18, -18,  18,  18, -18,  12,   6, -12,  -6, -12,  -6,  12,   6,   9,  -9,   9,  -9,  -9,   9,  -9,   9,  12, -12, -12,  12,   6,  -6,  -6,   6,
  6,   3,   6,   3,  -6,  -3,  -6,  -3,   6,  -6,   6,  -6,   3,  -3,   3,  -3,   8,   4,  -8,  -4,   4,   2,  -4,  -2,   4,   2,   4,   2,   2,   1,   2,   1,
-12,  12,  12, -12,  12, -12, -12,  12,  -6,  -6,   6,   6,   6,   6,  -6,  -6,  -6,   6,  -6,   6,   6,  -6,   6,  -6,  -8,   8,   8,  -8,  -4,   4,   4,  -4,
 -3,  -3,  -3,  -3,   3,   3,   3,   3,  -4,   4,  -4,   4,  -2,   2,  -2,   2,  -4,  -4,   4,   4,  -2,  -2,   2,   2,  -2,  -2,  -2,  -2,  -1,  -1,  -1,  -1,
  2,   0,   0,   0,  -2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   1,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   2,   0,   0,   0,  -2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
 -6,   6,   0,   0,   6,  -6,   0,   0,  -4,  -2,   0,   0,   4,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -3,   3,   0,   0,  -3,   3,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -2,  -1,   0,   0,  -2,  -1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  4,  -4,   0,   0,  -4,   4,   0,   0,   2,   2,   0,   0,  -2,  -2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   2,  -2,   0,   0,   2,  -2,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   1,   0,   0,   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   2,   0,   0,   0,  -2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  2,   0,   0,   0,  -2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   1,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -6,   6,   0,   0,   6,  -6,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
 -4,  -2,   0,   0,   4,   2,   0,   0,  -3,   3,   0,   0,  -3,   3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -2,  -1,   0,   0,  -2,  -1,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   4,  -4,   0,   0,  -4,   4,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  2,   2,   0,   0,  -2,  -2,   0,   0,   2,  -2,   0,   0,   2,  -2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   1,   0,   0,   1,   1,   0,   0,
 -6,   0,   6,   0,   6,   0,  -6,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -4,   0,  -2,   0,   4,   0,   2,   0,  -3,   0,   3,   0,  -3,   0,   3,   0,
  0,   0,   0,   0,   0,   0,   0,   0,  -2,   0,  -1,   0,  -2,   0,  -1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,  -6,   0,   6,   0,   6,   0,  -6,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
 -4,   0,  -2,   0,   4,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -3,   0,   3,   0,  -3,   0,   3,   0,  -2,   0,  -1,   0,  -2,   0,  -1,   0,
 18, -18, -18,  18, -18,  18,  18, -18,  12,   6, -12,  -6, -12,  -6,  12,   6,  12, -12,   6,  -6, -12,  12,  -6,   6,   9,  -9,  -9,   9,   9,  -9,  -9,   9,
  8,   4,   4,   2,  -8,  -4,  -4,  -2,   6,  -6,   3,  -3,   6,  -6,   3,  -3,   6,   3,  -6,  -3,   6,   3,  -6,  -3,   4,   2,   2,   1,   4,   2,   2,   1,
-12,  12,  12, -12,  12, -12, -12,  12,  -6,  -6,   6,   6,   6,   6,  -6,  -6,  -8,   8,  -4,   4,   8,  -8,   4,  -4,  -6,   6,   6,  -6,  -6,   6,   6,  -6,
 -4,  -4,  -2,  -2,   4,   4,   2,   2,  -4,   4,  -2,   2,  -4,   4,  -2,   2,  -3,  -3,   3,   3,  -3,  -3,   3,   3,  -2,  -2,  -1,  -1,  -2,  -2,  -1,  -1,
  4,   0,  -4,   0,  -4,   0,   4,   0,   0,   0,   0,   0,   0,   0,   0,   0,   2,   0,   2,   0,  -2,   0,  -2,   0,   2,   0,  -2,   0,   2,   0,  -2,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   1,   0,   1,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0,   0,   0,   0,   0,   4,   0,  -4,   0,  -4,   0,   4,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  2,   0,   2,   0,  -2,   0,  -2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   2,   0,  -2,   0,   2,   0,  -2,   0,   1,   0,   1,   0,   1,   0,   1,   0,
-12,  12,  12, -12,  12, -12, -12,  12,  -8,  -4,   8,   4,   8,   4,  -8,  -4,  -6,   6,  -6,   6,   6,  -6,   6,  -6,  -6,   6,   6,  -6,  -6,   6,   6,  -6,
 -4,  -2,  -4,  -2,   4,   2,   4,   2,  -3,   3,  -3,   3,  -3,   3,  -3,   3,  -4,  -2,   4,   2,  -4,  -2,   4,   2,  -2,  -1,  -2,  -1,  -2,  -1,  -2,  -1,
  8,  -8,  -8,   8,  -8,   8,   8,  -8,   4,   4,  -4,  -4,  -4,  -4,   4,   4,   4,  -4,   4,  -4,  -4,   4,  -4,   4,   4,  -4,  -4,   4,   4,  -4,  -4,   4,
  2,   2,   2,   2,  -2,  -2,  -2,  -2,   2,  -2,   2,  -2,   2,  -2,   2,  -2,   2,   2,  -2,  -2,   2,   2,  -2,  -2,   1,   1,   1,   1,   1,   1,   1,   1 };
```