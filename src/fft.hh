// ~\~ language=C++ filename=src/fft.hh
// ~\~ begin <<appendix.md|src/fft.hh>>[init]
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
// ~\~ end
