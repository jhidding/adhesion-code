// ~/~ begin <<adhesion_example.md#src/initial_conditions/apply_power_spectrum.cc>>[init]
#include "initial_conditions.hh"
#include "fft.hh"

void compute_potential(
    BoxParam const &box,
    std::vector<double> &field,
    PowerSpectrum const &P)
{
  RFFT3 rfft(box);
  std::copy(field.begin(), field.end(),
            rfft.real_space.begin());
  rfft.forward_transform();

  // ~/~ begin <<adhesion_example.md#apply-power-spectrum>>[init]
  auto f_shape = box.rfft_shape();
  std::array<size_t, 3> loc = {0, 0, 1};
  double v = pow(box.L / box.N, 3.0);

  for (size_t i = 1; i < box.rfft_size(); ++i) {
    double k = box.k_abs(loc);
    rfft.fourier_space[i] *= sqrt(P(k) / v) / (k * k);
    increment_index<3>(f_shape, loc);
  }
  // ~/~ end

  rfft.backward_transform();
  std::copy(rfft.real_space.begin(), rfft.real_space.end(),
            field.begin());
}
// ~/~ end