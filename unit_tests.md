# Testing

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
