
#pragma once
#include <vector>
#include "numerics.hpp"

struct State1D {
  // cell-centered conservative variables
  std::vector<Cons> U;
  // scratch
  std::vector<Cons> Utmp;
};
