
#pragma once
#include <cmath>
#include <algorithm>

// Placeholder numerics to keep project compiling & runnable in Codespaces.
// Swap with real MUSCL + HLLC and wire into the solver loop later.

inline double dummy_flux(double u_left, double u_right) {
  // Centered flux placeholder
  return 0.5 * (u_left + u_right);
}

inline double minmod(double a, double b){
  if (a*b <= 0.0) return 0.0;
  return (std::abs(a) < std::abs(b)) ? a : b;
}
