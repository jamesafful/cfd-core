
#pragma once
#include <vector>
struct State {
  int nvars{3}; // minimal 1D placeholder: rho, rhou, E (not used in placeholder calc)
  std::vector<double> q; // size = nvars * nx
};
