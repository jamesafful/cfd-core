#pragma once
#include <vector>

struct Cons { double r, ru, E; };   // conservative vars
struct Prim { double r, u, p; };    // primitive vars

struct State1D {
  std::vector<Cons> U;
  std::vector<Cons> Utmp;
};
