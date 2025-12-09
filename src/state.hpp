#pragma once
#include <vector>

struct Cons { double r, ru, E; };   // conservative vars
struct Prim { double r, u, p; };    // primitive vars

// Basic arithmetic for Cons
inline Cons operator+(const Cons& a, const Cons& b){ return {a.r+b.r, a.ru+b.ru, a.E+b.E}; }
inline Cons operator-(const Cons& a, const Cons& b){ return {a.r-b.r, a.ru-b.ru, a.E-b.E}; }
inline Cons operator*(double s, const Cons& a){ return {s*a.r, s*a.ru, s*a.E}; }
inline Cons operator*(const Cons& a, double s){ return s*a; }

struct State1D {
  std::vector<Cons> U;
  std::vector<Cons> Utmp;
};
