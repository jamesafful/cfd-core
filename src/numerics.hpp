#pragma once
#include <algorithm>
#include <cmath>
#include <utility>
#include "state.hpp"

// ---------- MUSCL helper ----------
inline double minmod(double a, double b){
  if (a*b <= 0.0) return 0.0;
  return (std::abs(a) < std::abs(b)) ? a : b;
}

// ---------- Cons <-> Prim ----------
inline Prim cons_to_prim(const Cons& U, double g){
  const double r = std::max(U.r, 1e-14);
  const double u = U.ru / r;
  const double e = std::max((U.E - 0.5 * r * u * u) / r, 0.0);
  const double p = std::max((g - 1.0) * r * e, 1e-14);
  return Prim{r, u, p};
}
inline Cons prim_to_cons(const Prim& W, double g){
  const double E = W.p/(g-1.0) + 0.5 * W.r * W.u * W.u;
  return Cons{W.r, W.r*W.u, E};
}

// ---------- Euler flux ----------
struct Flux { double Fr, Fru, FE; };
inline Flux euler_flux(const Cons& U, double g){
  auto W = cons_to_prim(U,g);
  const double H = (U.E + W.p)/W.r;
  return Flux{ U.ru, U.ru*W.u + W.p, U.ru*H };
}

// ---------- Robust HLL (two-wave) Riemann solver ----------
inline Flux hll(const Cons& UL, const Cons& UR, double g){
  auto WL = cons_to_prim(UL,g);
  auto WR = cons_to_prim(UR,g);
  const double aL = std::sqrt(g*WL.p/WL.r);
  const double aR = std::sqrt(g*WR.p/WR.r);
  // Davis waves (safe & simple)
  const double SL = std::min(WL.u - aL, WR.u - aR);
  const double SR = std::max(WL.u + aL, WR.u + aR);

  const Flux FL = euler_flux(UL,g);
  const Flux FR = euler_flux(UR,g);

  if (SL >= 0.0) return FL;   // supersonic rightward
  if (SR <= 0.0) return FR;   // supersonic leftward

  const double inv = 1.0 / (SR - SL + 1e-14);
  Flux FH;
  FH.Fr  = (SR*FL.Fr  - SL*FR.Fr  + SL*SR*(UR.r  - UL.r )) * inv;
  FH.Fru = (SR*FL.Fru - SL*FR.Fru + SL*SR*(UR.ru - UL.ru)) * inv;
  FH.FE  = (SR*FL.FE  - SL*FR.FE  + SL*SR*(UR.E  - UL.E )) * inv;
  return FH;
}

// ---------- Exact Sod (adequate for L2 checks) ----------
inline Prim sod_exact(double x, double t, double x0, Prim WL, Prim WR, double g){
  if (t <= 0.0) return (x < x0) ? WL : WR;

  auto f = [&](double p, const Prim& W){
    double a = std::sqrt(g*W.p/W.r);
    if (p > W.p){
      double A = 2.0/((g+1.0)*W.r);
      double B = (g-1.0)/(g+1.0)*W.p;
      double ff = (p - W.p) * std::sqrt(A/(p+B));
      double df = std::sqrt(A/(p+B)) * (1.0 - 0.5*(p - W.p)/(p+B));
      return std::pair<double,double>(ff, df);
    } else {
      double pr = std::max(p/W.p, 1e-14);
      double aPow = std::pow(pr, (g-1.0)/(2.0*g));
      double ff = 2.0*a/(g-1.0)*(aPow - 1.0);
      double df = (a/(g*W.p)) * std::pow(pr, -(g+1.0)/(2.0*g));
      return std::pair<double,double>(ff, df);
    }
  };

  double pL=WL.p, pR=WR.p, uL=WL.u, uR=WR.u;
  double p = std::max(1e-6, 0.5*(pL + pR));
  for(int it=0; it<50; ++it){
    auto [fL, dfL] = f(p, WL);
    auto [fR, dfR] = f(p, WR);
    double gF = fL + fR + (uR - uL);
    double gD = dfL + dfR + 1e-14;
    double pnew = p - gF/gD;
    if (std::abs(pnew - p) < 1e-12*(1.0+p)) { p = pnew; break; }
    p = std::max(1e-10, pnew);
  }
  auto [fL2, _d1] = f(p, WL);
  auto [fR2, _d2] = f(p, WR);
  const double u = 0.5*(uL + uR) + 0.5*(fR2 - fL2);

  const double aL = std::sqrt(g*pL/WL.r), aR = std::sqrt(g*pR/WR.r);
  double SL, SHL, SR, SHR;

  if (p > pL){ SL = uL - aL*std::sqrt(1.0 + (g+1.0)/(2.0*g)*(p/pL - 1.0)); SHL=SL; }
  else       { SL = uL - aL; double aS = aL*std::pow(p/pL, (g-1.0)/(2.0*g)); SHL=u-aS; }

  if (p > pR){ SR = uR + aR*std::sqrt(1.0 + (g+1.0)/(2.0*g)*(p/pR - 1.0)); SHR=SR; }
  else       { SR = uR + aR; double aS = aR*std::pow(p/pR, (g-1.0)/(2.0*g)); SHR=u+aS; }

  const double xi = (x - x0)/t;
  if (xi <= SL) return WL;
  if (p <= pL){
    if (xi <= SHL){
      double a = aL + 0.5*(g-1.0)*(xi - uL);
      double r = WL.r * std::pow(a/aL, 2.0/(g-1.0));
      double pxi = pL * std::pow(a/aL, 2.0*g/(g-1.0));
      return Prim{r, xi, pxi};
    } else if (xi <= u){
      double r = WL.r * std::pow(p/pL, 1.0/g);
      return Prim{r, u, p};
    }
  } else {
    if (xi <= u){
      double r = WL.r * ((p/pL + (g-1.0)/(g+1.0)) / ((g-1.0)/(g+1.0)*(p/pL) + 1.0));
      return Prim{r, u, p};
    }
  }
  if (xi <= SHR){
    double r = (p <= pR) ? WR.r * std::pow(p/pR, 1.0/g)
                         : WR.r * ((p/pR + (g-1.0)/(g+1.0)) / ((g-1.0)/(g+1.0)*(p/pR) + 1.0));
    return Prim{r, u, p};
  }
  if (p <= pR){
    if (xi <= SR){
      double a = aR - 0.5*(g-1.0)*(xi - uR);
      double r = WR.r * std::pow(a/aR, 2.0/(g-1.0));
      double pxi = pR * std::pow(a/aR, 2.0*g/(g-1.0));
      return Prim{r, xi, pxi};
    } else return WR;
  } else {
    return WR;
  }
}
