#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <tuple>
#include <vector>

struct Cons { double r, ru, E; };
struct Prim { double r, u, p; };

inline Prim cons_to_prim(const Cons& U, double g){
  const double r = std::max(U.r, 1e-14);
  const double u = U.ru / r;
  const double e = (U.E - 0.5 * r * u * u) / r;
  const double p = std::max((g - 1.0) * r * std::max(e, 0.0), 1e-14);
  return Prim{r, u, p};
}

inline Cons prim_to_cons(const Prim& W, double g){
  const double E = W.p/(g-1.0) + 0.5 * W.r * W.u * W.u;
  return Cons{W.r, W.r*W.u, E};
}

inline Cons operator+(const Cons& a, const Cons& b){ return {a.r+b.r, a.ru+b.ru, a.E+b.E}; }
inline Cons operator-(const Cons& a, const Cons& b){ return {a.r-b.r, a.ru-b.ru, a.E-b.E}; }
inline Cons operator*(double s, const Cons& a){ return {s*a.r, s*a.ru, s*a.E}; }

struct Flux { double Fr, Fru, FE; };

inline Flux euler_flux(const Cons& U, double g){
  auto W = cons_to_prim(U,g);
  double H = (U.E + W.p)/W.r;
  return Flux{ U.ru, U.ru*W.u + W.p, U.ru*H };
}

inline double minmod(double a, double b){
  if (a*b <= 0.0) return 0.0;
  return (std::abs(a) < std::abs(b)) ? a : b;
}

// Compact, standard HLLC for 1D Euler (Toro)
inline Flux hllc(const Cons& UL, const Cons& UR, double g){
  auto WL = cons_to_prim(UL,g);
  auto WR = cons_to_prim(UR,g);
  const double aL = std::sqrt(g*WL.p/WL.r);
  const double aR = std::sqrt(g*WR.p/WR.r);

  // Davis' estimates
  const double SL = std::min(WL.u - aL, WR.u - aR);
  const double SR = std::max(WL.u + aL, WR.u + aR);

  const Flux FL = euler_flux(UL,g);
  const Flux FR = euler_flux(UR,g);

  // Degenerate case
  if (SL >= 0.0) return FL;
  if (SR <= 0.0) return FR;

  // Contact speed S*
  const double num = (WR.p - WL.p) + UL.ru*(SL - WL.u) - UR.ru*(SR - WR.u);
  const double den = UL.r*(SL - WL.u) - UR.r*(SR - WR.u) + 1e-14;
  const double Sstar = num / den;

  auto UstarL = [&](Cons U, Prim W){
    const double fac = (SL - W.u)/(SL - Sstar + 1e-14);
    const double rS  = U.r * fac;
    const double ruS = rS * Sstar;
    // Energy star from Rankine–Hugoniot
    const double ES = ( (SL - W.u)*U.E - W.p*W.u + W.p*Sstar ) / (SL - Sstar + 1e-14);
    return Cons{rS, ruS, ES};
  };
  auto UstarR = [&](Cons U, Prim W){
    const double fac = (SR - W.u)/(SR - Sstar + 1e-14);
    const double rS  = U.r * fac;
    const double ruS = rS * Sstar;
    const double ES = ( (SR - W.u)*U.E - W.p*W.u + W.p*Sstar ) / (SR - Sstar + 1e-14);
    return Cons{rS, ruS, ES};
  };

  if (Sstar >= 0.0){
    Cons US = UstarL(UL, WL);
    return Flux{
      FL.Fr  + SL*(US.r  - UL.r ),
      FL.Fru + SL*(US.ru - UL.ru),
      FL.FE  + SL*(US.E  - UL.E )
    };
  } else {
    Cons US = UstarR(UR, WR);
    return Flux{
      FR.Fr  + SR*(US.r  - UR.r ),
      FR.Fru + SR*(US.ru - UR.ru),
      FR.FE  + SR*(US.E  - UR.E )
    };
  }
}

// Exact Sod (Toro) – same as before; used for L2 checks
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
      double pr = p/W.p;
      double ff = 2.0*a/(g-1.0)*(std::pow(pr, (g-1.0)/(2.0*g)) - 1.0);
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
  double SL, STL, SHL, SR, STR, SHR;
  if (p > pL){ SL = uL - aL*std::sqrt(1.0 + (g+1.0)/(2.0*g)*(p/pL - 1.0)); SHL=SL; STL=SL; }
  else       { SL = uL - aL; double aS = aL*std::pow(p/pL, (g-1.0)/(2.0*g)); SHL=u-aS; STL=SL; }
  if (p > pR){ SR = uR + aR*std::sqrt(1.0 + (g+1.0)/(2.0*g)*(p/pR - 1.0)); SHR=SR; STR=SR; }
  else       { SR = uR + aR; double aS = aR*std::pow(p/pR, (g-1.0)/(2.0*g)); SHR=u+aS; STR=SR; }

  double xi = (x - x0)/t;
  if (xi <= SL) return WL;
  if (p <= pL){
    if (xi <= SHL){
      double a = aL + 0.5*(g-1.0)*(xi - uL);
      double r = WL.r * std::pow(a/aL, 2.0/(g-1.0));
      double uxi = xi;
      double pxi = pL * std::pow(a/aL, 2.0*g/(g-1.0));
      return Prim{r, uxi, pxi};
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
      double uxi = xi;
      double pxi = pR * std::pow(a/aR, 2.0*g/(g-1.0));
      return Prim{r, uxi, pxi};
    } else return WR;
  } else {
    return WR;
  }
}
