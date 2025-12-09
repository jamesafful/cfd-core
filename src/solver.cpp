#include "solver.hpp"
#include <fstream>
#include <cmath>
#include <vector>
#include <numeric>

static inline double clamp(double x, double a, double b){ return std::max(a, std::min(b, x)); }

Solver1D::Solver1D(const GridDesc& g, const PhysicsParams& p, Problem1D prob)
  : grid_(g), phys_(p), problem_(prob) {
  S_.U.assign(grid_.nx, Cons{0,0,0});
  S_.Utmp.assign(grid_.nx, Cons{0,0,0});
  grid_.dx = (grid_.xr - grid_.xl) / static_cast<double>(grid_.nx);
  init_problem();
}

void Solver1D::set_bc(const BC1D& b){ bc_ = b; }

void Solver1D::init_problem(){
  if (problem_ == Problem1D::Sod){
    for(int i=0;i<grid_.nx;i++){
      double xc = grid_.xl + (i + 0.5)*grid_.dx;
      Prim W = (xc < 0.5) ? Prim{1.0, 0.0, 1.0} : Prim{0.125, 0.0, 0.1};
      S_.U[i] = prim_to_cons(W, phys_.gamma);
    }
  }
}

void Solver1D::apply_bc(std::vector<Cons>& U) const{
  (void)U; // transmissive handled by clamped indices
}

double Solver1D::max_wavespeed(const std::vector<Cons>& U) const{
  double amax = 1e-12;
  for(int i=0;i<grid_.nx;i++){
    auto W = cons_to_prim(U[i], phys_.gamma);
    amax = std::max(amax, std::abs(W.u) + std::sqrt(phys_.gamma*W.p/W.r));
  }
  return amax;
}

void Solver1D::step_ssprk3(double dt){
  const int N = grid_.nx;

  auto idx = [&](int j){ return std::min(std::max(j,0), N-1); };

  auto compute_fluxes = [&](const std::vector<Cons>& U, std::vector<Flux>& F){
    std::vector<Prim> W(N);
    for(int i=0;i<N;i++) W[i] = cons_to_prim(U[i], phys_.gamma);

    std::vector<Cons> UL(N+1), UR(N+1);
    for(int i=0;i<=N;i++){
      int il = idx(i-1), ir = idx(i);
      Prim WL = W[il], WR = W[ir];
      auto slope = [&](auto Prim::*m){
        double dl = (W[il].*m) - (W[idx(il-1)].*m);
        double dr = (W[idx(ir+1)].*m) - (W[ir].*m);
        return minmod(dl, dr);
      };
      Prim WLrec{ WL.r + 0.5*slope(&Prim::r),
                  WL.u + 0.5*slope(&Prim::u),
                  WL.p + 0.5*slope(&Prim::p) };
      Prim WRrec{ WR.r - 0.5*slope(&Prim::r),
                  WR.u - 0.5*slope(&Prim::u),
                  WR.p - 0.5*slope(&Prim::p) };
      UL[i] = prim_to_cons(WLrec, phys_.gamma);
      UR[i] = prim_to_cons(WRrec, phys_.gamma);
    }
    F.resize(N+1);
    for(int i=0;i<=N;i++) F[i] = hllc(UL[i], UR[i], phys_.gamma);
  };

  auto conservative_update = [&](std::vector<Cons>& Udst,
                                 const std::vector<Cons>& Usrc,
                                 const std::vector<Flux>& F,
                                 double factor){
    Udst = Usrc;
    const double invdx = 1.0 / grid_.dx;
    for(int i=0;i<N;i++){
      Cons rhs;
      rhs.r  = - (F[i+1].Fr  - F[i].Fr ) * invdx;
      rhs.ru = - (F[i+1].Fru - F[i].Fru) * invdx;
      rhs.E  = - (F[i+1].FE  - F[i].FE ) * invdx;
      Udst[i] = Cons{
        Usrc[i].r  + factor * dt * rhs.r,
        Usrc[i].ru + factor * dt * rhs.ru,
        Usrc[i].E  + factor * dt * rhs.E
      };
    }
  };

  std::vector<Flux> F;
  std::vector<Cons> U1(N), U2(N);

  compute_fluxes(S_.U, F);
  conservative_update(U1, S_.U, F, 1.0);

  compute_fluxes(U1, F);
  conservative_update(U2, U1, F, 0.25);
  for(int i=0;i<N;i++) U2[i] = 0.75*S_.U[i] + 0.25*U2[i];

  compute_fluxes(U2, F);
  conservative_update(S_.Utmp, U2, F, 2.0/3.0);
  for(int i=0;i<N;i++) S_.U[i] = (1.0/3.0)*S_.U[i] + S_.Utmp[i];
}

void Solver1D::compute_metrics(double t){
  mass_ = mom_ = energy_ = 0.0;
  for(int i=0;i<grid_.nx;i++){
    mass_   += S_.U[i].r;
    mom_    += S_.U[i].ru;
    energy_ += S_.U[i].E;
  }
  mass_   *= grid_.dx; mom_ *= grid_.dx; energy_ *= grid_.dx;

  double num=0.0, den=0.0, inf=0.0;
  for(int i=0;i<grid_.nx;i++){
    double xc = grid_.xl + (i + 0.5)*grid_.dx;
    Prim WL{1.0,0.0,1.0}, WR{0.125,0.0,0.1};
    Prim Wex = sod_exact(xc, t, 0.5, WL, WR, phys_.gamma);
    double rho = cons_to_prim(S_.U[i], phys_.gamma).r;
    double d = rho - Wex.r;
    num += d*d; den += Wex.r*Wex.r + 1e-30;
    inf = std::max(inf, std::abs(d));
  }
  l2_ = std::sqrt(num / den);
  linf_ = inf;
}

void Solver1D::run(const SimParams& sp, const RunFlags& rf){
  auto t_total_start = std::chrono::high_resolution_clock::now();
  double t = 0.0;
  while (t < sp.final_time){
    double amax = max_wavespeed(S_.U);
    double dt = sp.cfl * grid_.dx / (amax + 1e-14);
    if (t + dt > sp.final_time) dt = sp.final_time - t;
    auto t_kernel_start = std::chrono::high_resolution_clock::now();
    step_ssprk3(dt);
    auto t_kernel_end = std::chrono::high_resolution_clock::now();
    kernel_ms_ += std::chrono::duration<double, std::milli>(t_kernel_end - t_kernel_start).count();
    t += dt;
  }
  compute_metrics(t);
  auto t_total_end = std::chrono::high_resolution_clock::now();
  total_ms_ = std::chrono::duration<double, std::milli>(t_total_end - t_total_start).count();
}

void Solver1D::write_output(const std::string& path) const{
  std::ofstream f(path);
  f << "nx " << grid_.nx << "\n";
}
