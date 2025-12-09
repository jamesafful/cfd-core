
#include "solver.hpp"
#include <fstream>
#include <cmath>
#include <vector>

Solver::Solver(const GridDesc& g, const PhysicsParams& p) : grid_(g), phys_(p) {
  const int N = grid_.nx * 3; // rho, rhou, E (placeholder usage)
  U_.q.assign(N, 0.0);
  // simple initial condition: a step function in the first variable
  for(int i=0;i<grid_.nx;i++){
    double x = (i + 0.5) * grid_.dx;
    U_.q[i*3 + 0] = (x < 0.5) ? 1.0 : 0.125; // rho-like placeholder
    U_.q[i*3 + 1] = 0.0;
    U_.q[i*3 + 2] = 1.0;
  }
}

void Solver::set_bc(const BC& b){ bc_ = b; }

void Solver::step(const SimParams& sp, const RunFlags& rf){
  auto t_total_start = std::chrono::high_resolution_clock::now();
  auto t_kernel_start = std::chrono::high_resolution_clock::now();

  // Placeholder: do a few smoothing iterations on first component only
  std::vector<double> u(grid_.nx);
  for(int i=0;i<grid_.nx;i++) u[i] = U_.q[i*3];

  for(int n=0;n<sp.nsteps;n++){
    std::vector<double> u_new = u;
    for(int i=1;i<grid_.nx-1;i++){
      double fluxL = dummy_flux(u[i-1], u[i]);
      double fluxR = dummy_flux(u[i], u[i+1]);
      u_new[i] = u[i] - 0.4*(fluxR - fluxL); // toy update
    }
    // periodic BC for placeholder
    u_new[0] = u_new[grid_.nx-2];
    u_new[grid_.nx-1] = u_new[1];
    u.swap(u_new);
  }

  for(int i=0;i<grid_.nx;i++) U_.q[i*3] = u[i];

  auto t_kernel_end = std::chrono::high_resolution_clock::now();
  kernel_ms_ = std::chrono::duration<double, std::milli>(t_kernel_end - t_kernel_start).count();
  total_ms_  = std::chrono::duration<double, std::milli>(t_kernel_end - t_total_start).count();
}

void Solver::write_output(const std::string& path) const{
  std::ofstream f(path);
  f << "nx " << grid_.nx << "\n";
}

double Solver::l2_error() const {
  // Placeholder L2 against initial step (not a real analytic)
  double num = 0.0, den = 0.0;
  for(int i=0;i<grid_.nx;i++){
    double x = (i + 0.5) * grid_.dx;
    double ref = (x < 0.5) ? 1.0 : 0.125;
    double val = U_.q[i*3];
    num += (val - ref)*(val - ref);
    den += ref*ref;
  }
  return std::sqrt(num / std::max(den, 1e-12));
}
