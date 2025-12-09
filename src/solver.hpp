#pragma once
#include "grid.hpp"
#include "bc.hpp"
#include "state.hpp"
#include "numerics.hpp"
#include <string>
#include <chrono>
#include <vector>

struct SimParams { double cfl{0.5}; int rk_stages{3}; int nsteps{0}; double final_time{0.2}; };
struct PhysicsParams { double gamma{1.4}; };
struct RunFlags { bool deterministic{false}; };

enum class Problem1D { Sod };

class Solver1D {
 public:
  Solver1D(const GridDesc& g, const PhysicsParams& p, Problem1D prob);
  void set_bc(const BC1D&);
  void run(const SimParams&, const RunFlags&);
  void write_output(const std::string& path) const;
  // metrics
  double l2_density() const { return l2_; }
  double linf_density() const { return linf_; }
  double kernel_ms() const { return kernel_ms_; }
  double total_ms() const { return total_ms_; }
  double mass_total() const { return mass_; }
  double mom_total() const { return mom_; }
  double energy_total() const { return energy_; }
  const GridDesc& grid() const { return grid_; }
  const State1D& state() const { return S_; }

 private:
  GridDesc grid_; PhysicsParams phys_; BC1D bc_{}; Problem1D problem_;
  State1D S_;
  double kernel_ms_{0.0}, total_ms_{0.0};
  double l2_{0.0}, linf_{0.0};
  double mass_{0.0}, mom_{0.0}, energy_{0.0};
  void init_problem();
  void apply_bc(std::vector<Cons>& U) const;
  void step_ssprk3(double dt);
  double max_wavespeed(const std::vector<Cons>& U) const;
  void compute_metrics(double t);
};
