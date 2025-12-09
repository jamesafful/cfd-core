
#pragma once
#include "grid.hpp"
#include "state.hpp"
#include "bc.hpp"
#include "numerics.hpp"
#include <string>
#include <chrono>

struct SimParams { double cfl{0.5}; int rk_stages{3}; int nsteps{10}; };
struct PhysicsParams { double gamma{1.4}; bool viscous{false}; };
struct RunFlags { bool deterministic{false}; };

class Solver {
 public:
  Solver(const GridDesc& g, const PhysicsParams& p);
  void set_bc(const BC&);
  void step(const SimParams&, const RunFlags&);
  void write_output(const std::string& path) const;
  double l2_error() const;
  double kernel_ms() const { return kernel_ms_; }
  double total_ms() const { return total_ms_; }
 private:
  GridDesc grid_; PhysicsParams phys_; BC bc_{}; State U_;
  double kernel_ms_{0.0}, total_ms_{0.0};
};
