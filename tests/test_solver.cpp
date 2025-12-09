
#include <gtest/gtest.h>
#include "solver.hpp"

TEST(Solver1D, Sod_L2_Sane){
  GridDesc g{400, 0.0, 1.0, (1.0-0.0)/400.0};
  PhysicsParams phys{1.4}; BC1D bc{};
  Solver1D s(g, phys, Problem1D::Sod); s.set_bc(bc);
  RunFlags rf{true};
  SimParams sim{0.5, 3, 0, 0.2};
  s.run(sim, rf);
  // sanity: error not huge, conserved totals finite
  EXPECT_TRUE(s.l2_density() >= 0.0);
  EXPECT_TRUE(s.l2_density() < 0.2);
  EXPECT_TRUE(std::isfinite(s.mass_total()));
  EXPECT_TRUE(std::isfinite(s.energy_total()));
}
