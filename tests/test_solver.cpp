
#include <gtest/gtest.h>
#include "solver.hpp"

TEST(Solver, Smoke){
  GridDesc g{64,1,1, 1.0/64.0,1.0,1.0};
  PhysicsParams phys{}; SimParams sim{0.5,3,2}; BC bc{}; RunFlags rf{true};
  Solver s(g, phys); s.set_bc(bc); s.step(sim, rf);
  // L2 is finite and non-negative
  EXPECT_GE(s.l2_error(), 0.0);
  EXPECT_LT(s.l2_error(), 10.0);
}
