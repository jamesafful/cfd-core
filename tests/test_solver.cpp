#include <gtest/gtest.h>
#include "solver.hpp"

TEST(Solver1D, Sod_L2_Sane){
  GridDesc g{800, 0.0, 1.0, 1.0/800.0};
  PhysicsParams phys{1.4};
  Solver1D s(g, phys, Problem1D::Sod);
  s.set_bc(BC1D{});
  SimParams sp{0.5, 3, 0, 0.2};
  RunFlags rf{true};
  s.run(sp, rf);
  EXPECT_LT(s.l2_density(), 0.2);
}
