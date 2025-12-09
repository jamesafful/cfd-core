#include "solver.hpp"
#include <iostream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <iomanip>

int main(int argc, char** argv){
  int nx = 400; double cfl = 0.5; int rk = 3; int nsteps = 5; bool deterministic=false;
  for (int i=1;i<argc;i++){
    std::string a(argv[i]);
    auto nextd = [&](double& v){ v = std::atof(argv[++i]); };
    auto nexti = [&](int& v){ v = std::atoi(argv[++i]); };
    if(a=="--nx") nexti(nx);
    else if(a=="--cfl") nextd(cfl);
    else if(a=="--rk") nexti(rk);
    else if(a=="--nsteps") nexti(nsteps);
    else if(a=="--deterministic") deterministic=true;
  }

  GridDesc g{nx,1,1, 1.0/static_cast<double>(nx),1.0,1.0};
  PhysicsParams phys{}; SimParams sim{cfl, rk, nsteps}; BC bc{};
  Solver s(g, phys); s.set_bc(bc);
  RunFlags rf{deterministic};
  s.step(sim, rf);

  // JSON metrics for orchestrator
  std::ostringstream os;
  os << std::fixed << std::setprecision(6)
     << "{"
     << "\"nx\":"<<nx<<","
     << "\"l2\":"<< s.l2_error() << ","
     << "\"kernel_ms\":"<< s.kernel_ms() << ","
     << "\"total_ms\":"<< s.total_ms() << ","
     << "\"deterministic\":"<<(deterministic?1:0)
     << "}";

  std::cout << os.str() << std::endl;
  return 0;
}
