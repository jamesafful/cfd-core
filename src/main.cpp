#include "solver.hpp"
#include <iostream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <iomanip>

// Simple 64-bit FNV-1a hash -> 16-hex char string (portable, no OpenSSL)
static std::string fnv1a16(const std::string& s){
  const uint64_t FNV_OFFSET = 14695981039346656037ull;
  const uint64_t FNV_PRIME  = 1099511628211ull;
  uint64_t h = FNV_OFFSET;
  for(unsigned char c : s){ h ^= (uint64_t)c; h *= FNV_PRIME; }
  static const char* hex = "0123456789abcdef";
  std::string out; out.resize(16);
  for(int i=0;i<16;i++){
    int shift = (60 - 4*i);
    unsigned v = (unsigned)((h >> (shift < 0 ? 0 : shift)) & 0xF);
    out[i] = hex[v];
  }
  return out;
}

int main(int argc, char** argv){
  int nx = 800; double cfl = 0.5; int rk = 3; double final_time = 0.2; bool deterministic=false;
  double xl=0.0, xr=1.0; double gamma=1.4;
  std::string problem = "sod1d";

  for (int i=1;i<argc;i++){
    std::string a(argv[i]);
    auto nextd = [&](double& v){ if(i+1<argc){ v = std::atof(argv[++i]); } };
    auto nexti = [&](int& v){ if(i+1<argc){ v = std::atoi(argv[++i]); } };
    if(a=="--nx") nexti(nx);
    else if(a=="--cfl") nextd(cfl);
    else if(a=="--rk") nexti(rk);
    else if(a=="--final_time") nextd(final_time);
    else if(a=="--xl") nextd(xl);
    else if(a=="--xr") nextd(xr);
    else if(a=="--gamma") nextd(gamma);
    else if(a=="--problem" && i+1<argc) problem = argv[++i];
    else if(a=="--deterministic") deterministic=true;
  }

  GridDesc g{nx, xl, xr, (xr-xl)/static_cast<double>(nx)};
  PhysicsParams phys{gamma}; BC1D bc{};
  Solver1D s(g, phys, Problem1D::Sod); s.set_bc(bc);
  RunFlags rf{deterministic};
  SimParams sim{cfl, rk, 0, final_time};
  s.run(sim, rf);

  std::ostringstream key;
  key << std::fixed << std::setprecision(7)
      << s.l2_density() << "," << s.linf_density()
      << std::setprecision(10)
      << "," << s.mass_total() << "," << s.mom_total() << "," << s.energy_total();
  std::string phash = fnv1a16(key.str());

  std::ostringstream os;
  os << std::fixed << std::setprecision(6)
     << "{"
     << "\"nx\":"<<nx<<","
     << "\"l2\":"<< s.l2_density() << ","
     << "\"linf\":"<< s.linf_density() << ","
     << "\"kernel_ms\":"<< s.kernel_ms() << ","
     << "\"total_ms\":"<< s.total_ms() << ","
     << "\"deterministic\":"<<(deterministic?1:0) << ","
     << "\"mass_total\":"<< s.mass_total() << ","
     << "\"momentum_total\":"<< s.mom_total() << ","
     << "\"energy_total\":"<< s.energy_total() << ","
     << "\"physics_hash\":\""<< phash << "\""
     << "}";
  std::cout << os.str() << std::endl;
  return 0;
}
