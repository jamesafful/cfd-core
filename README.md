
# Agentic CFD — 1D MUSCL+HLLC (Sod) MVP

This repo is **Codespaces-ready** and implements a **1D finite-volume Euler solver** with **MUSCL (minmod)** reconstruction and **HLLC** flux. It adds:

- `--problem sod1d`, `--final_time T` (CFL-based `dt`), `--gamma`.
- Deterministic, conservative update on a uniform grid with transmissive BCs.
- L2/Linf error vs **embedded Sod exact solution**.
- Conserved totals (`mass_total`, `momentum_total`, `energy_total`) and a `physics_hash` in JSON.

## Quickstart

```bash
# Build (CPU)
bash scripts/build.sh -G Ninja    # or omit -G Ninja if Ninja unavailable
# Run Sod (nx=800, t=0.2)
./build/agentic_cfd --problem sod1d --nx 800 --final_time 0.2 --cfl 0.5 --rk 3 --gamma 1.4 --deterministic
# Run tests
ctest --test-dir build --output-on-failure
# Tuner (optional)
pip install -r orchestrator/requirements.txt
KRUNS=2 python orchestrator/study.py --trials 3
```

### Output JSON (example)
```json
{
  "nx": 800,
  "l2": 0.0103,
  "linf": 0.0589,
  "kernel_ms": 4.12,
  "total_ms": 4.17,
  "deterministic": 1,
  "mass_total": 0.56,
  "momentum_total": 0.00,
  "energy_total": 1.80,
  "physics_hash": "e3a1b7a2c1d4e6f8"
}
```

> Note: This is a compact educational implementation. For production, replace arrays with Kokkos Views and extend to 2D/3D.
