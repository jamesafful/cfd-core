
# Agentic CFD — MVP Scaffold (Codespaces-ready)

This repository is a minimal, correctness-first CFD core scaffold with an autotuning orchestrator.
It's designed to run **out of the box on GitHub Codespaces** (CPU-only by default).

## Features
- C++20 core with a tiny 1D structured-grid stepper placeholder (ready to swap in MUSCL+HLLC).
- JSON metrics output (L2 placeholder, kernel/total timing).
- Python orchestrator with Optuna searching **params/build flags**; baseline-normalized speedup; per-trial build dirs.
- Tier-1 CI workflow (CPU): build + unit tests + orchestrator smoke.
- Devcontainer for Codespaces with all dependencies preinstalled (cmake, ninja, Python libs).

> NOTE: The numerics shipped here are **placeholders** to ensure turnkey run in Codespaces.
> Replace `src/numerics.hpp` and wire the solver loop for MUSCL+HLLC + MMS once ready.

## Quickstart (Codespaces)
1. Open in GitHub Codespaces.
2. Build and run:
   ```bash
   bash scripts/build.sh -G Ninja
   ./build/agentic_cfd --nx 400 --cfl 0.5 --rk 3 --nsteps 5 --deterministic
   ```
3. Run unit tests:
   ```bash
   ctest --test-dir build --output-on-failure
   ```
4. Try the autotuner (2 quick trials):
   ```bash
   pip install -r orchestrator/requirements.txt
   KRUNS=2 python orchestrator/study.py --trials 2
   ```

## Structure
```
agentic-cfd/
├─ CMakeLists.txt
├─ src/                # core C++
├─ tests/              # gtest unit tests
├─ orchestrator/       # Python autotuner
├─ scripts/            # helpers
├─ .github/workflows/  # CI
└─ .devcontainer/      # Codespaces config
```

## Next steps
- Swap in **MUSCL + HLLC** for 1D Euler in `src/numerics.hpp` and wire calls in `Solver::step`.
- Add MMS/Sod exact checks; emit real L2/Linf metrics.
- Enable Kokkos and GPU backends (CUDA/HIP) once you're off Codespaces.

License: MIT
