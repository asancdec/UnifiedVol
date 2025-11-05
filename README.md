[![C++20](https://img.shields.io/badge/C%2B%2B-20-blue.svg)]()
[![License](https://img.shields.io/badge/license-Apache--2.0-lightgrey.svg)]()
[![Build](https://img.shields.io/badge/build-passing-success.svg)]()

# UnifiedVol

**UnifiedVol** is a modern **C++20 library** for volatility-surface construction, option pricing, and model calibration.  
Designed for quantitative developers, it emphasizes **numerical precision**, **arbitrage-free consistency**, and **modular extensibility**.

---

## Core Modules

### Arbitrage-Free SVI Calibration
- Full enforcement of static no-arbitrage conditions:  
  - Calendar arbitrage  
  - Butterfly arbitrage  
  - Positive minimum variance  
  - Roger Lee’s asymptotic wing-slope bounds  
- Calibration methodology:  
  - Sequential Quadratic Programming (SQP) via [NLopt](https://nlopt.readthedocs.io/)  
  - Closed-form analytical gradients for improved speed and stability  
  - Post-calibration validation to confirm absence of arbitrage violations  

---

### Heston Stochastic Volatility Model
- Characteristic-function-based pricing for European options using the Fourier-transform approach (Carr–Madan framework).  
- Implements the **Andersen–Lake (2018)** contour-shift formulation for enhanced numerical stability in the complex plane.  
- High-precision numerical integration performed via a custom **TanH–Sinh quadrature** scheme.  
- Supports:
  - Analytical gradient computation for calibration.  
  - Calibration via the Ceres Solver backend.  
  - Long-double precision arithmetic for improved numerical robustness.
    
---

## References

- Roger W. Lee (2003), *The Moment Formula for Implied Volatility at Extreme Strikes*.  
- Jim Gatheral & Antoine Jacquier (2014), *Arbitrage-Free SVI Volatility Surfaces*.  
- Zeliade Systems (2012), *Quasi-Explicit Calibration of Gatheral’s SVI Model* (White Paper).  
- Tahar Ferhati (2020), *Robust Calibration for SVI Model Arbitrage-Free*.  
- Steven L. Heston (1993), *A Closed-Form Solution for Options with Stochastic Volatility with Applications to Bond and Currency Options*, *The Review of Financial Studies*.  
- Hansjörg Albrecher, Philipp Mayer, Wim Schoutens & Jurgen Tistaert (2006), *The Little Heston Trap*.  
- Yiran Cui, Sebastian del Baño Rollin & Guido Germano (2016), *Full and Fast Calibration of the Heston Stochastic Volatility Model*, *arXiv:1511.08718*.  
- Leif Andersen & Mark Lake (2018), *Robust High-Precision Option Pricing by Fourier Transforms: Contour Deformations and Double-Exponential Quadrature*, *Bank of America Merrill Lynch Working Paper*.

---

## Data

Sample surfaces calibrated in this project are derived from publicly available option data:  
- [CBOE Volatility Surface Data](https://datashop.cboe.com/volatility-surfaces)

---

## Dependencies

- [CMake ≥ 3.22](https://cmake.org/download/) — build system (BSD 3-Clause License)    
- [vcpkg](https://github.com/microsoft/vcpkg) (manifest mode enabled) — package manager (MIT License)  
- [NLopt](https://github.com/stevengj/nlopt) — nonlinear optimization (LGPL v2.1)  
- [Ceres Solver](https://github.com/ceres-solver/ceres-solver) — nonlinear least-squares optimization (Apache-2.0 License)
  
---

## Build

### Standard Build (Recommended)

```bash
rm -rf build
cmake -S . -B build -G "Ninja Multi-Config" \
  -DCMAKE_CXX_COMPILER=g++-13 \
  -DCMAKE_TOOLCHAIN_FILE="$HOME/vcpkg/scripts/buildsystems/vcpkg.cmake" \
  -DVCPKG_MANIFEST_MODE=ON \
  -DUNIFIEDVOL_BUILD_EXAMPLE=ON
cmake --build build --config Release
./build/Release/unifiedvol_example
```

### Profile-Guided Optimization (Fun)

```bash
# Generate profiled execution
rm -rf build-pgo pgo-data
cmake -S . -B build-pgo -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=g++-13 \
  -DCMAKE_TOOLCHAIN_FILE="$HOME/vcpkg/scripts/buildsystems/vcpkg.cmake" \
  -DVCPKG_MANIFEST_MODE=ON -DVCPKG_TARGET_TRIPLET=x64-linux \
  -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=ON \
  -DCMAKE_CXX_FLAGS_RELEASE="-O3 -march=native -fprofile-generate=$(pwd)/pgo-data -DNDEBUG"
cmake --build build-pgo -j
./build-pgo/unifiedvol_example

# Use the profiled execution to re-compile optimally
cmake -S . -B build-pgo -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=g++-13 \
  -DCMAKE_TOOLCHAIN_FILE="$HOME/vcpkg/scripts/buildsystems/vcpkg.cmake" \
  -DVCPKG_MANIFEST_MODE=ON -DVCPKG_TARGET_TRIPLET=x64-linux \
  -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=ON \
  -DCMAKE_CXX_FLAGS_RELEASE="-O3 -march=native -fprofile-use=$(pwd)/pgo-data -fprofile-correction -DNDEBUG"
cmake --build build-pgo -j
./build-pgo/unifiedvol_example
```

---

## Status 

This repository is under active development. Interfaces, APIs, and model implementations are subject to further changes.
