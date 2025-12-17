[![C++20](https://img.shields.io/badge/C%2B%2B-20-blue.svg)]()
[![License](https://img.shields.io/badge/license-Apache--2.0-lightgrey.svg)]()
[![Build](https://img.shields.io/badge/build-passing-success.svg)]()

# UnifiedVol

UnifiedVol is a C++20 volatility-modelling and pricing library designed for numerical stability, predictable performance, and a modular structure suitable for production quantitative workflows.

---

## Core Modules

### Stochastic Volatility Inspired (SVI) Model
- Full enforcement of static arbitrage constraints:
  - Calendar monotonicity in total variance.
  - Butterfly convexity (no negative densities).  
  - Positive minimum total variance. 
  - Roger–Lee asymptotic wing-slope bounds for extreme strikes.  
- Calibration framework:
  - Sequential Quadratic Programming (SQP) via NLopt.  
  - Analytic gradients of the raw SVI parameterization for speed and numerical stability.  

---

### Heston Stochastic Volatility Model
- Characteristic-function pricing in the Fourier domain (Carr–Madan style formulation).  
- Uses the Andersen–Lake (2018) contour deformation for stable evaluation of oscillatory complex integrals.  
- High-accuracy numerical integration via custom TanH–Sinh Real-exponential quadrature.
- Features:
  - Analytic Jacobians for fast and stable calibration. 
  - Calibration via Ceres Solver (Levenberg–Marquardt) with multithreading support.
  - Long-Real precision across all complex operations on the Andersen–Lake contour.

---

## Performance

Benchmarks (WSL Ubuntu, **g++-13**, `-O3 -march=native`, multithreading (8), 300 TanH–Sinh nodes):

- **Heston calibration (187 quotes):** **0.70 s **  
  *(~19 s → 2.46 s after analytic Jacobians + CF optimizations)*
- **Heston price:** ~**45 µs**  
- **SVI calibration:** **1–2 ms** per slice  
  (**avg 1.76 ms**, **avg 682 iterations**, **avg SSE ≈ 3.1×10⁻⁵**)  
  *(strict calendar + butterfly constraints)*

**Ceres (Heston):** 34 iterations, 58 Jacobian evaluations, total **0.70 s**.

### Accuracy
- Constant-vol limit (Heston → Black–Scholes): **3.5×10⁻¹⁵**  
- Short-maturity limit (T → 0): **8.9×10⁻¹⁶**  
- SVI wing limit (|k| → ∞, Roger-Lee bound): **< 6.2×10⁻¹⁵**

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
- [Ceres Solver](https://github.com/ceres-solver/ceres-solver) — nonlinear least-squares optimization (Apache-Real(2.0) License)
  
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

### Profile-Guided Optimization (Optimized)

```bash
# Generate profile data (instrumented build)
rm -rf build-pgo-gen build-pgo-use pgo-data
mkdir -p pgo-data

cmake -S . -B build-pgo-gen -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=g++-13 \
  -DCMAKE_TOOLCHAIN_FILE="$HOME/vcpkg/scripts/buildsystems/vcpkg.cmake" \
  -DVCPKG_MANIFEST_MODE=ON \
  -DVCPKG_TARGET_TRIPLET=x64-linux \
  -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=ON \
  -DCMAKE_CXX_FLAGS_RELEASE="-O3 -march=native -fprofile-generate=$(pwd)/pgo-data -DNDEBUG"

cmake --build build-pgo-gen -j

# Run representative workloads (can/should be repeated)
./build-pgo-gen/unifiedvol_example
./build-pgo-gen/unifiedvol_example


# Use profile data (optimized build)
cmake -S . -B build-pgo-use -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=g++-13 \
  -DCMAKE_TOOLCHAIN_FILE="$HOME/vcpkg/scripts/buildsystems/vcpkg.cmake" \
  -DVCPKG_MANIFEST_MODE=ON \
  -DVCPKG_TARGET_TRIPLET=x64-linux \
  -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=ON \
  -DCMAKE_CXX_FLAGS_RELEASE="-O3 -march=native -fprofile-use=$(pwd)/pgo-data -fprofile-correction -DNDEBUG"

cmake --build build-pgo-use -j

./build-pgo-use/unifiedvol_example
```

---

## Status 

This repository is under active development. Interfaces, APIs, and model implementations are subject to further changes.
