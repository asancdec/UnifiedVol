[![C++20](https://img.shields.io/badge/C%2B%2B-20-blue.svg)]()
[![License](https://img.shields.io/badge/license-Apache--2.0-lightgrey.svg)]()
[![Build](https://img.shields.io/badge/build-passing-success.svg)]()

# UnifiedVol

**UnifiedVol** is a modern **C++20 quantitative finance library** for volatility surface construction, model calibration, and derivatives pricing.

---

## Design Philosophy

UnifiedVol emphasizes:

- Strong emphasis on accurate replication of at-the-money market quotes during calibration
- Performance-aware design with minimal abstraction overhead
- Efficient memory access and cache-friendly data layouts suitable for production environments
- Numerical stability and reproducibility
- Explicit validation and structured error handling

---

## Currently Supported

### Calibrations
- SVI parametrization
- Dupire Local Volatility (Andreasen–Huge)
- Heston stochastic volatility (Andersen–Lake pricing)

### Numerical Methods
- Finite-difference PDE solvers
- Cubic Hermite interpolation
- PCHIP interpolation
- Sinh–tanh (TanH–SinH) quadrature
- Tridiagonal (Thomas) solvers

### Optimization
- NLopt (constrained optimization)
- Ceres Solver (nonlinear least squares)
- Parameter bounds and tolerance-controlled convergence
- Warm-start calibration support
- Calendar and butterfly arbitrage constraints
- ATM-weighted cost functions
- Calibration to implied volatilites as opposed to raw option prices

### Infrastructure
- Support for extended floating-point precision 
- Structured logging, timing, and error utilities
- Explicit input validation and invariant checks
- Explicit data-lifetime management

---

## Example Usage

Below is a **minimal excerpt** illustrating the structure of a typical calibration pipeline.  
See the `examples/` directory for complete working programs.

```cpp
// -------------- Market data -------------

// Build
const core::MarketState<Real> marketState{io::load::marketState(path, marketData)
};

// Inspect
io::report::volatility(marketState);

// --------------  SVI calibration -------------

// Calibrate
const core::VolSurface<Real> sviVolSurface{models::svi::buildSurface(marketState)
};

// Inspect
io::report::volatility(sviVolSurface);

// --------------  Heston calibration --------------

// Calibrate
const core::VolSurface<Real> hestonVolSurface{
    models::heston::buildSurface<Real>(sviVolSurface, marketState.interestCurve)
};

// Inspect
io::report::volatility(hestonVolSurface);
```

---

## Documentation

- Build instructions: [docs/BUILD.md](docs/BUILD.md)
- File tree: [docs/TREE.md](docs/TREE.md)
- Data sources: [docs/DATA.md](docs/DATA.md)
- Dependencies: [docs/DEPENDENCIES.md](docs/DEPENDENCIES.md)
- Bibliography: [docs/citations.bib](docs/citations.bib)

---

## Disclaimer

This repository is under active development. Interfaces, APIs, and model
implementations are subject to change.

This software is provided for **research and educational purposes only**. It is **not investment advice** and must not be used in trading systems.

All model implementations and market data used in examples and tests have been
derived from **publicly available sources** and are **explicitly cited**
where applicable (see `citations.bib` and `DATA.md`).

This project is **Apache License 2.0 compliant**, and all third-party
dependencies are used and distributed in accordance with their respective
licenses.