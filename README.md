![C++20](https://img.shields.io/badge/C%2B%2B-20-blue.svg)
[![License: Apache 2.0](https://img.shields.io/badge/License-Apache%202.0-lightgrey.svg)](LICENSE.txt)
[![PRs welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg)](docs/CONTRIBUTING.md)
[![C++ CI](https://github.com/asancdec/UnifiedVol/actions/workflows/cpp-ci.yml/badge.svg?branch=master)](https://github.com/asancdec/UnifiedVol/actions/workflows/cpp-ci.yml)
[![CodeQL](https://github.com/asancdec/UnifiedVol/actions/workflows/codeql.yml/badge.svg)](https://github.com/asancdec/UnifiedVol/actions/workflows/codeql.yml)
[![codecov](https://codecov.io/gh/asancdec/UnifiedVol/branch/master/graph/badge.svg)](https://codecov.io/gh/asancdec/UnifiedVol)


# UnifiedVol

**UnifiedVol** is a **C++20 quantitative finance library** for volatility surface modelling.

The project focuses on front-office style numerical engineering: explicit model
validation, reproducible public-data fixtures, regression coverage, performance
guardrails, and documented references for the implemented methods.

---

## Milestones

- C++20  library with CMake/vcpkg presets for reproducible Linux GCC builds
- SVI and Heston calibration pipelines with NLopt and Ceres Solver integration
- Numerical methods including tanh-sinh integration, Thomas solver, B-splines, PCHIP/Fritsch-Carlson interpolation, non-uniform grids, and root-finding utilities
- Core market-data objects for curves, volatility surfaces, and market states
- Unit, integration, regression, and performance tests with CI and coverage guardrails

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
    models::heston::buildSurface<Real>(sviVolSurface, marketState)
};

// Inspect
io::report::volatility(hestonVolSurface);
```

---

## Documentation

- Contributing: [docs/CONTRIBUTING.md](docs/CONTRIBUTING.md)
- Build instructions: [docs/BUILD.md](docs/BUILD.md)
- File tree: [docs/TREE.md](docs/TREE.md)
- Data sources: [docs/DATA.md](docs/DATA.md)
- Dependencies: [docs/DEPENDENCIES.md](docs/DEPENDENCIES.md)

For model and numerical-method references, see [References](docs/references.md).
Full BibTeX entries are available in [`citations.bib`](citations.bib).

---

## Disclaimer

This repository is under active development. Interfaces, APIs, and model
implementations are subject to change.

This software is provided for **research and educational purposes only**. It is **not investment advice** and must not be used in trading systems.

Model implementations and market data used in examples and tests are documented
where applicable in [docs/references.md](docs/references.md),
[citations.bib](citations.bib), and [docs/DATA.md](docs/DATA.md).

This project is **Apache License 2.0 compliant**, and all third-party
dependencies are used and distributed in accordance with their respective
licenses.
