[![C++20](https://img.shields.io/badge/C%2B%2B-20-blue.svg)]()
[![License: Apache 2.0](https://img.shields.io/badge/License-Apache%202.0-lightgrey.svg)](LICENSE)
[![C++ CI](https://github.com/asancdec/UnifiedVol/actions/workflows/cpp-ci.yml/badge.svg?branch=master)](https://github.com/asancdec/UnifiedVol/actions/workflows/cpp-ci.yml)

# UnifiedVol

**UnifiedVol** is a **C++20 quantitative finance library** for volatility surface modelling.

---

## Supported Volatility Models

- Arbitrage-free SVI parametrization
- Dupire Local Volatility (Andreasen–Huge)
- Heston stochastic volatility (Andersen–Lake)

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
