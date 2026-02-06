# Dependencies

This project relies on the following external tools and libraries.

## Build System

- **CMake >= 3.22**  
  https://cmake.org  
  License: BSD 3-Clause

- **vcpkg** (manifest mode)  
  https://github.com/microsoft/vcpkg  
  License: MIT

## Optimization Libraries

- **NLopt** — Nonlinear optimization  
  https://github.com/stevengj/nlopt  
  License: LGPL v2.1

- **Ceres Solver** — Nonlinear least-squares optimization  
  https://github.com/ceres-solver/ceres-solver  
  License: Apache License 2.0

## Implied Volatility

- **Let’s Be Rational** — Robust Black–Scholes implied volatility  
  Author: Peter Jäckel  
  https://github.com/vollib/lets_be_rational  
  License: Permissive (custom, attribution required)

### Notes
This project includes the *Let’s Be Rational* implementation by Peter Jäckel,
used to compute Black–Scholes implied volatility with full double-precision
accuracy and guaranteed convergence.

The library is included as a git submodule and built as part of this project.
The original copyright notice is preserved in accordance with the license.

## Notes

All third-party libraries remain the property of their respective authors.
Users are responsible for complying with the corresponding license terms.
