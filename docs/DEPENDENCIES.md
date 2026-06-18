# Dependencies

This project relies on the following direct external tools and libraries.
Transitive dependencies are resolved by vcpkg and by the upstream packages.

## Build System

- **CMake >= 3.22**  
  https://cmake.org  
  License: BSD 3-Clause  
  Citation: `CMake` in [citations.bib](../citations.bib)

- **vcpkg** (manifest mode)  
  https://github.com/microsoft/vcpkg  
  License: MIT  
  Citation: `Vcpkg` in [citations.bib](../citations.bib)

The package manifest is [vcpkg.json](../vcpkg.json).

## Mathematical Libraries

- **Boost.Math** — Special mathematical functions  
  https://www.boost.org/doc/libs/release/libs/math/doc/html/index.html  
  License: Boost Software License 1.0  
  Used for special functions in numerical integration code.  
  Citation: `BoostMath` in [citations.bib](../citations.bib)

## Optimization Libraries

- **NLopt** — Nonlinear optimization  
  https://github.com/stevengj/nlopt  
  License: GNU LGPL  
  Citation: `JohnsonNLopt` in [citations.bib](../citations.bib)

- **Ceres Solver** — Nonlinear least-squares optimization  
  https://github.com/ceres-solver/ceres-solver  
  License: See upstream LICENSE file  
  Citation: `Agarwal2023Ceres` in [citations.bib](../citations.bib)

## Implied Volatility

- **Let’s Be Rational** — Robust Black–Scholes implied volatility  
  Author: Peter Jäckel  
  https://github.com/vollib/lets_be_rational  
  License: Custom permissive license; copyright notice must be preserved  
  Citation: `Jaeckel2015LetsBeRational` in [citations.bib](../citations.bib)

## Test Dependencies

- **GoogleTest** — C++ unit, integration, regression, and performance tests  
  https://github.com/google/googletest  
  License: BSD 3-Clause  
  Used only by the test targets under `tests/`.  
  Citation: `GoogleTest` in [citations.bib](../citations.bib)

## Notes

This project includes the *Let’s Be Rational* implementation by Peter Jäckel,
used to compute Black–Scholes implied volatility with full double-precision
accuracy and guaranteed convergence.

The library is included as a git submodule and built as part of this project.
The original copyright notice is preserved in accordance with the license.

All third-party libraries remain the property of their respective authors.
Users are responsible for complying with the corresponding license terms.
