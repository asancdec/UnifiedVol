# Numerical References

This file maps the main numerical methods in the codebase to references in
[citations.bib](citations.bib). Use it as a quick guide when reviewing or
extending model, interpolation, integration, optimization, and pricing code.

## Volatility Models

- SVI calibration and static-arbitrage constraints:
  `GatheralJacquier2014SVI`, `MartiniMingone2020NoArbitrageSVI`,
  `Ferhati2020SVI`, `Lee2004Moment`
- Heston pricing and calibration:
  `Heston1993SV`, `Albrecher2006LittleHeston`, `Cui2016FastHeston`,
  `AndersenLake2018Fourier`
- Local volatility and volatility interpolation:
  `Dupire1994Smile`, `AndreasenHuge2011Interpolation`

## Pricing And Implied Volatility

- Black-76 forward option pricing:
  `Black1976Commodity`
- Black--Scholes implied volatility inversion:
  `Jaeckel2015LetsBeRational`

## Interpolation

- B-spline basis and de Boor evaluation:
  `DeBoor2001Splines`
- Monotone cubic interpolation:
  `FritschCarlson1980Monotone`, `FritschButland1984Monotone`

## Integration And Linear Algebra

- Tanh-sinh / double-exponential quadrature:
  `TakahasiMori1974DoubleExponential`, `AndersenLake2018Fourier`
- Boost.Math special functions used by integration code:
  `BoostMath`
- Tridiagonal linear solves and stability background:
  `Higham2002Accuracy`

## Optimization Libraries

- Ceres Solver:
  `Agarwal2023Ceres`
- NLopt:
  `JohnsonNLopt`
