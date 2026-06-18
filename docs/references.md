# Numerical References

This file maps the main numerical methods in the codebase to entries in
[citations.bib](../citations.bib). It is intended as a quick guide for review and
maintenance, not as a substitute for the source papers.

## Models and Calibration

- SVI parameterization, calibration, and arbitrage constraints:
  `GatheralJacquier2014SVI`, `MartiniMingone2020NoArbitrageSVI`,
  `Ferhati2020SVI`, `Lee2004Moment`

- Heston pricing and calibration:
  `Heston1993SV`, `Albrecher2007LittleHeston`, `Cui2016FastHeston`,
  `AndersenLake2018Fourier`

- Local volatility and volatility interpolation:
  `Dupire1994Smile`, `AndreasenHuge2010Interpolation`

## Numerical Methods

- Black--Scholes implied volatility inversion:
  `Jaeckel2015LetsBeRational`

- Black / Black-76 option pricing:
  `Black1976Commodity`

- Cubic splines and monotone interpolation:
  `DeBoor2001Splines`, `FritschCarlson1980Monotone`,
  `FritschButland1984Monotone`

- Double-exponential / tanh-sinh quadrature:
  `TakahasiMori1974DoubleExponential`

- Numerical linear algebra and floating-point stability:
  `Higham2002Accuracy`
