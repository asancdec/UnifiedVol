# UnifiedVol

UnifiedVol is a work-in-progress quantitative finance library for volatility surface modeling.  
Currently, it provides robust no-arbitrage calibration of the raw SVI parameterization.  

---

## Available Features

### 1. SVI Calibration
- Enforces no-arbitrage conditions:
  1. No calendar arbitrage  
  2. No butterfly spread arbitrage  
  3. Positive minimum variance constraint  
  4. Roger Lee’s asymptotic wing slope bounds  

- Calibration performed with Sequential Quadratic Programming (SQP). 
- Analytical gradients implemented for improved speed and numerical stability.
- Calibration results are tested afterward.

## References
- Roger W. Lee (2003), *The Moment Formula for Implied Volatility at Extreme Strikes*  
- Jim Gatheral & Antoine Jacquier (2014), *Arbitrage-Free SVI Volatility Surfaces*  
- Zeliade Systems (2012), *Quasi-Explicit Calibration of Gatheral’s SVI Model* (White Paper)  
- Tahar Ferhati (2020), *Robust Calibration for SVI Model Arbitrage-Free*
  
## Data
Sample option surface data used for testing has been obtained from the public sample data found in:  

- [CBOE Volatility Surface Data](https://datashop.cboe.com/volatility-surfaces)  

---

## Dependencies
- [NLopt](https://nlopt.readthedocs.io/) — nonlinear optimization library used as the back-end for SQP calibration, constraint handling, and stopping criteria  
- Standard C++20 STL  

---

## Status
This repository is under active development. Interfaces and implementation details are subject to change without notice.  
Future releases are planned to include SSVI, Heston, SABR, and stochastic local volatility (SLV) models.
