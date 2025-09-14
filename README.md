# UnifiedVol  

**UnifiedVol** is a C++20 quantitative finance library for volatility surface modeling and calibration.  

---

## Current Capabilities  

### Arbitrage-Free SVI Calibration  
- Full enforcement of static no-arbitrage conditions:  
  - Calendar arbitrage  
  - Butterfly spread arbitrage  
  - Positive minimum variance  
  - Roger Lee’s asymptotic wing slope bounds  

- Calibration methodology:  
  - Sequential Quadratic Programming (SQP) with constraint handling via [NLopt](https://nlopt.readthedocs.io/)  
  - Closed-form analytical gradients for improved speed and numerical stability  
  - Post-calibration validation to confirm absence of arbitrage violations  

---

## Examples

- [SVI Calibration Fits (PDF)](docs/SVI_calibration_fits.pdf)  
- Market implied volatilities vs fitted SVI curves (multi-tenor surface)  :

<img width="1867" height="1157" alt="SVI calibration" src="https://github.com/user-attachments/assets/ed85c8bf-f065-494e-b6ba-d7b19175c5f2" />


---

## Roadmap  

Planned extensions include:  
- SSVI (Stochastic Volatility Inspired) parameterization  
- SABR stochastic volatility model  
- Heston stochastic variance model  
- Stochastic Local Volatility (SLV) frameworks  
- Integration with PDE solvers for local volatility pricing  

---

## References  

- Roger W. Lee (2003), *The Moment Formula for Implied Volatility at Extreme Strikes*  
- Jim Gatheral & Antoine Jacquier (2014), *Arbitrage-Free SVI Volatility Surfaces*  
- Zeliade Systems (2012), *Quasi-Explicit Calibration of Gatheral’s SVI Model* (White Paper)  
- Tahar Ferhati (2020), *Robust Calibration for SVI Model Arbitrage-Free*  

---

## Data  

Sample surfaces calibrated in this project are derived from publicly available option data:  

- [CBOE Volatility Surface Data](https://datashop.cboe.com/volatility-surfaces)  

---

## Dependencies  

- C++20 STL  
- [NLopt](https://nlopt.readthedocs.io/) — nonlinear optimization library for SQP calibration, constraint enforcement, and stopping criteria  

---

## Status  

This repository is under active development. Interfaces, APIs, and model implementations are subject to further changes.
