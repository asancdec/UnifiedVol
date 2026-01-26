# Project Structure

```text
UnifiedVol/
├── .clang-format
├── .gitattributes
├── .githooks/
│   ├── pre-commit
├── .gitignore
├── CMakeLists.txt
├── CMakePresets.json
├── LICENSE.txt
├── README.md
├── data/
│   ├── VolSurface_SPY_04072011.csv
├── docs/
│   ├── BUILD.md
│   ├── DATA.md
│   ├── DEPENDENCIES.md
│   ├── TREE.md
│   ├── citations.bib
├── examples/
│   ├── main.cpp
├── include/
│   ├── Core/
│   │   ├── Functions.hpp
│   │   ├── Functions.inl
│   │   ├── MarketData.hpp
│   │   ├── Matrix/
│   │   │   ├── Functions.hpp
│   │   │   ├── Functions.inl
│   │   │   ├── Matrix.hpp
│   │   │   ├── Matrix.inl
│   │   ├── Types.hpp
│   │   ├── VolSurface.hpp
│   │   ├── VolSurface.inl
│   ├── Math/
│   │   ├── Functions.hpp
│   │   ├── Functions.inl
│   │   ├── Integration/
│   │   │   ├── TanHSinH.hpp
│   │   │   ├── TanHSinH.inl
│   │   ├── Interpolation/
│   │   │   ├── Interpolator.hpp
│   │   │   ├── Interpolator.inl
│   │   │   ├── Policies.hpp
│   │   │   ├── Policies.inl
│   │   ├── Optimization/
│   │   │   ├── Ceres/
│   │   │   │   ├── Config.hpp
│   │   │   │   ├── Optimizer.hpp
│   │   │   │   ├── Optimizer.inl
│   │   │   │   ├── Policy.hpp
│   │   │   ├── Functions.hpp
│   │   │   ├── Functions.inl
│   │   │   ├── NLopt/
│   │   │   │   ├── Config.hpp
│   │   │   │   ├── Optimizer.hpp
│   │   │   │   ├── Optimizer.inl
│   │   ├── PDE/
│   │   │   ├── AHCache.hpp
│   │   │   ├── Functions.hpp
│   │   │   ├── Functions.inl
│   ├── Models/
│   │   ├── Heston/
│   │   │   ├── Calibrator.hpp
│   │   │   ├── Calibrator.inl
│   │   │   ├── Config.hpp
│   │   │   ├── Params.hpp
│   │   │   ├── Pricer.hpp
│   │   │   ├── Pricer.inl
│   │   ├── LocalVol/
│   │   │   ├── Calibrator.hpp
│   │   │   ├── Calibrator.inl
│   │   │   ├── Pricer.hpp
│   │   │   ├── Pricer.inl
│   │   │   ├── Surface.hpp
│   │   │   ├── Surface.inl
│   │   │   ├── VarianceView.hpp
│   │   ├── SVI/
│   │   │   ├── Functions.hpp
│   │   │   ├── Functions.inl
│   │   │   ├── Params.hpp
│   ├── Utils/
│   │   ├── Aux/
│   │   │   ├── Errors.hpp
│   │   │   ├── StopWatch.hpp
│   │   │   ├── StopWatch.inl
│   │   ├── IO/
│   │   │   ├── CSV/
│   │   │   │   ├── Functions.hpp
│   │   │   │   ├── Functions.inl
│   │   │   ├── ConsoleRedirect.hpp
│   │   │   ├── Functions.hpp
│   │   │   ├── Functions.inl
│   │   │   ├── Log.hpp
│   │   ├── PCH.hpp
├── src/
│   ├── Math/
│   │   ├── Optimization/
│   │   │   ├── Functions.cpp
│   ├── Models/
│   │   ├── Heston/
│   │   │   ├── Calibrator.cpp
│   │   ├── SVI/
│   │   │   ├── Functions.cpp
│   ├── Utils/
│   │   ├── Aux/
│   │   │   ├── Errors.cpp
│   │   ├── IO/
│   │   │   ├── Log.cpp
│   │   │   ├── StopWatch.cpp
├── tools/
│   ├── gen_tree.sh
├── vcpkg.json
```

