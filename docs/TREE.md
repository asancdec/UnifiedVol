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
├── external/
│   ├── lets_be_rational/
│   │   ├── .gitignore
│   │   ├── MANIFEST.in
│   │   ├── README.MD
│   │   ├── clean.sh
│   │   ├── lets_be_rational.egg-info/
│   │   │   ├── PKG-INFO
│   │   │   ├── SOURCES.txt
│   │   │   ├── dependency_links.txt
│   │   │   ├── top_level.txt
│   │   ├── lets_be_rational/
│   │   │   ├── LetsBeRational.py
│   │   │   ├── __init__.py
│   │   ├── make.bat
│   │   ├── make_linux.sh
│   │   ├── make_osx.sh
│   │   ├── setup.py
│   │   ├── src/
│   │   │   ├── LetsBeRational.cpp
│   │   │   ├── LetsBeRational.i
│   │   │   ├── LetsBeRational.py
│   │   │   ├── erf_cody.cpp
│   │   │   ├── importexport.h
│   │   │   ├── normaldistribution.cpp
│   │   │   ├── normaldistribution.h
│   │   │   ├── rationalcubic.cpp
│   │   │   ├── rationalcubic.h
│   │   │   ├── version.h
├── include/
│   ├── excluded/
│   │   ├── Math/
│   │   │   ├── Integration/
│   │   │   │   ├── TanHSinH.hpp
│   │   │   │   ├── TanHSinH.inl
│   │   │   ├── Optimization/
│   │   │   │   ├── Ceres/
│   │   │   │   │   ├── Config.hpp
│   │   │   │   │   ├── Optimizer.hpp
│   │   │   │   │   ├── Optimizer.inl
│   │   │   │   │   ├── Policy.hpp
│   │   │   ├── PDE/
│   │   │   │   ├── AHCache.hpp
│   │   │   │   ├── Functions.hpp
│   │   │   │   ├── Functions.inl
│   │   ├── Models/
│   │   │   ├── Heston/
│   │   │   │   ├── Calibrator.hpp
│   │   │   │   ├── Calibrator.inl
│   │   │   │   ├── Config.hpp
│   │   │   │   ├── Params.hpp
│   │   │   │   ├── Pricer.hpp
│   │   │   │   ├── Pricer.inl
│   │   │   ├── LocalVol/
│   │   │   │   ├── Calibrator.hpp
│   │   │   │   ├── Calibrator.inl
│   │   │   │   ├── Pricer.hpp
│   │   │   │   ├── Pricer.inl
│   │   │   │   ├── Surface.hpp
│   │   │   │   ├── Surface.inl
│   │   │   │   ├── VarianceView.hpp
│   │   ├── Utils/
│   │   │   ├── IO/
│   │   │   │   ├── ConsoleRedirect.hpp
│   │   ├── src/
│   │   │   ├── Models/
│   │   │   │   ├── Heston/
│   │   │   │   │   ├── Calibrator.cpp
│   ├── uv/
│   │   ├── Base/
│   │   │   ├── Alias.hpp
│   │   │   ├── Config.hpp
│   │   │   ├── Errors/
│   │   │   │   ├── Detail/
│   │   │   │   │   ├── Validate.inl
│   │   │   │   ├── Errors.hpp
│   │   │   │   ├── Validate.hpp
│   │   │   ├── Macros/
│   │   │   │   ├── DevStatus.hpp
│   │   │   │   ├── Inform.hpp
│   │   │   │   ├── Require.hpp
│   │   │   │   ├── Unreachable.hpp
│   │   │   │   ├── Warn.hpp
│   │   │   ├── Utils/
│   │   │   │   ├── Detail/
│   │   │   │   │   ├── Log.hpp
│   │   │   │   │   ├── StopWatch.inl
│   │   │   │   ├── ScopedTimer.hpp
│   │   │   │   ├── StopWatch.hpp
│   │   ├── Core/
│   │   │   ├── Curve.hpp
│   │   │   ├── Detail/
│   │   │   │   ├── Curve.inl
│   │   │   │   ├── Generate.inl
│   │   │   │   ├── Matrix.inl
│   │   │   │   ├── VolSurface.inl
│   │   │   ├── Generate.hpp
│   │   │   ├── MarketData.hpp
│   │   │   ├── MarketState.hpp
│   │   │   ├── Matrix.hpp
│   │   │   ├── VolSurface.hpp
│   │   ├── IO/
│   │   │   ├── CSV/
│   │   │   │   ├── Detail/
│   │   │   │   │   ├── Read.inl
│   │   │   │   ├── Read.hpp
│   │   │   ├── Detail/
│   │   │   │   ├── Load.inl
│   │   │   │   ├── Print.hpp
│   │   │   │   ├── Print.inl
│   │   │   │   ├── Report.inl
│   │   │   ├── Load.hpp
│   │   │   ├── Report.hpp
│   │   ├── Math/
│   │   │   ├── Functions/
│   │   │   │   ├── Black.hpp
│   │   │   │   ├── Detail/
│   │   │   │   │   ├── Black.inl
│   │   │   │   │   ├── JackelDeclare.hpp
│   │   │   │   │   ├── Primitive.inl
│   │   │   │   │   ├── Volatility.inl
│   │   │   │   ├── Primitive.hpp
│   │   │   │   ├── Volatility.hpp
│   │   │   ├── Integration/
│   │   │   │   ├── Detail/
│   │   │   │   │   ├── TanHSinH.inl
│   │   │   │   ├── TanHSinH.hpp
│   │   │   ├── Interpolation/
│   │   │   │   ├── Detail/
│   │   │   │   │   ├── Interpolator.inl
│   │   │   │   │   ├── Policies.inl
│   │   │   │   ├── Interpolator.hpp
│   │   │   │   ├── Policies.hpp
│   │   │   ├── LinearAlgebra/
│   │   │   │   ├── Detail/
│   │   │   │   │   ├── MatrixOps.inl
│   │   │   │   │   ├── VectorOps.inl
│   │   │   │   ├── MatrixOps.hpp
│   │   │   │   ├── VectorOps.hpp
│   │   ├── Models/
│   │   │   ├── SVI/
│   │   │   │   ├── BuildSurface.hpp
│   │   │   │   ├── Calibrate/
│   │   │   │   │   ├── Calibrate.hpp
│   │   │   │   │   ├── Detail/
│   │   │   │   │   │   ├── Calibrate.inl
│   │   │   │   │   │   ├── Constraints.hpp
│   │   │   │   │   │   ├── Constraints.inl
│   │   │   │   │   │   ├── Contexts.hpp
│   │   │   │   │   │   ├── Contexts.inl
│   │   │   │   │   │   ├── Initialize.hpp
│   │   │   │   │   │   ├── Initialize.inl
│   │   │   │   │   │   ├── Objective.hpp
│   │   │   │   │   │   ├── Objective.inl
│   │   │   │   │   │   ├── SliceData.hpp
│   │   │   │   ├── Config.hpp
│   │   │   │   ├── Detail/
│   │   │   │   │   ├── BuildSurface.inl
│   │   │   │   │   ├── Math.inl
│   │   │   │   │   ├── Params.inl
│   │   │   │   ├── Math.hpp
│   │   │   │   ├── Params.hpp
│   │   ├── Optimization/
│   │   │   ├── Cost.hpp
│   │   │   ├── Detail/
│   │   │   │   ├── Cost.inl
│   │   │   ├── Helpers.hpp
│   │   │   ├── NLopt/
│   │   │   │   ├── Algorithm.hpp
│   │   │   │   ├── Config.hpp
│   │   │   │   ├── Detail/
│   │   │   │   │   ├── MapAlgorithm.hpp
│   │   │   │   │   ├── Optimizer.inl
│   │   │   │   ├── Optimizer.hpp
│   │   ├── UnifiedVol.hpp
├── src/
│   ├── Base/
│   │   ├── Config.cpp
│   │   ├── Errors/
│   │   │   ├── Errors.cpp
│   │   │   ├── Validate.cpp
│   │   ├── Utils/
│   │   │   ├── Detail/
│   │   │   │   ├── Log.cpp
│   │   │   │   ├── StopWatch.cpp
│   ├── IO/
│   │   ├── CSV/
│   │   │   ├── Read.cpp
│   ├── Math/
│   │   ├── Functions/
│   │   │   ├── Volatility.cpp
│   ├── Models/
│   │   ├── SVI/
│   │   │   ├── Calibrate/
│   │   │   │   ├── Detail/
│   │   │   │   │   ├── Constraints.cpp
│   │   │   │   │   ├── Contexts.cpp
│   │   │   │   │   ├── Initialize.cpp
│   │   │   │   │   ├── Objective.cpp
│   │   │   │   │   ├── SliceData.cpp
│   │   │   ├── Math.cpp
│   ├── Optimization/
│   │   ├── Helpers.cpp
├── vcpkg.json
```

