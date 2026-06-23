// SPDX-License-Identifier: Apache-2.0

#include "Core/Matrix.hpp"
#include "Math/Functions/Black.hpp"
#include "Math/LinearAlgebra/MatrixOps.hpp"
#include "Models/Heston/Calibrate/CeresAdapter.hpp"
#include "Models/Heston/Calibrate/Config.hpp"
#include "Models/Heston/Calibrate/Detail/Initialize.hpp"
#include "Models/Heston/Calibrate/Detail/MaturitySlice.hpp"
#include "Models/Heston/Calibrate/Detail/ResidualCost.hpp"

#include <algorithm>
#include <span>

namespace uv::models::heston::calibrate::detail
{
template <std::floating_point T> struct CalibrationData
{
    std::span<const T> maturities;
    std::span<const T> discountFactors;
    std::span<const T> forwards;
    std::span<const T> strikes;
    const core::Matrix<T>& vol;
};

template <
    std::floating_point T,
    std::size_t N,
    opt::ceres::GradientMode Mode,
    typename Policy>
Params<T> calibrate(
    const core::VolSurface<T>& volSurface,
    const core::Curve<T>& curve,
    opt::ceres::Optimizer<Policy>& optimizer,
    const opt::cost::WeightATM<double>& weightATM,
    price::Pricer<T, N>& pricer
);

template <
    std::floating_point T,
    std::size_t N,
    opt::ceres::GradientMode Mode,
    typename Policy>
Params<T> calibrate(
    const CalibrationData<T>& data,
    opt::ceres::Optimizer<Policy>& optimizer,
    const opt::cost::WeightATM<double>& weightATM,
    price::Pricer<T, N>& pricer
);

template <
    std::floating_point CalcT,
    std::size_t N,
    opt::ceres::GradientMode Mode,
    typename Policy>
Params<double> calibrateDouble(
    const CalibrationData<double>& data,
    opt::ceres::Optimizer<Policy>& optimizer,
    const opt::cost::WeightATM<double>& weightATM,
    price::Pricer<CalcT, N>& pricer
);
} // namespace uv::models::heston::calibrate::detail

namespace uv::models::heston::calibrate
{

template <
    std::floating_point T,
    std::size_t N,
    opt::ceres::GradientMode Mode,
    typename Policy>
Params<T> calibrate(
    const core::VolSurface<T>& volSurface,
    const core::Curve<T>& curve,
    const Config& config
)
{
    price::Pricer<T, N> pricer{};
    opt::ceres::Optimizer<Policy> optimizer{detail::makeOptimizer(config)};

    return detail::calibrate<T, N, Mode, Policy>(
        volSurface,
        curve,
        optimizer,
        config.weightATM,
        pricer
    );
}

template <
    std::floating_point T,
    std::size_t N,
    opt::ceres::GradientMode Mode,
    typename Policy>
Params<T> calibrate(
    const core::VolSurface<T>& volSurface,
    const core::Curve<T>& curve,
    const Config& config,
    price::Pricer<T, N>& pricer
)
{
    opt::ceres::Optimizer<Policy> optimizer{detail::makeOptimizer(config)};

    return detail::calibrate<T, N, Mode, Policy>(
        volSurface,
        curve,
        optimizer,
        config.weightATM,
        pricer
    );
}

} // namespace uv::models::heston::calibrate

namespace uv::models::heston::calibrate::detail
{

template <
    std::floating_point T,
    std::size_t N,
    opt::ceres::GradientMode Mode,
    typename Policy>
Params<T> calibrate(
    const core::VolSurface<T>& volSurface,
    const core::Curve<T>& curve,
    opt::ceres::Optimizer<Policy>& optimizer,
    const opt::cost::WeightATM<double>& weightATM,
    price::Pricer<T, N>& pricer
)
{
    std::span<const T> maturities{volSurface.maturities()};

    const Vector<T> discountFactors{curve.interpolateDF(maturities)};
    const CalibrationData<T> data{
        .maturities = volSurface.maturities(),
        .discountFactors = discountFactors,
        .forwards = volSurface.forwards(),
        .strikes = volSurface.strikes(),
        .vol = volSurface.vol()
    };

    return calibrate<T, N, Mode, Policy>(data, optimizer, weightATM, pricer);
}

template <
    std::floating_point T,
    std::size_t N,
    opt::ceres::GradientMode Mode,
    typename Policy>
Params<T> calibrate(
    const CalibrationData<T>& data,
    opt::ceres::Optimizer<Policy>& optimizer,
    const opt::cost::WeightATM<double>& weightATM,
    price::Pricer<T, N>& pricer
)
{
    if constexpr (std::is_same_v<T, double>)
    {
        return calibrateDouble<T, N, Mode, Policy>(data, optimizer, weightATM, pricer);
    }

    const Vector<double> maturities{convertVector<double>(data.maturities)};
    const Vector<double> discountFactors{convertVector<double>(data.discountFactors)};
    const Vector<double> forwards{convertVector<double>(data.forwards)};
    const Vector<double> strikes{convertVector<double>(data.strikes)};
    const core::Matrix<double> vol{data.vol.template as<double>()};
    const CalibrationData<double> converted{
        .maturities = maturities,
        .discountFactors = discountFactors,
        .forwards = forwards,
        .strikes = strikes,
        .vol = vol
    };

    return calibrateDouble<T, N, Mode, Policy>(converted, optimizer, weightATM, pricer)
        .template as<T>();
}

template <
    std::floating_point CalcT,
    std::size_t N,
    opt::ceres::GradientMode Mode,
    typename Policy>
Params<double> calibrateDouble(
    const CalibrationData<double>& data,
    opt::ceres::Optimizer<Policy>& optimizer,
    const opt::cost::WeightATM<double>& weightATM,
    price::Pricer<CalcT, N>& pricer
)
{
    setGuessBounds(optimizer);

    optimizer.beginRun();

    Vector<MaturitySlice> slices{makeSlices(
        data.maturities,
        data.discountFactors,
        data.forwards,
        data.strikes,
        data.vol,
        weightATM
    )};

    for (const auto& s : slices)
    {
        optimizer.addResidualBlock(makeSliceCost<Mode, CalcT, N>(s, pricer));
    }

    return Params<double>{optimizer.solve()};
}
} // namespace uv::models::heston::calibrate::detail
