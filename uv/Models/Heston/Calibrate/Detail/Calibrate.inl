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
    std::span<const T> maturities,
    std::span<const T> discountFactors,
    std::span<const T> forwards,
    std::span<const T> strikes,
    const core::Matrix<T>& vol,
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
    std::span<const double> maturities,
    std::span<const double> discountFactors,
    std::span<const double> forwards,
    std::span<const double> strikes,
    const core::Matrix<double>& vol,
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

    return calibrate<T, N, Mode, Policy>(
        volSurface.maturities(),
        curve.interpolateDF(maturities),
        volSurface.forwards(),
        volSurface.strikes(),
        volSurface.vol(),
        optimizer,
        weightATM,
        pricer
    );
}

template <
    std::floating_point T,
    std::size_t N,
    opt::ceres::GradientMode Mode,
    typename Policy>
Params<T> calibrate(
    const std::span<const T> maturities,
    const std::span<const T> discountFactors,
    const std::span<const T> forwards,
    const std::span<const T> strikes,
    const core::Matrix<T>& vol,
    opt::ceres::Optimizer<Policy>& optimizer,
    const opt::cost::WeightATM<double>& weightATM,
    price::Pricer<T, N>& pricer
)
{
    if constexpr (std::is_same_v<T, double>)
    {
        return calibrateDouble<T, N, Mode, Policy>(
            maturities,
            discountFactors,
            forwards,
            strikes,
            vol,
            optimizer,
            weightATM,
            pricer
        );
    }

    return calibrateDouble<T, N, Mode, Policy>(
               convertVector<double>(maturities),
               convertVector<double>(discountFactors),
               convertVector<double>(forwards),
               convertVector<double>(strikes),
               vol.template as<double>(),
               optimizer,
               weightATM,
               pricer
    )
        .template as<T>();
}

template <
    std::floating_point CalcT,
    std::size_t N,
    opt::ceres::GradientMode Mode,
    typename Policy>
Params<double> calibrateDouble(
    std::span<const double> maturities,
    std::span<const double> discountFactors,
    std::span<const double> forwards,
    std::span<const double> strikes,
    const core::Matrix<double>& vol,
    opt::ceres::Optimizer<Policy>& optimizer,
    const opt::cost::WeightATM<double>& weightATM,
    price::Pricer<CalcT, N>& pricer
)
{
    setGuessBounds(optimizer);

    optimizer.beginRun();

    Vector<MaturitySlice> slices{
        makeSlices(maturities, discountFactors, forwards, strikes, vol, weightATM)
    };

    for (const auto& s : slices)
    {
        optimizer.addResidualBlock(makeSliceCost<Mode, CalcT, N>(s, pricer));
    }

    return Params<double>{optimizer.solve()};
}
} // namespace uv::models::heston::calibrate::detail
