// SPDX-License-Identifier: Apache-2.0

#include "Base/Macros/Unreachable.hpp"
#include "Math/Functions/Black.hpp"
#include "Math/Functions/Volatility.hpp"

#include <span>
#include <utility>

#include <ceres/dynamic_numeric_diff_cost_function.h>

namespace uv::models::heston::calibrate::detail
{

template <std::floating_point T, std::size_t N> struct SliceCommon
{
    const MaturitySlice* s_{nullptr};
    const price::Pricer<T, N>* pricer_{nullptr};

    SliceCommon(const MaturitySlice& s, const price::Pricer<T, N>& pricer) noexcept;

    void residualOnly(const double* p, double* residuals) const;

    void residualAndJac(const double* p, double* residuals, double* jacobian) const;
};

template <std::floating_point T, std::size_t N> class SliceJacobian final
    : public ::ceres::CostFunction
{
  public:
    explicit SliceJacobian(SliceCommon<T, N> common) noexcept;

    bool Evaluate(double const* const* parameters, double* residuals, double** jacobians)
        const override;

  private:
    SliceCommon<T, N> common_;
};

template <std::floating_point T, std::size_t N> struct SliceFunctor
{
    SliceCommon<T, N> common_;

    bool operator()(double const* const* parameters, double* residuals) const;
};

template <std::floating_point T, std::size_t N> SliceCommon<T, N>::SliceCommon(
    const MaturitySlice& s,
    const price::Pricer<T, N>& pricer
) noexcept
    : s_(&s),
      pricer_(&pricer)
{
}

template <std::floating_point T, std::size_t N>
void SliceCommon<T, N>::residualOnly(const double* p, double* residuals) const
{
    const double t{s_->t};
    const double dF{s_->dF};
    const double F{s_->F};
    std::span<const double> strikes{s_->K};

    for (std::size_t i = 0; i < strikes.size(); ++i)
    {

        const double K{strikes[i]};

        const double model{static_cast<double>(
            pricer_->callPrice(p[0], p[1], p[2], p[3], p[4], t, dF, F, K)
        )};

        const double vol{math::vol::impliedVol(model, t, dF, F, K)};

        residuals[i] = (vol - s_->vol[i]) * s_->w[i];
    }
}

template <std::floating_point T, std::size_t N> void
SliceCommon<T, N>::residualAndJac(const double* p, double* residuals, double* J) const
{
    const double t{s_->t};
    const double dF{s_->dF};
    const double F{s_->F};
    std::span<const double> strikes{s_->K};

    for (std::size_t i = 0; i < strikes.size(); ++i)
    {
        const double K{strikes[i]};

        const auto pg = pricer_->callPriceWithGradient(
            p[0],
            p[1],
            p[2],
            p[3],
            p[4],
            s_->t,
            s_->dF,
            s_->F,
            s_->K[i]
        );

        const double vol{math::vol::impliedVol<double>(pg[0], t, dF, F, K)};

        const double wi{s_->w[i]};

        residuals[i] = (vol - s_->vol[i]) * wi;

        const double vegaInv{1.0 / math::black::vegaB76(t, dF, F, vol, K)};

        double* row = &J[i * 5];
        row[0] = pg[1] * wi * vegaInv;
        row[1] = pg[2] * wi * vegaInv;
        row[2] = pg[3] * wi * vegaInv;
        row[3] = pg[4] * wi * vegaInv;
        row[4] = pg[5] * wi * vegaInv;
    }
}

template <std::floating_point T, std::size_t N>
SliceJacobian<T, N>::SliceJacobian(SliceCommon<T, N> c) noexcept
    : common_(std::move(c))
{
    set_num_residuals(static_cast<int>(common_.s_->K.size()));
    mutable_parameter_block_sizes()->push_back(5);
}

template <std::floating_point T, std::size_t N> bool SliceJacobian<T, N>::Evaluate(
    double const* const* parameters,
    double* residuals,
    double** jacobians
) const
{
    const double* p = parameters[0];

    if (!jacobians || !jacobians[0])
    {
        common_.residualOnly(p, residuals);
        return true;
    }

    common_.residualAndJac(p, residuals, jacobians[0]);
    return true;
}

template <std::floating_point T, std::size_t N> bool
SliceFunctor<T, N>::operator()(double const* const* parameters, double* residuals) const
{
    common_.residualOnly(parameters[0], residuals);
    return true;
}

template <std::floating_point T, std::size_t N> std::unique_ptr<::ceres::CostFunction>
makeSliceCostAnalytic(const MaturitySlice& slice, const price::Pricer<T, N>& pricer)
{
    return std::make_unique<SliceJacobian<T, N>>(SliceCommon<T, N>{slice, pricer});
}

template <std::floating_point T, std::size_t N> std::unique_ptr<::ceres::CostFunction>
makeSliceCostNumericForward(const MaturitySlice& slice, const price::Pricer<T, N>& pricer)
{
    using Functor = SliceFunctor<T, N>;
    using Cost = ::ceres::DynamicNumericDiffCostFunction<Functor, ::ceres::FORWARD>;

    auto functor = std::make_unique<Functor>(SliceCommon<T, N>{slice, pricer});
    auto cost = std::make_unique<Cost>(functor.get(), ::ceres::TAKE_OWNERSHIP);
    static_cast<void>(functor.release()); // Ownership transferred to Cost.
    cost->AddParameterBlock(5);
    cost->SetNumResiduals(static_cast<int>(slice.K.size()));
    return cost;
}

template <std::floating_point T, std::size_t N> std::unique_ptr<::ceres::CostFunction>
makeSliceCostNumericCentral(const MaturitySlice& slice, const price::Pricer<T, N>& pricer)
{
    using Functor = SliceFunctor<T, N>;
    using Cost = ::ceres::DynamicNumericDiffCostFunction<Functor, ::ceres::CENTRAL>;

    auto functor = std::make_unique<Functor>(SliceCommon<T, N>{slice, pricer});
    auto cost = std::make_unique<Cost>(functor.get(), ::ceres::TAKE_OWNERSHIP);
    static_cast<void>(functor.release()); // Ownership transferred to Cost.
    cost->AddParameterBlock(5);
    cost->SetNumResiduals(static_cast<int>(slice.K.size()));
    return cost;
}

template <opt::ceres::GradientMode Mode, std::floating_point T, std::size_t N>
std::unique_ptr<::ceres::CostFunction>
makeSliceCost(const MaturitySlice& slice, const price::Pricer<T, N>& pricer)
{
    using enum opt::ceres::GradientMode;

    if constexpr (Mode == Analytic)
    {
        return makeSliceCostAnalytic<T, N>(slice, pricer);
    }
    else if constexpr (Mode == NumericForward)
    {
        return makeSliceCostNumericForward<T, N>(slice, pricer);
    }
    else if constexpr (Mode == NumericCentral)
    {
        return makeSliceCostNumericCentral<T, N>(slice, pricer);
    }
    else
    {
        UNREACHABLE(opt::ceres::GradientMode, Mode);
    }
}

} // namespace uv::models::heston::calibrate::detail
