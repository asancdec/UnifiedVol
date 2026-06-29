// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Base/Types.hpp"
#include "Core/Curve.hpp"
#include "Core/Matrix.hpp"
#include "Core/VolSurface.hpp"
#include "Math/Integration/TanHSinH.hpp"
#include "Models/Heston/Params.hpp"
#include "Models/Heston/Price/Config.hpp"

#include <array>
#include <concepts>
#include <memory>
#include <optional>
#include <span>

namespace uv::models::heston::price
{

template <std::floating_point T, std::size_t N = defaultNodes> class Pricer
{
  private:
    std::optional<Params<T>> params_;
    std::shared_ptr<const math::integration::TanHSinH<T, N>> quad_;

    T alphaItm_;
    T alphaOtm_;

    void validateAlphaDomain() const;
    void validateCallPrice(T t, T dF, T F, T K) const;

    void setAlphas(const Config<T>& config);

  public:
    Pricer();

    explicit Pricer(
        std::shared_ptr<const math::integration::TanHSinH<T, N>> quad,
        const Config<T>& config = {}
    );

    [[gnu::hot]] T
    callPrice(T kappa, T theta, T sigma, T rho, T v0, T t, T dF, T F, T K) const noexcept;

    [[gnu::hot]] T callPrice(T t, T dF, T F, T K, bool doValidate = true) const;

    [[gnu::hot]] void callPrice(
        std::span<T> out,
        T t,
        T dF,
        T F,
        std::span<const T> strikes,
        bool doValidate = true
    ) const;

    core::Matrix<T> callPrice(
        const core::VolSurface<T>& volSurface,
        const core::Curve<T>& curve,
        bool doValidate = true
    ) const;

    [[gnu::hot]] std::array<T, 6>
    callPriceWithGradient(T kappa, T theta, T sigma, T rho, T v0, T t, T dF, T F, T K)
        const noexcept;

    void setParams(const Params<T>& params) noexcept;
};

} // namespace uv::models::heston::price

#include "Models/Heston/Price/Detail/Pricer.inl"
