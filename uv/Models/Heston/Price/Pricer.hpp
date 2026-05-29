// SPDX-License-Identifier: Apache-2.0
/*
 * Copyright (c) 2025 Alvaro Sanchez de Carlos
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at:
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under this License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the LICENSE for the specific language governing permissions and
 * limitations under this License.
 */

#pragma once

#include "Base/Types.hpp"
#include "Core/Curve.hpp"
#include "Core/Matrix.hpp"
#include "Core/VolSurface.hpp"
#include "Math/Integration/TanHSinH.hpp"
#include "Models/Heston/Params.hpp"
#include "Models/Heston/Price/Config.hpp"
#include "Models/Heston/Price/Detail/CharFunction.hpp"

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

    void validateAlphaDomain_() const;
    void validateCallPrice_(T t, T dF, T F, T K) const;

    void setAlphas_(const Config<T>& config);

    [[gnu::hot]] T getAlpha_(T w) const noexcept;

    [[gnu::hot]] static T getResidues_(T alpha, T F, T K) noexcept;

    [[gnu::hot]] static T
    getPhi_(T kappa, T theta, T sigma, T rho, T v0, T t, T w) noexcept;

  public:
    Pricer();

    Pricer(
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
