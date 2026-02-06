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

#include "Core/VolSurface.hpp"
#include "Math/Integration/TanHSinH.hpp"
#include "Models/Heston/Config.hpp"
#include "Models/Heston/Params.hpp"

#include <array>
#include <concepts>
#include <cstddef>
#include <memory>
#include <optional>
#include <tuple>

namespace uv::models::heston
{

template <std::floating_point T, std::size_t N> class Pricer
{
  private:
    struct CharFunCache
    {
        Complex<T> psi;
        Complex<T> A;
        Complex<T> B;
        Complex<T> beta;
        Complex<T> D;
        Complex<T> DT;
        Complex<T> betaPlusD;
        Complex<T> betaMinusD;
        Complex<T> ui;
        Complex<T> kFac;
        T invSigma2;
        T kappaTheta;
        T sigma2;

        Complex<T> uu;
        Complex<T> eDT;
        Complex<T> g;
        Complex<T> Q;
        Complex<T> invQ;
        Complex<T> invQ2;
        Complex<T> R;
        Complex<T> S;
        Complex<T> fracB;
        Complex<T> denomG;
        Complex<T> betaMinusDinvSigma2;
    };

    std::optional<Params<T>> params_;
    std::shared_ptr<const math::integration::TanHSinH<T, N>> quad_;
    const Config<T> config_;

    static T getResidues(T alpha, T F, T K) noexcept;

    T getAlpha(T w) const noexcept;

    static T getPhi(T kappa, T theta, T sigma, T rho, T v0, T t, T w) noexcept;

    static Complex<T> charFunction(
        T kappa,
        T theta,
        T sigma,
        T rho,
        T v0,
        T t,
        const Complex<T>& u
    ) noexcept;

    static CharFunCache charFunctionCal(
        T kappa,
        T theta,
        T sigma,
        T rho,
        T v0,
        T t,
        const Complex<T>& u
    ) noexcept;

  public:
    Pricer() = delete;

    explicit Pricer(
        std::shared_ptr<const math::integration::TanHSinH<T, N>> quad,
        const Config<T>& config
    );

    T callPrice(T kappa, T theta, T sigma, T rho, T v0, T t, T F, T r, T K)
        const noexcept;

    T callPrice(T t, T F, T r, T K) const;

    std::array<T, 6>
    callPriceWithGradient(T kappa, T theta, T sigma, T rho, T v0, T t, T F, T r, T K)
        const noexcept;

    void setParams(const Params<T>& params) noexcept;
};

} // namespace uv::models::heston

#include "Pricer.inl"
