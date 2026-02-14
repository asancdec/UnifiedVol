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
#include "Models/Heston/Price/Detail/CharFunction.hpp"

#include <cmath>
#include <complex>

namespace uv::models::heston::price::detail
{
template <std::floating_point T> T Integrand<T>::operator()(T x) const noexcept
{
    constexpr Complex<T> i{T{0}, T{1}};

    const Complex<T> h{iAlpha + x * onePlusITanPhi};

    const Complex<T> hMinusI{h - i};

    const Complex<T> psi{charFunction(
        kappa,
        kappaThetaDivSigma2,
        sigma2,
        v0,
        t,
        tDivTwo,
        sigmaRho,
        hMinusI
    )};

    return std::real(std::exp(x * c) * psi / (hMinusI * h) * onePlusITanPhi);
};

} // namespace uv::models::heston::price::detail