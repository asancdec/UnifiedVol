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

#include <Base/Types.hpp>

#include <concepts>

namespace uv::models::heston::price::detail
{
template <std::floating_point T> struct Integrand
{
    Complex<T> iAlpha;
    Complex<T> onePlusITanPhi;
    Complex<T> c;
    Complex<T> tDivTwo;
    Complex<T> sigmaRho;

    T kappa;
    T kappaThetaDivSigma2;
    T sigma2;

    T v0;
    T t;

    [[gnu::hot]] T operator()(T x) const noexcept;
};

} // namespace uv::models::heston::price::detail

#include "Models/Heston/Price/Detail/Integrand.inl"