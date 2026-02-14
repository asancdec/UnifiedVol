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

#include <concepts>

namespace uv::models::heston::price::detail
{

template <std::floating_point T> struct CharFunCache
{
    Complex<T> psi;
    Complex<T> A;
    Complex<T> B;
    Complex<T> beta;
    Complex<T> D;
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
    Complex<T> R;
    Complex<T> S;
    Complex<T> fracB;
    Complex<T> denomG;
};
} // namespace uv::models::heston::price::detail