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

#include "Utils/IO/Log.hpp"

#include <bitset>
#include <boost/math/special_functions/lambert_w.hpp>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numbers>
#include <sstream>
#include <utility>

namespace uv::math::integration
{
template <std::floating_point T, std::size_t N>
TanHSinH<T, N>::TanHSinH()
    : h_(

          boost::math::lambert_w0(T(2) * std::numbers::pi_v<T> * T(N)) / T(N)
      )
{

    static_assert(N > 0, "TanHSinH<N>: N must be greater than zero.");
    static_assert(N % 2 == 0, "TanHSinH<N>: N must be even for unroll-by-2 integration.");

    for (unsigned int n = 0; n < N; ++n)
    {
        nodes_[n] = generateNode(T(n) * h_);
    }
}

template <std::floating_point T, std::size_t N>
template <std::size_t M, typename F>
std::array<T, M> TanHSinH<T, N>::integrateZeroToInfMulti(F&& f) const noexcept
{

    constexpr T eps{std::numeric_limits<T>::epsilon()};

    std::array<T, M> sR0{T(0.0)}, sR1{T(0.0)}, sL0{T(0.0)}, sL1{T(0.0)};

    std::bitset<M> actR0, actR1, actL0, actL1;
    actR0.set();
    actR1.set();
    actL0.set();
    actL1.set();

    auto&& func = std::forward<F>(f);

    for (std::size_t i = 0; i + 1 < N; i += 2)
    {

        const bool anyA{actR0.any()};
        const bool anyB{actR1.any()};
        if (!anyA && !anyB)
            break;

        std::array<T, M> ta;
        std::array<T, M> tb;
        T fa{T(0.0)}, fb{T(0.0)};
        if (anyA)
        {
            const Node& a{nodes_[i]};
            ta = func(a.inputRight);
            fa = a.factorRight;
        };
        if (anyB)
        {
            const Node& b{nodes_[i + 1]};
            tb = func(b.inputRight);
            fb = b.factorRight;
        };

        for (std::size_t m = 0; m < M; ++m)
        {
            if (anyA && actR0.test(m))
            {

                T term{fa * ta[m]};
                if (std::fabs(term) <= std::fabs(sR0[m] * eps))
                    actR0.reset(m);
                else
                    sR0[m] += term;
            }
            if (anyB && actR1.test(m))
            {

                T term{fb * tb[m]};
                if (std::fabs(term) <= std::fabs(sR1[m] * eps))
                    actR1.reset(m);
                else
                    sR1[m] += term;
            }
        }
    }

    for (std::size_t i = 1; i + 1 < N; i += 2)
    {

        const bool anyA{actL0.any()};
        const bool anyB{actL1.any()};
        if (!anyA && !anyB)
            break;

        std::array<T, M> ta;
        std::array<T, M> tb;
        T fa{T(0.0)}, fb{T(0.0)};
        if (anyA)
        {
            const Node& a{nodes_[i]};
            ta = func(a.inputLeft);
            fa = a.factorLeft;
        };
        if (anyB)
        {
            const Node& b{nodes_[i + 1]};
            tb = func(b.inputLeft);
            fb = b.factorLeft;
        };

        for (std::size_t m = 0; m < M; ++m)
        {
            if (anyA && actL0.test(m))
            {

                T term{fa * ta[m]};
                if (std::fabs(term) <= std::fabs(sL0[m] * eps))
                    actL0.reset(m);
                else
                    sL0[m] += term;
            }
            if (anyB && actL1.test(m))
            {

                T term{fb * tb[m]};
                if (std::fabs(term) <= std::fabs(sL1[m] * eps))
                    actL1.reset(m);
                else
                    sL1[m] += term;
            }
        }
    }

    std::array<T, M> out{};
    for (std::size_t m = 0; m < M; ++m)
        out[m] = sR0[m] + sR1[m] + sL0[m] + sL1[m];
    return out;
}

template <std::floating_point T, std::size_t N>
template <typename F>
T TanHSinH<T, N>::integrateZeroToInf(F&& f) const noexcept
{

    constexpr T eps{std::numeric_limits<T>::epsilon()};

    T sR0{T(0.0)}, sR1{T(0.0)}, sL0{T(0.0)}, sL1{T(0.0)};

    auto&& func = std::forward<F>(f);

    for (std::size_t i = 0; i + 1 < N; i += 2)
    {

        const Node& a{nodes_[i]};
        const Node& b{nodes_[i + 1]};

        const T ta{a.factorRight * func(a.inputRight)};

        if (std::fabs(ta) <= std::fabs(sR0 * eps))
            break;

        sR0 += ta;

        const T tb{b.factorRight * func(b.inputRight)};

        if (std::fabs(tb) <= std::fabs(sR1 * eps))
            break;

        sR1 += tb;
    }

    for (std::size_t i = 1; i + 1 < N; i += 2)
    {
        const Node& a{nodes_[i]};
        const Node& b{nodes_[i + 1]};

        const T ta{a.factorLeft * func(a.inputLeft)};

        if (std::fabs(ta) <= std::fabs(sL0 * eps))
            break;

        sL0 += ta;

        const T tb{b.factorLeft * func(b.inputLeft)};

        if (std::fabs(tb) <= std::fabs(sL1 * eps))
            break;

        sL1 += tb;
    }

    return sR0 + sR1 + sL0 + sL1;
}

template <std::floating_point T, std::size_t N>
void TanHSinH<T, N>::printGrid() const noexcept
{
    constexpr int idxW = 6;
    constexpr int colW = 24;

    std::ostringstream oss;

    oss << std::left << "\nFixed Tanh-Sinh Grid\n";

    oss << std::right;

    oss << std::setw(colW) << "x_n (node)" << ' ' << std::setw(colW) << "w_n (weight)"
        << '\n';

    oss << std::string(idxW + 1 + colW + 1 + colW, '-') << '\n';

    oss << std::scientific << std::setprecision(16);
    for (std::size_t i = 0; i < N; ++i)
    {
        oss << std::setw(idxW) << std::left << i << std::setw(colW) << std::right
            << nodes_[i].x << ' ' << std::setw(colW) << std::right << nodes_[i].w << '\n';
    }

    UV_INFO(oss.str());
}

template <std::floating_point T, std::size_t N>
TanHSinH<T, N>::Node TanHSinH<T, N>::generateNode(T nh) const noexcept
{

    const T q{std::exp(-std::numbers::pi_v<T> * std::sinh(nh))};

    const T qInv{T(1) / (T(1) + q)};

    const T y{T(2) * q * qInv};

    const T w{qInv * y * std::numbers::pi_v<T> * std::cosh(nh)};

    const T twoMinusY{T(2) - y};

    const T wh{w * h_};

    return {
        w,
        y,
        (T(1) - q) * qInv,
        wh * T(2) / (y * y),
        twoMinusY / y,
        wh * T(2) / (twoMinusY * twoMinusY),
        y / twoMinusY
    };
}
} // namespace uv::math::integration
