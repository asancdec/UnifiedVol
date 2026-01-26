// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Functions.hpp
 * Author:      Álvaro Sánchez de Carlos
 * Created:     2025-01-26
 *
 * Description:
 *   Optimisation helper utilities: ATM weighting, bounds diagnostics,
 *   and logging.
 *
 * Copyright (c) 2025 Álvaro Sánchez de Carlos
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

#include "Utils/Aux/Errors.hpp"
#include "Utils/IO/Log.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <string>

namespace uv::math::opt
{
template <std::floating_point T>
void weightsATM(
    std::span<const T> logKF,
    const WeightATM<T>& params,
    std::span<T> out,
    bool doValidate
)
{
    // ---------- Validate ----------

    if (doValidate)
    {

        detail::validateWeightsATM<T>(logKF, params, out);
    }

    // ---------- Precompute ----------

    const T wATMMinusOne{params.wATM - 1.0};
    const T invk0{1.0 / params.k0};

    // ---------- Calculate ----------

    for (std::size_t i{0}; i < logKF.size(); ++i)
    {
        const T z{logKF[i] * invk0};

        out[i] = std::sqrt(1.0 + wATMMinusOne * std::exp(-(z * z)));
    }
}

template <std::size_t N>
void clamp(
    std::array<double, N>& initGuess,
    const std::array<double, N>& lowerBounds,
    const std::array<double, N>& upperBounds,
    const std::array<std::string_view, N>& paramNames
) noexcept
{
    for (std::size_t i = 0; i < initGuess.size(); ++i)
    {
        const double before = initGuess[i];
        const double after = std::clamp(before, lowerBounds[i], upperBounds[i]);

        UV_WARN(
            after != before,
            std::format(
                "Calibrator: parameter '{}' initial guess = {:.4f} "
                "out of bounds -> clamped to {:.4f} "
                "(lb = {:.4f}, ub = {:.4f})",
                paramNames[i],
                before,
                after,
                lowerBounds[i],
                upperBounds[i]
            )
        );

        initGuess[i] = after;
    }
}

template <std::size_t N>
void warnBoundsHit(
    std::span<double> x,
    const std::array<double, N>& lowerBounds,
    const std::array<double, N>& upperBounds,
    const std::array<std::string_view, N>& paramNames
) noexcept
{
    // Define what touching the bound means
    auto near = [](double v, double bd) noexcept
    {
        constexpr double absEps{1e-12};
        constexpr double relEps{1e-7};
        return std::fabs(v - bd) <=
               (absEps + relEps * (std::max)(std::fabs(v), std::fabs(bd)));
    };

    // Check each parameter
    for (std::size_t i = 0; i < N; ++i)
    {
        const double v{x[i]};
        const double lb{lowerBounds[i]};
        const double ub{upperBounds[i]};

        UV_WARN(
            near(v, lb),
            std::format(
                "Calibrator: parameter '{}' hit LOWER bound: v = "
                "{:.4f} (lb = {:.4f})",
                paramNames[i],
                v,
                lb
            )
        );

        UV_WARN(
            near(v, ub),
            std::format(
                "Calibrator: parameter '{}' hit UPPER bound: v = "
                "{:.4f} (ub = {:.4f})",
                paramNames[i],
                v,
                ub
            )
        );
    }
}

template <std::size_t N>
void logResults(
    std::span<double> x,
    const std::array<std::string_view, N>& paramNames,
    double sse,
    unsigned iterCount,
    double elapsedMs,
    bool isSuccess
) noexcept
{
    std::string paramsLine;
    paramsLine.reserve(N * 24);
    for (std::size_t i = 0; i < N; ++i)
    {
        paramsLine +=
            std::format("{}={:.5f}{}", paramNames[i], x[i], (i + 1 < N ? "  " : ""));
    }

    UV_INFO(std::format(
        "[Calib] {}  SSE={:.7e} ({:.2f} ms, {} it, {})",
        paramsLine, // Parameter values
        sse,        // SSE
        elapsedMs,  // Time in ms
        iterCount,  // Number of iterations
        isSuccess ? "SUCCESS" : "FAIL"
    )); // Success
}
} // namespace uv::math::opt

namespace uv::math::opt::detail
{
template <std::floating_point T>
void validateWeightsATM(
    std::span<const T> logKF,
    const WeightATM<T>& params,
    std::span<T> out
)
{
    // ---------- Extract ----------

    std::size_t numStrikes{logKF.size()};

    const T wATM{params.wATM};
    const T k0{params.k0};

    // ---------- Size ----------

    UV_REQUIRE(
        out.size() == numStrikes,
        ErrorCode::InvalidArgument,
        "validateWeightsATM: size mismatch - out must have the same "
        "size as strikes"
    );

    // ---------- Values ----------

    UV_REQUIRE(
        std::isfinite(wATM) && wATM >= 1.0,
        ErrorCode::InvalidArgument,
        "validateWeightsATM: wATM must be finite and >= 1.0"
    );

    UV_REQUIRE(
        std::isfinite(k0) && k0 > 0.0,
        ErrorCode::InvalidArgument,
        "validateWeightsATM: k0 must be finite and > 0"
    );

    for (std::size_t i{0}; i < numStrikes; ++i)
    {
        UV_REQUIRE(
            std::isfinite(logKF[i]),
            ErrorCode::InvalidArgument,
            "validateWeightsATM: logKF must be finite"
        );
    }
}
} // namespace uv::math::opt::detail