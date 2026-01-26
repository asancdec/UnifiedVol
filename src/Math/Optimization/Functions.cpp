// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Functions.hpp
 * Author:      Álvaro Sánchez de Carlos
 * Created:     2025-01-26
 *
 * Description:
 *   Optimization helper utilities.
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

#include "Math/Optimization/Functions.hpp"
#include "Utils/Aux/Errors.hpp"
#include "Utils/IO/Log.hpp"

#include <string_view>
#include <string>
#include <format>
#include <algorithm>
#include <cstddef>
#include <cmath>

namespace uv::math::opt
{
void clamp(
    std::span<double> initGuess,
    std::span<const double> lowerBounds,
    std::span<const double> upperBounds,
    bool doValidate
)
{
    // ---------- Validate ----------

    if (doValidate)
    {
        detail::validateBounds(initGuess, lowerBounds, upperBounds);
    }

    // ---------- Clamp ----------

    for (std::size_t i = 0; i < initGuess.size(); ++i)
    {
        const double before{initGuess[i]};
        const double after{std::clamp(before, lowerBounds[i], upperBounds[i])};

        UV_WARN(
            after != before,
            std::format(
                "[Calib]: parameter [{}] initial guess = {:.4f} "
                "out of bounds -> clamped to {:.4f} "
                "(lb = {:.4f}, ub = {:.4f})",
                i,
                before,
                after,
                lowerBounds[i],
                upperBounds[i]
            )
        );

        initGuess[i] = after;
    }
}
void warnBoundsHit(
    std::span<const double> x,
    std::span<const double> lowerBounds,
    std::span<const double> upperBounds,
    bool doValidate
)
{
    // ---------- Validate ----------

    if (doValidate)
    {
        detail::validateBounds(x, lowerBounds, upperBounds);
    }

    // ---------- Criteria ----------

    constexpr double absEps{1e-8};
    constexpr double relEps{1e-8};

    // Lambda to evaluate if near the bounds
    const auto near = [absEps, relEps](double v, double bd) noexcept
    {
        return std::fabs(v - bd) <=
               (absEps + relEps * (std::max)(std::fabs(v), std::fabs(bd)));
    };

    // ---------- Check ----------

    for (std::size_t i = 0; i < x.size(); ++i)
    {
        const double v{x[i]};
        const double lb{lowerBounds[i]};
        const double ub{upperBounds[i]};

        UV_WARN(
            near(v, lb),
            std::format(
                "[Calib]: parameter [{}] hit LOWER bound: v = {:.4f} (lb = {:.4f})",
                i,
                v,
                lb
            )
        );

        UV_WARN(
            near(v, ub),
            std::format(
                "[Calib]: parameter [{}] hit UPPER bound: v = {:.4f} (ub = {:.4f})",
                i,
                v,
                ub
            )
        );
    }
}

void logResults(
    std::span<const double> x,
    std::span<const std::string_view> paramNames,
    double sse,
    unsigned iterCount,
    double elapsedMs,
    bool isSuccess
)
{
    // ---------- Check ----------

    if (paramNames.empty())
    {
        // ---------- Log ----------

        UV_INFO(std::format(
            "[Calib] SSE={:.7e} ({:.2f} ms, {} it, {})",
            sse,
            elapsedMs,
            iterCount,
            isSuccess ? "SUCCESS" : "FAIL"
        )); 
        return;
    }

    // ---------- Validate ----------

    const std::size_t n{x.size()};

    UV_REQUIRE(
        paramNames.size() == n,
        ErrorCode::InvalidArgument,
        std::format(
            "logResults: paramNames.size() = {} does not match x.size() = {}",
            paramNames.size(),
            n
        )
    );

    // ---------- Allocate ----------

    std::string paramsLine;
    paramsLine.reserve(n * 24);

    // ---------- Fill ----------

    for (std::size_t i = 0; i < n; ++i)
    {
        paramsLine += std::format("{}={:.5f}", paramNames[i], x[i]);

        if (i + 1 < n)
        {
            paramsLine += "  ";
        }
    }

    // ---------- Log ----------

    UV_INFO(std::format(
        "[Calib] {}  SSE={:.7e} ({:.2f} ms, {} it, {})",
        paramsLine,
        sse,
        elapsedMs,
        iterCount,
        isSuccess ? "SUCCESS" : "FAIL"
    ));
}
} // namespace uv::math::opt

namespace uv::math::opt::detail
{
    using ErrorCode::InvalidArgument;

    void validateBounds(
        std::span<const double> x,
        std::span<const double> lowerBounds,
        std::span<const double> upperBounds
    )
    {
        // ---------- Size ----------

        const std::size_t n{x.size()};

        UV_REQUIRE(
            n == lowerBounds.size(),
            InvalidArgument,
            std::format(
                "validateBounds: size mismatch: x.size() = {}, lowerBounds.size() = {}",
                n,
                lowerBounds.size()
            )
        );

        UV_REQUIRE(
            n == upperBounds.size(),
            InvalidArgument,
            std::format(
                "validateBounds: size mismatch: x.size() = {}, upperBounds.size() = {}",
                n,
                upperBounds.size()
            )
        );

        // ---------- Values ----------

        for (std::size_t i = 0; i < n; ++i)
        {
            UV_REQUIRE(
                std::isfinite(x[i]),
                InvalidArgument,
                std::format("validateBounds: x[{}] is not finite (value = {})", i, x[i])
            );

            UV_REQUIRE(
                std::isfinite(lowerBounds[i]) && std::isfinite(upperBounds[i]),
                InvalidArgument,
                std::format(
                    "validateBounds: bounds[{}] not finite (lb = {}, ub = {})",
                    i,
                    lowerBounds[i],
                    upperBounds[i]
                )
            );

            UV_REQUIRE(
                lowerBounds[i] <= upperBounds[i],
                InvalidArgument,
                std::format(
                    "validateBounds: invalid bounds[{}] (lb = {}, ub = {})",
                    i,
                    lowerBounds[i],
                    upperBounds[i]
                )
            );
        }
    }

} // namespace uv::math::opt::detail