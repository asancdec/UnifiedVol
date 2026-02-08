// SPDX-License-Identifier: Apache-2.0
/*
 * Copyright (c) 2025 �lvaro S�nchez de Carlos
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

#include "Optimization/Helpers.hpp"
#include "Base/Macros/Inform.hpp"
#include "Base/Macros/Require.hpp"
#include "Base/Macros/Warn.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <format>
#include <string>
#include <string_view>

namespace uv::opt
{
void clampBounds(
    std::span<double> initGuess,
    std::span<const double> lowerBounds,
    std::span<const double> upperBounds,
    bool doValidate
)
{

    if (doValidate)
    {
        validateBounds(initGuess, lowerBounds, upperBounds);
    }

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

void clampLowerBounds(
    std::span<double> initGuess,
    std::span<const double> lowerBounds,
    bool doValidate
)
{
    if (doValidate)
    {
        validateLowerBounds(initGuess, lowerBounds);
    }

    for (std::size_t i = 0; i < initGuess.size(); ++i)
    {
        const double before{initGuess[i]};
        const double bound{lowerBounds[i]};

        const double after{bound >= before ? bound : before};

        UV_WARN(
            after != before,
            std::format(
                "[Calib]: parameter [{}] initial guess = {:.6f} "
                "below lower bound -> clamped to {:.6f} "
                "(lb = {:.6f})",
                i,
                before,
                after,
                bound
            )
        );

        initGuess[i] = after;
    }
}

void clampUpperBounds(
    std::span<double> initGuess,
    std::span<const double> upperBounds,
    bool doValidate
)
{
    if (doValidate)
    {
        validateUpperBounds(initGuess, upperBounds);
    }

    for (std::size_t i = 0; i < initGuess.size(); ++i)
    {
        const double before{initGuess[i]};
        const double bound{upperBounds[i]};

        const double after{bound <= before ? bound : before};

        UV_WARN(
            after != before,
            std::format(
                "[Calib]: parameter [{}] initial guess = {:.6f} "
                "above upper bound -> clamped to {:.6f} "
                "(ub = {:.6f})",
                i,
                before,
                after,
                bound
            )
        );

        initGuess[i] = after;
    }
}

void warnBoundsHit(
    std::span<const double> x,
    const std::optional<std::vector<double>>& lowerBounds,
    const std::optional<std::vector<double>>& upperBounds,
    bool doValidate
)
{
    warnBoundsHit(
        x,
        lowerBounds ? std::span<const double>(*lowerBounds) : std::span<const double>{},
        upperBounds ? std::span<const double>(*upperBounds) : std::span<const double>{},
        doValidate
    );
}

void warnBoundsHit(
    std::span<const double> x,
    std::span<const double> lowerBounds,
    std::span<const double> upperBounds,
    bool doValidate
)
{
    const bool hasLB{!lowerBounds.empty()};
    const bool hasUB{!upperBounds.empty()};

    if (doValidate)
    {
        if (hasLB)
        {
            validateLowerBounds(x, lowerBounds);
        }
        if (hasUB)
        {
            validateUpperBounds(x, upperBounds);
        }
    }

    constexpr double absEps{1e-8};
    constexpr double relEps{1e-8};

    const auto near = [absEps, relEps](double v, double bd) noexcept
    {
        return std::fabs(v - bd) <=
               (absEps + relEps * (std::max)(std::fabs(v), std::fabs(bd)));
    };

    for (std::size_t i = 0; i < x.size(); ++i)
    {
        const double v{x[i]};

        if (hasLB)
        {
            const double lb{lowerBounds[i]};

            UV_WARN(
                near(v, lb),
                std::format(
                    "[Calib]: parameter [{}] hit LOWER bound: v = {:.4f} (lb = {:.4f})",
                    i,
                    v,
                    lb
                )
            );
        }
        if (hasUB)
        {

            const double ub{upperBounds[i]};

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

    if (paramNames.empty())
    {

        UV_INFO(std::format(
            "[Calib] SSE={:.7e} ({:.2f} ms, {} it, {})",
            sse,
            elapsedMs,
            iterCount,
            isSuccess ? "SUCCESS" : "FAIL"
        ));
        return;
    }

    const std::size_t n{x.size()};

    UV_REQUIRE_SAME_SIZE(paramNames, n);

    std::string paramsLine;
    paramsLine.reserve(n * 24);

    for (std::size_t i = 0; i < n; ++i)
    {
        paramsLine += std::format("{}={:.5f}", paramNames[i], x[i]);

        if (i + 1 < n)
        {
            paramsLine += "  ";
        }
    }

    UV_INFO(std::format(
        "[Calib] {}  SSE={:.7e} ({:.2f} ms, {} it, {})",
        paramsLine,
        sse,
        elapsedMs,
        iterCount,
        isSuccess ? "SUCCESS" : "FAIL"
    ));
}

void validateBounds(
    std::span<const double> x,
    std::span<const double> lowerBounds,
    std::span<const double> upperBounds
)
{
    UV_REQUIRE_NON_EMPTY(x);
    UV_REQUIRE_FINITE(x);

    validateBoundsSpec(x.size(), lowerBounds, upperBounds);
}

void validateBoundsSpec(
    std::size_t n,
    std::span<const double> lowerBounds,
    std::span<const double> upperBounds
)
{
    UV_REQUIRE_SAME_SIZE(lowerBounds, n);
    UV_REQUIRE_SAME_SIZE(upperBounds, n);

    UV_REQUIRE_FINITE(lowerBounds);
    UV_REQUIRE_FINITE(upperBounds);

    UV_REQUIRE_EQUAL_OR_GREATER(upperBounds, lowerBounds);
}

void validateLowerBounds(std::span<const double> x, std::span<const double> lowerBounds)
{
    UV_REQUIRE_NON_EMPTY(x);
    UV_REQUIRE_FINITE(x);
    validateLowerBoundsSpec(x.size(), lowerBounds);
}

void validateLowerBoundsSpec(std::size_t n, std::span<const double> lowerBounds)
{
    UV_REQUIRE_NON_EMPTY(lowerBounds);
    UV_REQUIRE_FINITE(lowerBounds);
    UV_REQUIRE_SAME_SIZE(lowerBounds, n);
}

void validateUpperBounds(std::span<const double> x, std::span<const double> upperBounds)
{
    UV_REQUIRE_NON_EMPTY(x);
    UV_REQUIRE_FINITE(x);
    validateUpperBoundsSpec(x.size(), upperBounds);
}

void validateUpperBoundsSpec(std::size_t n, std::span<const double> upperBounds)
{
    UV_REQUIRE_NON_EMPTY(upperBounds);
    UV_REQUIRE_FINITE(upperBounds);
    UV_REQUIRE_SAME_SIZE(upperBounds, n);
}

} // namespace uv::opt
