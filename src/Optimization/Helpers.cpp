// SPDX-License-Identifier: Apache-2.0

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

        WARN(
            before < lowerBounds[i] || before > upperBounds[i],
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

        WARN(
            before < bound,
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

        WARN(
            before > bound,
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

            WARN(
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

            WARN(
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
    bool isSuccess,
    std::string_view extraInfo
)
{
    const std::string_view successMessage{isSuccess ? "SUCCESS" : "FAIL"};

    const std::string statusMessage{
        extraInfo.empty() ? std::format("{}", successMessage)
                          : std::format("{} {}", successMessage, extraInfo)
    };

    if (paramNames.empty())
    {

        INFO(std::format(
            "[Calib] SSE={:.8e} ({:.2f} ms, {} it, {}",
            sse,
            elapsedMs,
            iterCount,
            statusMessage
        ));
        return;
    }

    const std::size_t n{x.size()};

    REQUIRE_SAME_SIZE(paramNames, n);

    std::string paramsLine;
    paramsLine.reserve(n * 24);

    for (std::size_t i = 0; i < n; ++i)
    {
        paramsLine += std::format("{}={:.6f}", paramNames[i], x[i]);

        if (i + 1 < n)
        {
            paramsLine += "  ";
        }
    }

    INFO(std::format(
        "[Calib] {}  SSE={:.10e} ({:.2f} ms, {} it, {})",
        paramsLine,
        sse,
        elapsedMs,
        iterCount,
        statusMessage
    ));
}

void validateBounds(
    std::span<const double> x,
    std::span<const double> lowerBounds,
    std::span<const double> upperBounds
)
{
    REQUIRE_NON_EMPTY(x);
    REQUIRE_FINITE(x);

    validateBoundsSpec(x.size(), lowerBounds, upperBounds);
}

void validateBoundsSpec(
    std::size_t n,
    std::span<const double> lowerBounds,
    std::span<const double> upperBounds
)
{
    REQUIRE_SAME_SIZE(lowerBounds, n);
    REQUIRE_SAME_SIZE(upperBounds, n);

    REQUIRE_FINITE(lowerBounds);
    REQUIRE_FINITE(upperBounds);

    REQUIRE_EQUAL_OR_GREATER(upperBounds, lowerBounds);
}

void validateLowerBounds(std::span<const double> x, std::span<const double> lowerBounds)
{
    REQUIRE_NON_EMPTY(x);
    REQUIRE_FINITE(x);
    validateLowerBoundsSpec(x.size(), lowerBounds);
}

void validateLowerBoundsSpec(std::size_t n, std::span<const double> lowerBounds)
{
    REQUIRE_NON_EMPTY(lowerBounds);
    REQUIRE_FINITE(lowerBounds);
    REQUIRE_SAME_SIZE(lowerBounds, n);
}

void validateUpperBounds(std::span<const double> x, std::span<const double> upperBounds)
{
    REQUIRE_NON_EMPTY(x);
    REQUIRE_FINITE(x);
    validateUpperBoundsSpec(x.size(), upperBounds);
}

void validateUpperBoundsSpec(std::size_t n, std::span<const double> upperBounds)
{
    REQUIRE_NON_EMPTY(upperBounds);
    REQUIRE_FINITE(upperBounds);
    REQUIRE_SAME_SIZE(upperBounds, n);
}

} // namespace uv::opt
