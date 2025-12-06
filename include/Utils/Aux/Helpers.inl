/*
* Helpers.hpp
* Author: Alvaro Sanchez de Carlos
*/

#include "Utils/IO/Log.hpp"

#include "Utils/Aux/Errors.hpp"   

#include <format>
#include <algorithm>
#include <string> 
#include <cmath>

namespace uv::utils
{
	template <typename T>
	std::vector<std::vector<T>> transposeMatrix(const std::vector<std::vector<T>>& input)
	{
        // ---------- Sanity checks ------

        // Extract dimensions
        const std::size_t numRows{ input.size() };
        const std::size_t numCols{ input[0].size() };

        // Validate that all rows have the same number of columns
        for (std::size_t i = 1; i < numRows; ++i)
        {
            const std::size_t actual{ input[i].size() };
            const std::size_t expected{ numCols };

            UV_REQUIRE(
                actual == expected,
                ErrorCode::InvalidArgument,
                "transposeVector: inconsistent row length — row " +
                std::to_string(i) + " has " + std::to_string(actual) +
                " elements, expected " + std::to_string(expected)
            );
        }

        // ---------- Transpose matrix ------
        
        // Allocate result with swapped dimensions
        // result[j][i] = input[i][j]
        std::vector<std::vector<T>> transposed(numCols, std::vector<T>(numRows));

        // Perform transpose
        for (std::size_t i = 0; i < numRows; ++i)
        {
            for (std::size_t j = 0; j < numCols; ++j)
            {
                transposed[j][i] = input[i][j];
            }
        }

        return transposed;
	}

    template <std::size_t N>
    void clamp(std::array<double, N>& initGuess,
        const std::array<double, N>& lowerBounds,
        const std::array<double, N>& upperBounds,
        const std::array<std::string_view, N>& paramNames) noexcept
    {
        for (std::size_t i = 0; i < initGuess.size(); ++i)
        {
            const double before = initGuess[i];
            const double after = std::clamp(before, lowerBounds[i], upperBounds[i]);

            UV_WARN(after != before,
                std::format("Calibrator: parameter '{}' initial guess = {:.4f} "
                    "out of bounds -> clamped to {:.4f} "
                    "(lb = {:.4f}, ub = {:.4f})",
                    paramNames[i], before, after,
                    lowerBounds[i], upperBounds[i]));

            initGuess[i] = after;
        }
    }

    template <std::size_t N>
    void warnBoundsHit(std::span<double> x,
        const std::array<double, N>& lowerBounds,
        const std::array<double, N>& upperBounds,
        const std::array<std::string_view, N>& paramNames) noexcept
    {
        // Define what touching the bound means
        auto near = [](double v, double bd) noexcept
            {
                constexpr double absEps{ 1e-12 };
                constexpr double relEps{ 1e-7 };
                return std::fabs(v - bd) <= (absEps + relEps * (std::max)(std::fabs(v), std::fabs(bd)));
            };

        // Check each parameter
        for (std::size_t i = 0; i < N; ++i)
        {
            const double v{ x[i] };
            const double lb{ lowerBounds[i] };
            const double ub{ upperBounds[i] };

            UV_WARN(near(v, lb),
                std::format("Calibrator: parameter '{}' hit LOWER bound: v = {:.4f} (lb = {:.4f})",
                    paramNames[i], v, lb));

            UV_WARN(near(v, ub),
                std::format("Calibrator: parameter '{}' hit UPPER bound: v = {:.4f} (ub = {:.4f})",
                    paramNames[i], v, ub));
        }
    }

    template <std::size_t N>
    void logResults(std::span<double> x,
        const std::array<std::string_view, N>& paramNames,
        double sse,
        unsigned iterCount,
        double elapsedMs,
        bool isSuccess) noexcept
    {
        std::string paramsLine;
        paramsLine.reserve(N * 24);
        for (std::size_t i = 0; i < N; ++i)
        {
            paramsLine += std::format("{}={:.4f}{}",
                paramNames[i],
                x[i],
                (i + 1 < N ? "  " : ""));
        }

        UV_INFO(std::format(
            "[Calib] {}  SSE={:.4e} ({:.2f} ms, {} it, {})",
            paramsLine,                                              // Parameter values
            sse,                                                     // SSE                                               
            elapsedMs,                                               // Time in ms
            iterCount,                                               // Number of iterations            
            isSuccess ? "SUCCESS" : "FAIL"));                        // Success
    }
}