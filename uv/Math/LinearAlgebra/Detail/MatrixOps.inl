// SPDX-License-Identifier: Apache-2.0

#include "Base/Macros/Require.hpp"
#include "Core/Matrix.hpp"
#include "Math/LinearAlgebra/VectorOps.hpp"

#include <cmath>
#include <concepts>
#include <cstddef>
#include <functional>
#include <span>
#include <utility>

namespace uv::math::linear_algebra
{
template <std::floating_point T, typename F>
requires std::invocable<F&, std::size_t, std::size_t>
core::Matrix<T> generateIndexed(std::size_t rows, std::size_t cols, F&& f)
{
    auto&& func = std::forward<F>(f);

    core::Matrix<T> out(rows, cols);

    for (std::size_t i = 0; i < rows; ++i)
    {
        std::span<T> row{out[i]};

        for (std::size_t j = 0; j < cols; ++j)
        {
            row[j] = std::invoke(func, i, j);
        }
    }
    return out;
}

template <std::floating_point T, typename F>
requires std::invocable<F&, std::size_t, std::size_t, T>
core::Matrix<T> transformIndexed(const core::Matrix<T>& m, F&& f)
{
    auto&& func = std::forward<F>(f);

    return generateIndexed<T>(
        m.rows(),
        m.cols(),
        [&func, &m](std::size_t i, std::size_t j)
        {
            return std::invoke(func, i, j, m[i][j]);
        }
    );
}

template <std::floating_point T, typename F>
requires std::invocable<F&, std::size_t, std::size_t, T>
void transformIndexedInplace(core::Matrix<T>& m, F&& f)
{
    auto&& func = std::forward<F>(f);

    for (std::size_t i = 0; i < m.rows(); ++i)
    {
        std::span<T> row{m[i]};

        for (std::size_t j = 0; j < m.cols(); ++j)
        {
            row[j] = std::invoke(func, i, j, row[j]);
        }
    }
}

template <std::floating_point T>
core::Matrix<T> hadamard(const core::Matrix<T>& lhs, const core::Matrix<T>& rhs)
{
    core::Matrix<T> out{lhs};
    hadamardInplace(out, rhs);
    return out;
}

template <std::floating_point T>
void hadamardInplace(core::Matrix<T>& lhs, const core::Matrix<T>& rhs)
{
    const std::size_t numRows{lhs.rows()};
    const std::size_t numCols{lhs.cols()};

    REQUIRE_SAME_SIZE(numRows, rhs.rows());
    REQUIRE_SAME_SIZE(numCols, rhs.cols());

    for (std::size_t i = 0; i < numRows; ++i)
    {
        std::span<T> lhsRow{lhs[i]};
        std::span<const T> rhsRow{rhs[i]};

        for (std::size_t j = 0; j < numCols; ++j)
        {
            lhsRow[j] *= rhsRow[j];
        }
    }
}

template <std::floating_point T>
core::Matrix<T> hadamard(const core::Matrix<T>& lhs, const Vector<T>& rhs)
{
    core::Matrix<T> out{lhs};
    hadamardInplace(out, rhs);
    return out;
}

template <std::floating_point T>
void hadamardInplace(core::Matrix<T>& lhs, const Vector<T>& rhs)
{
    const std::size_t numRows{lhs.rows()};
    const std::size_t numCols{lhs.cols()};

    REQUIRE_SAME_SIZE(numRows, rhs.size());

    for (std::size_t i = 0; i < numRows; ++i)
    {
        std::span<T> row{lhs[i]};
        const T scale{rhs[i]};

        for (std::size_t j = 0; j < numCols; ++j)
        {
            row[j] *= scale;
        }
    }
}

template <std::floating_point T>
core::Matrix<T> divide(const core::Matrix<T>& lhs, const core::Matrix<T>& rhs)
{
    core::Matrix<T> out{lhs};
    divideInplace(out, rhs);
    return out;
}

template <std::floating_point T>
void divideInplace(core::Matrix<T>& lhs, const core::Matrix<T>& rhs)
{
    const std::size_t numRows{lhs.rows()};
    const std::size_t numCols{lhs.cols()};

    REQUIRE_SAME_SIZE(numRows, rhs.rows());
    REQUIRE_SAME_SIZE(numCols, rhs.cols());

    for (std::size_t i = 0; i < numRows; ++i)
    {
        std::span<T> lhsRow{lhs[i]};
        std::span<const T> rhsRow{rhs[i]};

        for (std::size_t j = 0; j < numCols; ++j)
        {
            lhsRow[j] /= rhsRow[j];
        }
    }
}

template <std::floating_point T> core::Matrix<T> reciprocal(const core::Matrix<T>& m)
{
    core::Matrix<T> out{m};
    reciprocalInplace(out);
    return out;
}

template <std::floating_point T> void reciprocalInplace(core::Matrix<T>& m) noexcept
{
    const std::size_t numCols{m.cols()};

    for (std::size_t i = 0; i < m.rows(); ++i)
    {
        std::span<T> mRow{m[i]};

        for (std::size_t j = 0; j < numCols; ++j)
        {
            mRow[j] = T{1} / mRow[j];
        }
    }
}

template <std::floating_point T> core::Matrix<T> square(const core::Matrix<T>& m)
{
    core::Matrix<T> out{m};
    squareInplace(out);
    return out;
}

template <std::floating_point T> void squareInplace(core::Matrix<T>& m) noexcept
{
    const std::size_t numCols{m.cols()};

    for (std::size_t i = 0; i < m.rows(); ++i)
    {
        std::span<T> row{m[i]};

        for (std::size_t j = 0; j < numCols; ++j)
        {
            row[j] *= row[j];
        }
    }
}

template <std::floating_point T> core::Matrix<T> sqrt(const core::Matrix<T>& m)
{
    core::Matrix<T> out{m};
    sqrtInplace(out);
    return out;
}

template <std::floating_point T> void sqrtInplace(core::Matrix<T>& m)
{
    const std::size_t numRows{m.rows()};

    for (std::size_t i = 0; i < numRows; ++i)
    {
        std::span<T> row{m[i]};

        REQUIRE_NON_NEGATIVE(row);

        squareRootInplace<T>(row, row);
    }
}
} // namespace uv::math::linear_algebra
