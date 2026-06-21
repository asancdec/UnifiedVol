// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Base/Types.hpp"
#include "Core/Matrix.hpp"

#include <concepts>
#include <cstddef>
#include <istream>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

namespace uv::io::csv::detail
{
template <class T> using StdVector = std::vector<T>;

struct Options
{
    bool allowPercent{true};
    bool allowExtraCols{true};
    bool skipBlankLines{true};
};

std::string_view trimView(std::string_view) noexcept;

StdVector<std::string> splitComma(std::string_view);

template <std::floating_point T> T parseNumberCellOrThrow(
    std::string_view raw,
    std::string_view what,
    std::size_t lineNo,
    std::size_t colNo,
    Options opt = {}
);

template <std::floating_point T, template <class> class Vector = StdVector>
struct LabeledDense
{
    Vector<T> rowLabels;
    Vector<T> colLabels;
    Vector<T> values;
    std::size_t rows{0};
    std::size_t cols{0};
};

template <std::floating_point T, template <class> class Vector = StdVector>
LabeledDense<T, Vector> readLabeledDenseOrThrow(
    std::istream& is,
    std::string_view filenameForErrors,
    Options opt = {}
);

template <std::floating_point T> std::tuple<Vector<T>, Vector<T>, core::Matrix<T>>
readLabeledMatrixCsv(const std::string& filename, const Options& opt = {});
} // namespace uv::io::csv::detail

#include "IO/CSV/Detail/Read.inl"
