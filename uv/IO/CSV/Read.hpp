// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <concepts>
#include <cstddef>
#include <istream>
#include <string>
#include <string_view>
#include <vector>

namespace uv::io::csv
{
template <class T> using StdVector = std::vector<T>;

struct Options
{
    bool allowPercent{true};
    bool allowExtraCols{true};
    bool skipBlankLines{true};
};

std::string_view trimView(std::string_view s) noexcept;

StdVector<std::string> splitComma(std::string_view line);

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
} // namespace uv::io::csv

#include "IO/CSV/Detail/Read.inl"