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

template <std::floating_point T>
T parseNumberCellOrThrow(
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

#include <IO/CSV/Detail/Read.inl>