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

#include "Base/Errors/Errors.hpp"

#include <cctype>
#include <charconv>
#include <system_error>

namespace uv::io::csv
{

template <std::floating_point T>
T parseNumberCellOrThrow(
    std::string_view raw,
    std::string_view what,
    std::size_t lineNo,
    std::size_t colNo,
    Options opt
)
{
    auto s = trimView(raw);
    if (s.empty())
    {
        errors::raise(
            errors::ErrorCode::DataFormat,
            std::string(what) + " is empty at line " + std::to_string(lineNo) + ", col " +
                std::to_string(colNo)
        );
    }

    bool percent = false;
    if (opt.allowPercent && !s.empty() && s.back() == '%')
    {
        percent = true;
        s.remove_suffix(1);
        s = trimView(s);
        if (s.empty())
        {
            errors::raise(
                errors::ErrorCode::DataFormat,
                "Lonely % at line " + std::to_string(lineNo) + ", col " +
                    std::to_string(colNo)
            );
        }
    }

    double tmp = 0.0;
    const char* begin = s.data();
    const char* end = s.data() + s.size();

    auto [ptr, ec] = std::from_chars(begin, end, tmp, std::chars_format::general);

    if (ec != std::errc{} || ptr != end)
    {
        errors::raise(
            errors::ErrorCode::DataFormat,
            std::string("Non-numeric ") + std::string(what) + " \"" + std::string(raw) +
                "\" at line " + std::to_string(lineNo) + ", col " + std::to_string(colNo)
        );
    }

    T out = static_cast<T>(tmp);
    if (percent)
        out *= static_cast<T>(0.01);

    return out;
}

template <std::floating_point T, template <class> class Vector>
LabeledDense<T, Vector>
readLabeledDenseOrThrow(std::istream& is, std::string_view filenameForErrors, Options opt)
{
    std::string line;

    if (!std::getline(is, line))
    {
        errors::raise(
            errors::ErrorCode::DataFormat,
            "CSV file is empty: " + std::string(filenameForErrors)
        );
    }

    const auto headerCells = splitComma(line);
    if (headerCells.size() < 2)
    {
        errors::raise(
            errors::ErrorCode::DataFormat,
            "Header must have at least 2 columns (label + >=1 numeric col): " +
                std::string(filenameForErrors)
        );
    }

    LabeledDense<T, Vector> out;
    out.cols = headerCells.size() - 1;

    out.colLabels.reserve(out.cols);
    for (std::size_t j = 1; j < headerCells.size(); ++j)
    {
        out.colLabels.push_back(
            parseNumberCellOrThrow<T>(headerCells[j], "header value", 1, j + 1, opt)
        );
    }

    std::size_t lineNo = 1;
    while (std::getline(is, line))
    {
        ++lineNo;

        if (opt.skipBlankLines)
        {
            if (trimView(line).empty())
                continue;
        }

        const auto cells = splitComma(line);
        if (cells.size() < 2)
        {
            errors::raise(
                errors::ErrorCode::DataFormat,
                "Row has fewer than 2 columns at line " + std::to_string(lineNo) +
                    " in " + std::string(filenameForErrors)
            );
        }

        if (cells.size() < 1 + out.cols)
        {
            errors::raise(
                errors::ErrorCode::DataFormat,
                "Row " + std::to_string(lineNo) + " has only " +
                    std::to_string(cells.size() - 1) + " data cols; expected " +
                    std::to_string(out.cols)
            );
        }

        if (!opt.allowExtraCols && cells.size() != 1 + out.cols)
        {
            errors::raise(
                errors::ErrorCode::DataFormat,
                "Row " + std::to_string(lineNo) + " has extra columns (got " +
                    std::to_string(cells.size() - 1) + ", expected " +
                    std::to_string(out.cols) + ")"
            );
        }

        const T rowLabel =
            parseNumberCellOrThrow<T>(cells[0], "row label", lineNo, 1, opt);
        out.rowLabels.push_back(rowLabel);

        out.values.reserve(out.values.size() + out.cols);
        for (std::size_t j = 0; j < out.cols; ++j)
        {
            out.values.push_back(
                parseNumberCellOrThrow<T>(cells[1 + j], "cell", lineNo, (1 + j) + 1, opt)
            );
        }

        ++out.rows;
    }

    if (out.rows == 0)
    {
        errors::raise(
            errors::ErrorCode::DataFormat,
            "CSV file has no data rows: " + std::string(filenameForErrors)
        );
    }

    if (out.values.size() != out.rows * out.cols)
    {
        errors::raise(
            errors::ErrorCode::DataFormat,
            "Internal error: values size mismatch"
        );
    }

    return out;
}
} // namespace uv::io::csv