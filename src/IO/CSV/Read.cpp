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

#include <IO/CSV/Read.hpp>

#include <cctype>
#include <sstream>

namespace uv::io::csv
{

std::string_view trimView(std::string_view s) noexcept
{
    auto issp = [](unsigned char c)
    {
        return std::isspace(c) != 0;
    };

    while (!s.empty() && issp(static_cast<unsigned char>(s.front())))
        s.remove_prefix(1);
    while (!s.empty() && issp(static_cast<unsigned char>(s.back())))
        s.remove_suffix(1);
    return s;
}

StdVector<std::string> splitComma(std::string_view line)
{
    StdVector<std::string> out;
    std::size_t start = 0;

    while (true)
    {
        const std::size_t pos = line.find(',', start);
        if (pos == std::string_view::npos)
        {
            out.emplace_back(line.substr(start));
            break;
        }
        out.emplace_back(line.substr(start, pos - start));
        start = pos + 1;

        if (start == line.size())
        {
            out.emplace_back("");
            break;
        }
    }
    return out;
}

} // namespace uv::io::csv