// SPDX-License-Identifier: Apache-2.0

#include "IO/CSV/Detail/Read.hpp"

#include <cctype>
#include <sstream>

namespace uv::io::csv::detail
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

    while (start <= line.size())
    {
        const std::size_t pos = line.find(',', start);
        if (pos == std::string_view::npos)
        {
            out.emplace_back(line.substr(start));
            return out;
        }
        out.emplace_back(line.substr(start, pos - start));
        start = pos + 1;
    }

    return out;
}

} // namespace uv::io::csv::detail
