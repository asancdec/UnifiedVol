// SPDX-License-Identifier: Apache-2.0

#include "IO/Console/Report.hpp"

#include <algorithm>
#include <limits>

namespace uv::io::report::detail
{

int precision(unsigned int value) noexcept
{
    const unsigned int capped{std::min<unsigned int>(
        value,
        static_cast<unsigned int>(std::numeric_limits<int>::max())
    )};
    return static_cast<int>(capped);
}

} // namespace uv::io::report::detail
