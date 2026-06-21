// SPDX-License-Identifier: Apache-2.0

#include "Base/Execution/ThreadPolicy.hpp"
#include "Base/Errors/Errors.hpp"

#include <algorithm>
#include <format>
#include <limits>
#include <thread>

namespace uv::execution
{
namespace detail
{
int availableThreads() noexcept
{
    const unsigned int hw{std::thread::hardware_concurrency()};
    if (hw == 0u)
    {
        return 1;
    }

    const unsigned int capped{std::min<unsigned int>(
        hw,
        static_cast<unsigned int>(std::numeric_limits<int>::max())
    )};
    return static_cast<int>(capped);
}
} // namespace detail

int requestThreads(int numRequested)
{
    const int numAvailable{detail::availableThreads()};

    if (numRequested < 0)
    {
        return std::clamp(numAvailable + numRequested + 1, 1, numAvailable);
    }
    if (numRequested > 0)
    {

        return std::clamp(numRequested, 1, numAvailable);
    }

    errors::raise(
        errors::ErrorCode::InvalidArgument,
        std::format("Cannot request {} number of threads", numRequested)
    );
}

} // namespace uv::execution
