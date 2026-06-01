// SPDX-License-Identifier: Apache-2.0

#include "Base/Execution/ThreadPolicy.hpp"
#include "Base/Errors/Errors.hpp"

#include <algorithm>
#include <format>
#include <thread>

namespace uv::execution
{
int requestThreads(int numRequested)
{
    unsigned int hw{std::thread::hardware_concurrency()};
    int numAvailable = (hw == 0u) ? 1 : hw;

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