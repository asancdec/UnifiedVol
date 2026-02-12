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