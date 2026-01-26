// SPDX-License-Identifier: Apache-2.0
/*
 * File:        ConsoleRedirect.hpp
 * Author:      Alvaro Sanchez de Carlos
 * Created:     2025-12-08
 *
 * Description:
 *   [Brief description of what this file declares or implements.]
 *
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

#include "Utils/IO/Log.hpp"

#include <iostream>
#include <sstream>
#include <string>

namespace uv::utils
{
/**
 * @brief RAII helper to temporarily redirect std::cout to an internal buffer.
 *
 * Upon construction, the current std::cout stream buffer is replaced so that
 * all output written to std::cout is captured internally.
 *
 * Upon destruction, the original stream buffer is restored and the captured
 * output is forwarded to the UnifiedVol logging system.
 *
 * Typical use case:
 *   - Capturing verbose output from third-party libraries (e.g. Ceres)
 *     and redirecting it into the UnifiedVol logger.
 *
 * @note This class is intended for scoped use only.
 * @note Not thread-safe: redirects the global std::cout stream.
 */
struct ConsoleRedirect
{
    /**
     * @brief Redirect std::cout to an internal string buffer.
     *
     * Saves the original stream buffer so it can be restored on destruction.
     */
    explicit ConsoleRedirect()
    {
        oldBuf_ = std::cout.rdbuf(buffer_.rdbuf());
    }

    /**
     * @brief Restore std::cout and flush captured output to the logger.
     *
     * Restores the original std::cout stream buffer and logs the accumulated
     * output using the UnifiedVol logging facility.
     */
    ~ConsoleRedirect()
    {
        std::cout.rdbuf(oldBuf_);
        UV_INFO("\n" + buffer_.str());
    }

  private:
    std::stringstream buffer_; // Internal buffer capturing std::cout output
    std::streambuf* oldBuf_{}; // Original std::cout stream buffer
};
} // namespace uv::utils