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

#include "Base/Errors/Validate.hpp"
#include "Base/Errors/Errors.hpp"

#include <format>

namespace uv::errors
{

void validateSameSize(
    std::size_t a,
    std::size_t b,
    std::string_view what,
    std::source_location loc
)
{
    if (a != b) [[unlikely]]
    {
        raise(
            ErrorCode::InvalidArgument,
            std::format("{} size mismatch: {} != {}", what, a, b),
            loc
        );
    }
}

void validateState(bool ok, std::string_view message, std::source_location loc)
{
    if (!ok) [[unlikely]]
    {
        raise(ErrorCode::InvalidState, message, loc);
    }
}

void validateDirCreated(
    bool ok,
    const std::filesystem::path& dir,
    std::source_location loc
)
{
    if (!ok) [[unlikely]]
    {
        raise(
            ErrorCode::FileIO,
            std::format("Failed to create directory: {}", dir.string()),
            loc
        );
    }
}

void validateFileOpened(
    bool ok,
    const std::filesystem::path& file,
    std::source_location loc
)
{
    if (!ok) [[unlikely]]
    {
        raise(
            ErrorCode::FileIO,
            std::format("Unable to open file: {}", file.string()),
            loc
        );
    }
}
} // namespace uv::errors