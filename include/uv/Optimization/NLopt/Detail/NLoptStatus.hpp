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

#include <nlopt.hpp>
#include <string_view>

namespace uv::opt::nlopt::detail
{

enum class NLoptStatus : int
{
    Success = 1,
    StopvalReached = 2,
    FtolReached = 3,
    XtolReached = 4,
    MaxevalReached = 5,
    MaxtimeReached = 6,

    Failure = -1,
    InvalidArgs = -2,
    OutOfMemory = -3,
    RoundoffLimited = -4,
    ForcedStop = -5
};

NLoptStatus toStatus(const ::nlopt::result& r) noexcept;

std::string_view toString(NLoptStatus s) noexcept;

std::string_view toString(const ::nlopt::result& r) noexcept;

} // namespace uv::opt::nlopt::detail