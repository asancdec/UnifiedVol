// SPDX-License-Identifier: Apache-2.0

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

NLoptStatus toStatus(const ::nlopt::result&) noexcept;

std::string_view toString(NLoptStatus) noexcept;

std::string_view toString(const ::nlopt::result&) noexcept;

} // namespace uv::opt::nlopt::detail
