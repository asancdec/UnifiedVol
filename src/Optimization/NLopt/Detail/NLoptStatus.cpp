// SPDX-License-Identifier: Apache-2.0

#include "Optimization/NLopt/Detail/NLoptStatus.hpp"
#include "Base/Macros/Unreachable.hpp"

namespace uv::opt::nlopt::detail
{
NLoptStatus toStatus(const ::nlopt::result& r) noexcept
{
    return static_cast<NLoptStatus>(static_cast<int>(r));
}

std::string_view toString(NLoptStatus s) noexcept
{
    switch (s)
    {
    case NLoptStatus::Success:
        return "SUCCESS";
    case NLoptStatus::StopvalReached:
        return "STOPVAL_REACHED";
    case NLoptStatus::FtolReached:
        return "FTOL_REACHED";
    case NLoptStatus::XtolReached:
        return "XTOL_REACHED";
    case NLoptStatus::MaxevalReached:
        return "MAXEVAL_REACHED";
    case NLoptStatus::MaxtimeReached:
        return "MAXTIME_REACHED";

    case NLoptStatus::Failure:
        return "FAILURE";
    case NLoptStatus::InvalidArgs:
        return "INVALID_ARGS";
    case NLoptStatus::OutOfMemory:
        return "OUT_OF_MEMORY";
    case NLoptStatus::RoundoffLimited:
        return "ROUNDOFF_LIMITED";
    case NLoptStatus::ForcedStop:
        return "FORCED_STOP";
    }

    UNREACHABLE(NLoptStatus, s);
}

std::string_view toString(const ::nlopt::result& r) noexcept
{
    return toString(toStatus(r));
}

} // namespace uv::opt::nlopt::detail
