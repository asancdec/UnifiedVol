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
    using enum NLoptStatus;

    switch (s)
    {
    case Success:
        return "SUCCESS";
    case StopvalReached:
        return "STOPVAL_REACHED";
    case FtolReached:
        return "FTOL_REACHED";
    case XtolReached:
        return "XTOL_REACHED";
    case MaxevalReached:
        return "MAXEVAL_REACHED";
    case MaxtimeReached:
        return "MAXTIME_REACHED";

    case Failure:
        return "FAILURE";
    case InvalidArgs:
        return "INVALID_ARGS";
    case OutOfMemory:
        return "OUT_OF_MEMORY";
    case RoundoffLimited:
        return "ROUNDOFF_LIMITED";
    case ForcedStop:
        return "FORCED_STOP";
    }

    UNREACHABLE(NLoptStatus, s);
}

std::string_view toString(const ::nlopt::result& r) noexcept
{
    return toString(toStatus(r));
}

} // namespace uv::opt::nlopt::detail
