// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Base/Macros/Unreachable.hpp"
#include "Optimization/NLopt/Algorithm.hpp"

#include <nlopt.hpp>

namespace uv::opt::nlopt::detail
{
inline constexpr ::nlopt::algorithm toNlopt(Algorithm a) noexcept
{
    switch (a)
    {
    case Algorithm::LD_MMA:
        return ::nlopt::LD_MMA;
    case Algorithm::LD_SLSQP:
        return ::nlopt::LD_SLSQP;
    case Algorithm::LN_BOBYQA:
        return ::nlopt::LN_BOBYQA;
    }

    UNREACHABLE(Algorithm, a);
}
} // namespace uv::opt::nlopt::detail