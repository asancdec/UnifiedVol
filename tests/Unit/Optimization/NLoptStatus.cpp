// SPDX-License-Identifier: Apache-2.0

#include "Optimization/NLopt/Detail/NLoptStatus.hpp"

#include <gtest/gtest.h>
#include <nlopt.hpp>

namespace nlopt_detail = uv::opt::nlopt::detail;

TEST(UnitOptimizationNLoptStatus, MapsResultCodesToStatuses)
{
    EXPECT_EQ(
        nlopt_detail::toStatus(::nlopt::SUCCESS),
        nlopt_detail::NLoptStatus::Success
    );
    EXPECT_EQ(
        nlopt_detail::toStatus(::nlopt::STOPVAL_REACHED),
        nlopt_detail::NLoptStatus::StopvalReached
    );
    EXPECT_EQ(
        nlopt_detail::toStatus(::nlopt::FTOL_REACHED),
        nlopt_detail::NLoptStatus::FtolReached
    );
    EXPECT_EQ(
        nlopt_detail::toStatus(::nlopt::XTOL_REACHED),
        nlopt_detail::NLoptStatus::XtolReached
    );
    EXPECT_EQ(
        nlopt_detail::toStatus(::nlopt::MAXEVAL_REACHED),
        nlopt_detail::NLoptStatus::MaxevalReached
    );
    EXPECT_EQ(
        nlopt_detail::toStatus(::nlopt::MAXTIME_REACHED),
        nlopt_detail::NLoptStatus::MaxtimeReached
    );
    EXPECT_EQ(
        nlopt_detail::toStatus(::nlopt::FAILURE),
        nlopt_detail::NLoptStatus::Failure
    );
    EXPECT_EQ(
        nlopt_detail::toStatus(::nlopt::INVALID_ARGS),
        nlopt_detail::NLoptStatus::InvalidArgs
    );
    EXPECT_EQ(
        nlopt_detail::toStatus(::nlopt::OUT_OF_MEMORY),
        nlopt_detail::NLoptStatus::OutOfMemory
    );
    EXPECT_EQ(
        nlopt_detail::toStatus(::nlopt::ROUNDOFF_LIMITED),
        nlopt_detail::NLoptStatus::RoundoffLimited
    );
    EXPECT_EQ(
        nlopt_detail::toStatus(::nlopt::FORCED_STOP),
        nlopt_detail::NLoptStatus::ForcedStop
    );
}

TEST(UnitOptimizationNLoptStatus, ConvertsStatusesToStrings)
{
    EXPECT_EQ(nlopt_detail::toString(nlopt_detail::NLoptStatus::Success), "SUCCESS");
    EXPECT_EQ(
        nlopt_detail::toString(nlopt_detail::NLoptStatus::StopvalReached),
        "STOPVAL_REACHED"
    );
    EXPECT_EQ(
        nlopt_detail::toString(nlopt_detail::NLoptStatus::FtolReached),
        "FTOL_REACHED"
    );
    EXPECT_EQ(
        nlopt_detail::toString(nlopt_detail::NLoptStatus::XtolReached),
        "XTOL_REACHED"
    );
    EXPECT_EQ(
        nlopt_detail::toString(nlopt_detail::NLoptStatus::MaxevalReached),
        "MAXEVAL_REACHED"
    );
    EXPECT_EQ(
        nlopt_detail::toString(nlopt_detail::NLoptStatus::MaxtimeReached),
        "MAXTIME_REACHED"
    );
    EXPECT_EQ(nlopt_detail::toString(nlopt_detail::NLoptStatus::Failure), "FAILURE");
    EXPECT_EQ(
        nlopt_detail::toString(nlopt_detail::NLoptStatus::InvalidArgs),
        "INVALID_ARGS"
    );
    EXPECT_EQ(
        nlopt_detail::toString(nlopt_detail::NLoptStatus::OutOfMemory),
        "OUT_OF_MEMORY"
    );
    EXPECT_EQ(
        nlopt_detail::toString(nlopt_detail::NLoptStatus::RoundoffLimited),
        "ROUNDOFF_LIMITED"
    );
    EXPECT_EQ(
        nlopt_detail::toString(nlopt_detail::NLoptStatus::ForcedStop),
        "FORCED_STOP"
    );
}

TEST(UnitOptimizationNLoptStatus, ConvertsRawResultsToStrings)
{
    EXPECT_EQ(nlopt_detail::toString(::nlopt::SUCCESS), "SUCCESS");
    EXPECT_EQ(nlopt_detail::toString(::nlopt::ROUNDOFF_LIMITED), "ROUNDOFF_LIMITED");
}
