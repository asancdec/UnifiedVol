// SPDX-License-Identifier: Apache-2.0

#include "Golden.hpp"

#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <stdexcept>
#include <string>

namespace
{
std::filesystem::path writeFixture(const std::string& name, const std::string& contents)
{
    const auto path = std::filesystem::temp_directory_path() / name;
    std::ofstream out{path};
    out << contents;
    return path;
}
} // namespace

TEST(RegressionGoldenFixtures, CommittedFixturesHaveValidSchemas)
{
    EXPECT_NO_THROW(static_cast<void>(
        uv::tests::golden::readExamplePipeline("tests/Golden/example_pipeline.json")
    ));
    EXPECT_NO_THROW(static_cast<void>(uv::tests::golden::readSyntheticSVICalibration(
        "tests/Golden/synthetic_svi_calibration.json"
    )));
    EXPECT_NO_THROW(static_cast<void>(uv::tests::golden::readBSplineKnownValues(
        "tests/Golden/bspline_known_values.json"
    )));
    EXPECT_NO_THROW(static_cast<void>(
        uv::tests::golden::readBlackKnownValue("tests/Golden/black_known_value.json")
    ));
}

TEST(RegressionGoldenFixtures, RejectsMalformedBSplineKnownValues)
{
    const auto path = writeFixture(
        "unifiedvol_bad_bspline_golden.json",
        R"({
          "tolerance": 1e-12,
          "controlPoints": [0.0, 1.0, 2.0],
          "knots": [0.0, 0.0, 1.0, 1.0],
          "x": [0.0, 1.0],
          "y": [0.0]
        })"
    );

    EXPECT_THROW(
        static_cast<void>(uv::tests::golden::readBSplineKnownValues(path)),
        std::runtime_error
    );
}

TEST(RegressionGoldenFixtures, RejectsInvalidBlackKnownValue)
{
    const auto path = writeFixture(
        "unifiedvol_bad_black_golden.json",
        R"({
          "tolerance": 1e-12,
          "input": {
            "t": 1.0,
            "discountFactor": 0.95,
            "forward": 100.0,
            "volatility": -0.20,
            "strike": 100.0
          },
          "price": 7.5
        })"
    );

    EXPECT_THROW(
        static_cast<void>(uv::tests::golden::readBlackKnownValue(path)),
        std::runtime_error
    );
}

TEST(RegressionGoldenFixtures, RejectsFractionalVolPointIndices)
{
    const auto path = writeFixture(
        "unifiedvol_bad_pipeline_golden.json",
        R"({
          "tolerance": 1e-8,
          "sviParams": [
            {"t": 0.5, "a": 0.01, "b": 0.2, "rho": -0.2, "m": 0.0, "sigma": 0.3}
          ],
          "hestonParams": {
            "kappa": 2.0,
            "theta": 0.04,
            "sigma": 0.3,
            "rho": -0.5,
            "v0": 0.04
          },
          "marketVols": [{"maturity": 0.5, "strike": 0, "value": 0.2}],
          "sviVols": [{"maturity": 0, "strike": 0, "value": 0.2}],
          "hestonVols": [{"maturity": 0, "strike": 0, "value": 0.2}],
          "errors": {
            "meanAbsVol": 0.0,
            "maxAbsVol": 0.0,
            "rmseVol": 0.0
          }
        })"
    );

    EXPECT_THROW(
        static_cast<void>(uv::tests::golden::readExamplePipeline(path)),
        std::runtime_error
    );
}
