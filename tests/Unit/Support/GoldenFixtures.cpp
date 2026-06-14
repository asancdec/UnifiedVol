// SPDX-License-Identifier: Apache-2.0

#include "Support/Golden.hpp"
#include "Support/TempFile.hpp"

#include <gtest/gtest.h>
#include <stdexcept>

TEST(GoldenFixtures, CommittedFixturesHaveValidSchemas)
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
        uv::tests::golden::readBlackKnownValue("tests/Golden/black_known_values.json")
    ));
}

TEST(GoldenFixtures, CommittedFixturesExposeExpectedReferenceFields)
{
    const auto black =
        uv::tests::golden::readBlackKnownValue("tests/Golden/black_known_values.json");
    EXPECT_DOUBLE_EQ(black.forward, 100.0);
    EXPECT_DOUBLE_EQ(black.strike, 100.0);

    const auto bspline =
        uv::tests::golden::readBSplineKnownValues("tests/Golden/bspline_known_values.json"
        );
    ASSERT_FALSE(bspline.x.empty());
    EXPECT_DOUBLE_EQ(bspline.x.front(), 0.0);

    const auto pipeline =
        uv::tests::golden::readExamplePipeline("tests/Golden/example_pipeline.json");
    ASSERT_FALSE(pipeline.sviParams.empty());
    EXPECT_DOUBLE_EQ(pipeline.sviParams.front().t, 0.083333332999999996);

    const auto synthetic = uv::tests::golden::readSyntheticSVICalibration(
        "tests/Golden/synthetic_svi_calibration.json"
    );
    ASSERT_FALSE(synthetic.calibratedParams.empty());
    EXPECT_DOUBLE_EQ(synthetic.calibratedParams.front().t, 0.5);
}

TEST(GoldenFixtures, RejectsMalformedBSplineKnownValues)
{
    const auto path = uv::tests::writeTempFile(
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

TEST(GoldenFixtures, RejectsInvalidBlackKnownValue)
{
    const auto path = uv::tests::writeTempFile(
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

TEST(GoldenFixtures, RejectsFractionalVolPointIndices)
{
    const auto path = uv::tests::writeTempFile(
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

TEST(GoldenFixtures, RejectsEmptySVIParams)
{
    const auto path = uv::tests::writeTempFile(
        "unifiedvol_empty_svi_params_golden.json",
        R"({
          "tolerance": 1e-8,
          "sviParams": [],
          "hestonParams": {
            "kappa": 2.0,
            "theta": 0.04,
            "sigma": 0.3,
            "rho": -0.5,
            "v0": 0.04
          },
          "marketVols": [{"maturity": 0, "strike": 0, "value": 0.2}],
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

TEST(GoldenFixtures, RejectsNonIncreasingSVIMaturities)
{
    const auto path = uv::tests::writeTempFile(
        "unifiedvol_non_increasing_svi_golden.json",
        R"({
          "tolerance": 1e-8,
          "sviParams": [
            {"t": 1.0, "a": 0.01, "b": 0.2, "rho": -0.2, "m": 0.0, "sigma": 0.3},
            {"t": 0.5, "a": 0.02, "b": 0.2, "rho": -0.2, "m": 0.0, "sigma": 0.3}
          ],
          "hestonParams": {
            "kappa": 2.0,
            "theta": 0.04,
            "sigma": 0.3,
            "rho": -0.5,
            "v0": 0.04
          },
          "marketVols": [{"maturity": 0, "strike": 0, "value": 0.2}],
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

TEST(GoldenFixtures, RejectsInvalidSVIRho)
{
    const auto path = uv::tests::writeTempFile(
        "unifiedvol_invalid_svi_rho_golden.json",
        R"({
          "tolerance": 1e-8,
          "sviParams": [
            {"t": 0.5, "a": 0.01, "b": 0.2, "rho": 1.0, "m": 0.0, "sigma": 0.3}
          ],
          "hestonParams": {
            "kappa": 2.0,
            "theta": 0.04,
            "sigma": 0.3,
            "rho": -0.5,
            "v0": 0.04
          },
          "marketVols": [{"maturity": 0, "strike": 0, "value": 0.2}],
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

TEST(GoldenFixtures, RejectsMissingErrorsObject)
{
    const auto path = uv::tests::writeTempFile(
        "unifiedvol_missing_errors_golden.json",
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
          "marketVols": [{"maturity": 0, "strike": 0, "value": 0.2}],
          "sviVols": [{"maturity": 0, "strike": 0, "value": 0.2}],
          "hestonVols": [{"maturity": 0, "strike": 0, "value": 0.2}]
        })"
    );

    EXPECT_THROW(
        static_cast<void>(uv::tests::golden::readExamplePipeline(path)),
        std::runtime_error
    );
}

TEST(GoldenFixtures, RejectsEmptyVolPointArrays)
{
    const auto path = uv::tests::writeTempFile(
        "unifiedvol_empty_vol_points_golden.json",
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
          "marketVols": [],
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
