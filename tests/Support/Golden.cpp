// SPDX-License-Identifier: Apache-2.0

#include "Support/Golden.hpp"

#include "IO/JSON/Read.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

namespace uv::tests::golden
{
namespace
{
namespace json = uv::io::json;

void requireObject(const json::Value& value, const std::string& path)
{
    if (value.type != json::Value::Type::Object)
        throw std::runtime_error("Expected JSON object at: " + path);
}

void requireArray(const json::Value& value, const std::string& path)
{
    if (value.type != json::Value::Type::Array)
        throw std::runtime_error("Expected JSON array at: " + path);
}

void requireFinite(const double value, const std::string& path)
{
    if (!std::isfinite(value))
        throw std::runtime_error("Expected finite JSON number at: " + path);
}

void requirePositiveFinite(const double value, const std::string& path)
{
    requireFinite(value, path);
    if (value <= 0.0)
        throw std::runtime_error("Expected positive JSON number at: " + path);
}

std::size_t readIndex(const json::Value& tree, const std::string& path)
{
    const double value{tree.asNumber()};
    requireFinite(value, path);
    const double nearestInteger{std::round(value)};
    if (value < 0.0 || std::abs(value - nearestInteger) > 0.0)
        throw std::runtime_error("Expected non-negative integer JSON index at: " + path);
    return static_cast<std::size_t>(value);
}

models::svi::Params<double> readSVIParams(const json::Value& tree)
{
    requireObject(tree, "sviParams[]");
    return {
        tree.at("t").asNumber(),
        tree.at("a").asNumber(),
        tree.at("b").asNumber(),
        tree.at("rho").asNumber(),
        tree.at("m").asNumber(),
        tree.at("sigma").asNumber()
    };
}

Vector<models::svi::Params<double>>
readSVIParamsArray(const json::Value& tree, const std::string& path)
{
    requireArray(tree.at(path), path);
    Vector<models::svi::Params<double>> params;
    for (const auto& child : tree.at(path).array)
        params.emplace_back(readSVIParams(child));
    if (params.empty())
        throw std::runtime_error("Expected non-empty SVI params array at: " + path);
    return params;
}

Vector<VolPoint> readVolPoints(const json::Value& tree, const std::string& path)
{
    requireArray(tree.at(path), path);
    Vector<VolPoint> points;
    std::size_t i{};
    for (const auto& child : tree.at(path).array)
    {
        requireObject(child, path + "[]");
        points.push_back(
            {readIndex(
                 child.at("maturity"),
                 path + "[" + std::to_string(i) + "].maturity"
             ),
             readIndex(child.at("strike"), path + "[" + std::to_string(i) + "].strike"),
             child.at("value").asNumber()}
        );
        requireFinite(points.back().value, path + "[" + std::to_string(i) + "].value");
        ++i;
    }
    if (points.empty())
        throw std::runtime_error("Expected non-empty vol point array at: " + path);
    return points;
}

Vector<double> readDoubleArray(const json::Value& tree, const std::string& path)
{
    requireArray(tree.at(path), path);
    Vector<double> values;
    std::size_t i{};
    for (const auto& child : tree.at(path).array)
    {
        values.emplace_back(child.asNumber());
        requireFinite(values.back(), path + "[" + std::to_string(i) + "]");
        ++i;
    }
    if (values.empty())
        throw std::runtime_error("Expected non-empty numeric array at: " + path);
    return values;
}

void validateSVIParams(
    const Vector<models::svi::Params<double>>& params,
    const std::string& path
)
{
    double previousT{0.0};
    bool first{true};
    for (std::size_t i = 0; i < params.size(); ++i)
    {
        const auto& p = params[i];
        const std::string prefix{path + "[" + std::to_string(i) + "]"};
        requirePositiveFinite(p.t, prefix + ".t");
        requireFinite(p.a, prefix + ".a");
        requirePositiveFinite(p.b, prefix + ".b");
        requireFinite(p.rho, prefix + ".rho");
        requireFinite(p.m, prefix + ".m");
        requirePositiveFinite(p.sigma, prefix + ".sigma");
        if (p.rho <= -1.0 || p.rho >= 1.0)
            throw std::runtime_error("Expected rho in (-1, 1) at: " + prefix + ".rho");
        if (!first && p.t <= previousT)
            throw std::runtime_error("Expected increasing SVI maturities at: " + path);
        first = false;
        previousT = p.t;
    }
}

void validateHestonParams(const models::heston::Params<double>& params)
{
    requirePositiveFinite(params.kappa, "hestonParams.kappa");
    requirePositiveFinite(params.theta, "hestonParams.theta");
    requirePositiveFinite(params.sigma, "hestonParams.sigma");
    requireFinite(params.rho, "hestonParams.rho");
    requirePositiveFinite(params.v0, "hestonParams.v0");
    if (params.rho <= -1.0 || params.rho >= 1.0)
        throw std::runtime_error("Expected rho in (-1, 1) at: hestonParams.rho");
}
} // namespace

ExamplePipeline readExamplePipeline(const std::filesystem::path& path)
{
    const json::Value tree{json::read(path)};
    requireObject(tree, path.string());

    const auto& heston = tree.at("hestonParams");
    requireObject(heston, "hestonParams");
    ExamplePipeline fixture{
        .tolerance = tree.at("tolerance").asNumber(),
        .sviParams = readSVIParamsArray(tree, "sviParams"),
        .hestonParams =
            {heston.at("kappa").asNumber(),
             heston.at("theta").asNumber(),
             heston.at("sigma").asNumber(),
             heston.at("rho").asNumber(),
             heston.at("v0").asNumber()},
        .marketVols = readVolPoints(tree, "marketVols"),
        .sviVols = readVolPoints(tree, "sviVols"),
        .hestonVols = readVolPoints(tree, "hestonVols"),
        .meanAbsVolError = tree.at("errors").at("meanAbsVol").asNumber(),
        .maxAbsVolError = tree.at("errors").at("maxAbsVol").asNumber(),
        .rmseVolError = tree.at("errors").at("rmseVol").asNumber()
    };
    requirePositiveFinite(fixture.tolerance, "tolerance");
    validateSVIParams(fixture.sviParams, "sviParams");
    validateHestonParams(fixture.hestonParams);
    requireFinite(fixture.meanAbsVolError, "errors.meanAbsVol");
    requireFinite(fixture.maxAbsVolError, "errors.maxAbsVol");
    requireFinite(fixture.rmseVolError, "errors.rmseVol");
    return fixture;
}

SyntheticSVICalibration readSyntheticSVICalibration(const std::filesystem::path& path)
{
    const json::Value tree{json::read(path)};
    requireObject(tree, path.string());

    SyntheticSVICalibration fixture{
        .tolerance = tree.at("tolerance").asNumber(),
        .calibratedParams = readSVIParamsArray(tree, "calibratedParams"),
        .meanAbsTotalVarianceError =
            tree.at("errors").at("meanAbsTotalVariance").asNumber(),
        .maxAbsTotalVarianceError =
            tree.at("errors").at("maxAbsTotalVariance").asNumber(),
        .rmseTotalVarianceError = tree.at("errors").at("rmseTotalVariance").asNumber()
    };
    requirePositiveFinite(fixture.tolerance, "tolerance");
    validateSVIParams(fixture.calibratedParams, "calibratedParams");
    requireFinite(fixture.meanAbsTotalVarianceError, "errors.meanAbsTotalVariance");
    requireFinite(fixture.maxAbsTotalVarianceError, "errors.maxAbsTotalVariance");
    requireFinite(fixture.rmseTotalVarianceError, "errors.rmseTotalVariance");
    return fixture;
}

BSplineKnownValues readBSplineKnownValues(const std::filesystem::path& path)
{
    const json::Value tree{json::read(path)};
    requireObject(tree, path.string());

    BSplineKnownValues fixture{
        .tolerance = tree.at("tolerance").asNumber(),
        .controlPoints = readDoubleArray(tree, "controlPoints"),
        .knots = readDoubleArray(tree, "knots"),
        .x = readDoubleArray(tree, "x"),
        .y = readDoubleArray(tree, "y")
    };
    requirePositiveFinite(fixture.tolerance, "tolerance");
    if (fixture.x.size() != fixture.y.size())
        throw std::runtime_error("Expected B-spline x and y arrays to have the same size"
        );
    if (fixture.knots.size() <= fixture.controlPoints.size())
        throw std::runtime_error(
            "Expected B-spline knot count to exceed control point count"
        );
    for (std::size_t i = 1; i < fixture.knots.size(); ++i)
    {
        if (fixture.knots[i] < fixture.knots[i - 1])
            throw std::runtime_error("Expected non-decreasing B-spline knots");
    }
    return fixture;
}

BlackKnownValue readBlackKnownValue(const std::filesystem::path& path)
{
    const json::Value tree{json::read(path)};
    requireObject(tree, path.string());
    requireObject(tree.at("input"), "input");

    BlackKnownValue fixture{
        .tolerance = tree.at("tolerance").asNumber(),
        .t = tree.at("input").at("t").asNumber(),
        .discountFactor = tree.at("input").at("discountFactor").asNumber(),
        .forward = tree.at("input").at("forward").asNumber(),
        .volatility = tree.at("input").at("volatility").asNumber(),
        .strike = tree.at("input").at("strike").asNumber(),
        .price = tree.at("price").asNumber()
    };
    requirePositiveFinite(fixture.tolerance, "tolerance");
    requirePositiveFinite(fixture.t, "input.t");
    requirePositiveFinite(fixture.discountFactor, "input.discountFactor");
    requirePositiveFinite(fixture.forward, "input.forward");
    requirePositiveFinite(fixture.volatility, "input.volatility");
    requirePositiveFinite(fixture.strike, "input.strike");
    requireFinite(fixture.price, "price");
    return fixture;
}
} // namespace uv::tests::golden
