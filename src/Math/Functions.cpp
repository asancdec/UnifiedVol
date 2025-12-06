/**
* Functions.inl
* Author: Alvaro Sanchez de Carlos
*/

#include "Math/Functions.hpp"
#include "Utils/IO/Log.hpp"
#include "Utils/Aux/Errors.hpp"

#include <algorithm>
#include <format>

namespace uv::math
{
    Real impliedVolBS(Real mktPriceBS,
        Real T,
        Real r,
        Real q,
        Real S,
        Real K,
        bool isCall)
    {

        // ---------- Optimization parameters ----------

        const Real TOL{ 1e-14 };
        const unsigned int EVAL{ 100 };

        // ---------- Sanity checks ----------

        // Spot must be positive
        UV_REQUIRE(S > Real(0.0),
            ErrorCode::InvalidArgument,
            std::format("impliedVolBS: spot S = {:.6f} must be > 0", S)
        );

        // Strike must be positive
        UV_REQUIRE(K > Real(0.0),
            ErrorCode::InvalidArgument,
            std::format("impliedVolBS: strike K = {:.6f} must be > 0", K)
        );

        // Tenor must be positive
        UV_REQUIRE(T > Real(0.0),
            ErrorCode::InvalidArgument,
            std::format("impliedVolBS: tenor T = {:.6f} must be > 0", T)
        );

        // Market price must be non-negative
        UV_REQUIRE(mktPriceBS >= Real(0.0),
            ErrorCode::InvalidArgument,
            std::format("impliedVolBS: market price = {:.6f} must be >= 0", mktPriceBS)
        );

        // ---------- Initial guess ----------

        // Log(F/K)
        const Real logFM{ std::log((S * std::exp((r - q) * T)) / K) };

        // Heuristic guess
        const Real volGuess
        {
            (std::fabs(logFM) < Real(1e-6))
                ? Real(0.3)
                : std::sqrt(Real(2.0) * std::fabs(logFM) / T)
        };

        // Clamp bounds
        const Real volMin{ Real(1e-4) };
        const Real volMax{ Real(5.0) };

        // Clamp
        Real vol{ std::clamp(volGuess, volMin, volMax) };

        // Warn if clamped
        UV_WARN(vol != volGuess,
            std::format(
                "impliedVolBS: initial guess = {:.6f} clamped to {:.6f} (lb = {:.4f}, ub = {:.4f})",
                volGuess, vol, volMin, volMax
            )
        );

        // ---------- Halley's method ----------

        for (unsigned int i = 0U; i < EVAL; ++i)
        {

            const Real priceBS{ blackScholes(T, r, q, vol, S, K, isCall) };
            const Real objEval{ priceBS - mktPriceBS };

            // Check Absolute and relative tolerance threshold
            if (std::fabs(objEval) < TOL * (Real(1.0) + priceBS)) break;

            // Derivatives
            const Real d1{ d1BS(T, r, q, vol, S, K) };
            const Real vega{ vegaBS(d1, T, q, S) };
            const Real volga{ volgaBS(vega, d1, T, vol) };

            // Update volatility using Halley's method
            vol -= (Real(2.0) * objEval * vega) / (Real(2.0) * (vega * vega) - objEval * volga);
        }

        // ---------- Evaluate calibration ----------

        // Throw if volatility is negative or larger than 100
        UV_REQUIRE(
            (vol > Real(0.0)) && (vol < Real(100.0)),
            ErrorCode::CalibrationError,
            std::format(
                "impliedVolBS: resulting volatility out of bounds: vol = {:.6f} "
                "(lb = {:.4f}, ub = {:.4f})",
                vol, Real(0.0), Real(100.0)
            )
        );

        const Real finalPrice{ blackScholes(T, r, q, vol, S, K, isCall) };
        const Real finalResidual{ finalPrice - mktPriceBS };
        const Real finalTol{ TOL * (Real(1.0) + std::fabs(finalPrice)) };

        // Warn if no convergence
        UV_WARN(std::fabs(finalResidual) > finalTol,
            std::format(
                "impliedVolBS: no convergence after {} iterations "
                "(|f| = {:.3e} > tol = {:.3e}, vol = {:.6f})",
                EVAL, std::fabs(finalResidual), finalTol, vol
            )
        );

        return vol;
    }
} // namespace uv::math
