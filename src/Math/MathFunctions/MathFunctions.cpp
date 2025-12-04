/**
* MathFunctions.inl
* Author: Alvaro Sanchez de Carlos
*/

#include "Math/MathFunctions/MathFunctions.hpp"

#include <algorithm>
#include <format>

namespace uv
{
    double impliedVolBS(double mktPriceBS,
        double t,
        double r,
        double q,
        double S,
        double K,
        bool isCall,
        double ftolAbs,
        unsigned int maxEval) noexcept
    {
        // Convert all inputs to long double
        const long double mktPriceBS_{ static_cast<long double>(mktPriceBS) };
        const long double t_{ static_cast<long double>(t) };
        const long double r_{ static_cast<long double>(r) };
        const long double q_{ static_cast<long double>(q) };
        const long double S_{ static_cast<long double>(S) };
        const long double K_{ static_cast<long double>(K) };
        const long double ftolAbs_{ static_cast<long double>(ftolAbs) };

        // Set initial guess using forward-moneyness heuristic
        const long double m{ std::log((S_ * std::exp((r_ - q_) * t_)) / K_) };
        long double vol
        {
            std::clamp
            (
                (std::fabs(m) < 1e-6L) ? 0.3L :
                std::sqrt(2.0L * std::fabs(m) / t_),
                1e-3L,
                5.0L
            )
        };

        // Solve for implied volatility using Halley's method
        for (unsigned int i = 0U; i < maxEval; ++i)
        {
            // Calculate Black-Scholes price
            const long double priceBS{ blackScholes(t_, r_, q_, vol, S_, K_, isCall) };

            // Evaluate objective function
            const long double objEval{ priceBS - mktPriceBS_ };

            // If target tolerance is achieved, exit optimization
            if (std::fabs(objEval) < ftolAbs_) break;

            // Calculate d1
            const long double d1{ d1BS(t_, r_, q_, vol, S_, K_) };

            // Calculate first derivative
            const long double vega{ vegaBS(d1, t_, q_, S_) };

            // Calculate second derivative
            const long double volga{ volgaBS(vega, d1, t_, vol) };

            // Update volatility
            vol -= (2.0L * objEval * vega) / (2.0L * (vega * vega) - objEval * volga);

            // Volatility safequard
            vol = std::clamp(vol, std::numeric_limits<long double>::epsilon(), 100.0L);
        }

        // Calculate residual
        const long double objEval{ blackScholes(t_, r_, q_, vol, S_, K_, isCall) - mktPriceBS_ };
        UV_WARN(std::fabs(objEval) > ftolAbs_,
            std::format("impliedVolBS: no convergence after {} iterations "
                "(|f| = {:.3e} > tol = {:.3e}, vol = {:.6Lf})",
                maxEval, std::fabs(objEval), ftolAbs_, vol));

        return vol;
    }
}
