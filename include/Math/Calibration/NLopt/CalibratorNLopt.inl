/**
* CalibratorNLopt.inl
* Author: Alvaro Sanchez de Carlos
*/

#include "Errors/Errors.hpp"
#include "Utils/Log.hpp"
#include "Math/Calibration/Utils/CalibratorUtils.hpp"

#include <cstddef>
#include <format>
#include <cmath>
#include <chrono>
#include <algorithm> 

namespace uv
{
    template <::std::size_t N>
    CalibratorNLopt<N>::CalibratorNLopt(const CalibratorConfig<N>& config,
        ::nlopt::algorithm algo) :
        config_(config),
        opt_(algo, N),
        algo_(algo)
    {
        opt_.set_ftol_rel(config_.ftolRel);  // Stop when objective stops improving
        opt_.set_maxeval(config_.maxEval);   // Maximum number of evaluations
    }

    template <::std::size_t N>
    CalibratorNLopt<N> CalibratorNLopt<N>::fresh() const noexcept
    {
        return CalibratorNLopt{ config_, algo_ };
    }

    template <::std::size_t N>
    void CalibratorNLopt<N>::setGuessBounds(::std::array<double, N> initGuess,
        ::std::array<double, N> lowerBounds,
        ::std::array<double, N> upperBounds) noexcept
    {
        // Clamp initial guess within upper and lower bounds
        uv::clamp<N>(initGuess, lowerBounds, upperBounds, config_.paramNames);

        // Assign to class members
        initGuess_.assign(initGuess.begin(), initGuess.end());
        lowerBounds_.assign(lowerBounds.begin(), lowerBounds.end());
        upperBounds_.assign(upperBounds.begin(), upperBounds.end());

        // Configure NLopt with bounds and tolerances
        opt_.set_lower_bounds(lowerBounds_);
        opt_.set_upper_bounds(upperBounds_);
    }

    template <::std::size_t N>
    double CalibratorNLopt<N>::ObjectiveThunk(unsigned n, const double* x, double* grad, void* p) noexcept
    {
        auto* c = static_cast<ObjWrapCtx*>(p);
        if (c->iter) ++(*c->iter);
        return c->fn(n, x, grad, c->user);
    }

    template <::std::size_t N>
    void CalibratorNLopt<N>::addInequalityConstraint(
        NloptFunction c,
        void* data) noexcept
    {
        opt_.add_inequality_constraint(c, data, config_.tol);
    }

    template <::std::size_t N>
    void CalibratorNLopt<N>::setMinObjective(NloptFunction f, void* data) noexcept
    {
        iterCount_ = 0;
        objCtx_.fn = f;
        objCtx_.user = data;
        objCtx_.iter = &iterCount_;
        opt_.set_min_objective(&CalibratorNLopt<N>::ObjectiveThunk, &objCtx_);
    }

    template <::std::size_t N>
    ::std::vector<double> CalibratorNLopt<N>::optimize() noexcept
    {
        // Start time
        const auto t0 = ::std::chrono::high_resolution_clock::now();

        // Copy initial guess to working vector
        ::std::vector<double> x(initGuess_.begin(), initGuess_.end());
        double sse{ 0.0 };

        // Run optimizer
        ::nlopt::result code = opt_.optimize(x, sse);

        // End time
        const auto t1 = ::std::chrono::high_resolution_clock::now();

        // Warn if upper or lower bounds are touched
        auto near = [](double v, double bd) noexcept
            {
                constexpr double abs_eps{ 1e-12 };
                constexpr double rel_eps{ 1e-7 };
                return ::std::fabs(v - bd) <= (abs_eps + rel_eps * (::std::max)(::std::fabs(v), ::std::fabs(bd)));
            };

        for (::std::size_t i = 0; i < x.size(); ++i)
        {
            const double v{ x[i] };
            const double lb{ lowerBounds_[i] };
            const double ub{ upperBounds_[i] };

            UV_WARN(near(v, lb),
                ::std::format("Calibrator: parameter '{}' hit LOWER bound: v = {:.4f} (lb = {:.4f})",
                    config_.paramNames[i], v, lb));

            UV_WARN(near(v, ub),
                ::std::format("Calibrator: parameter '{}' hit UPPER bound: v = {:.4f} (ub = {:.4f})",
                    config_.paramNames[i], v, ub));
        }

        // Log results
        ::std::string paramsLine;
        paramsLine.reserve(N * 24);
        for (::std::size_t i = 0; i < N; ++i)
        {
            paramsLine += ::std::format("{}={:.4f}{}",
                config_.paramNames[i],
                x[i],
                (i + 1 < N ? "  " : ""));
        }

        UV_INFO(::std::format(
            "[Calib] {}  SSE={:.4e} ({:.2f} ms, {} it, {})",
            paramsLine,                                                                         // Parameter values
            sse,                                                                                // SSE                                               
            ::std::chrono::duration_cast<::std::chrono::microseconds>(t1 - t0).count() / 1000.0,    // Time in ms
            iterCount_,                                                                         // Number of iterations            
            (code > 0 ? "SUCCESS" : "FAIL")));                                                  // Success message

        return x;
    }

    template <::std::size_t N>
    const double& CalibratorNLopt<N>::eps() const noexcept
    {
        return config_.eps;
    }

    template <::std::size_t N>
    double CalibratorNLopt<N>::tol() const noexcept
    {
        return config_.tol;
    }
}
