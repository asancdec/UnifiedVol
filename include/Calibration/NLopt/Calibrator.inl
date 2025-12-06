/**
* Calibrator.inl
* Author: Alvaro Sanchez de Carlos
*/

#include "Utils/Aux/Helpers.hpp"

namespace uv::cal::nlopt
{
    template <std::size_t N, ::nlopt::algorithm Algo>
    Calibrator<N, Algo>::Calibrator(const Config<N>& config) :
        config_(config), opt_(Algo, N), timer_(),
        lowerBounds_(), upperBounds_(), initGuess_(),
        userFn_(nullptr), userData_(nullptr), iterCount_(0U) 
    {
        opt_.set_ftol_rel(config_.ftolRel); 
        opt_.set_maxeval(config_.maxEval);
    }

    template <std::size_t N, ::nlopt::algorithm Algo>
    Calibrator<N, Algo> Calibrator<N, Algo>::fresh() const noexcept
    {
        return Calibrator<N, Algo>{ config_ };
    }

    template <std::size_t N, ::nlopt::algorithm Algo>
    void Calibrator<N, Algo>::setGuessBounds(std::array<double, N> initGuess,
        std::array<double, N> lowerBounds,
        std::array<double, N> upperBounds) noexcept
    {
        // Clamp initial guess within upper and lower bounds
        utils::clamp<N>(initGuess, lowerBounds, upperBounds, config_.paramNames);

        // Store arrays
        initGuess_ = initGuess;
        lowerBounds_ = lowerBounds;
        upperBounds_ = upperBounds;

        // Configure NLopt with bounds and tolerances
        opt_.set_lower_bounds(std::vector<double>(lowerBounds_.begin(), lowerBounds_.end()));
        opt_.set_upper_bounds(std::vector<double>(upperBounds_.begin(), upperBounds_.end()));
    }

    template <std::size_t N, ::nlopt::algorithm Algo>
    double Calibrator<N, Algo>::ObjectiveThunk(unsigned n, const double* x, double* grad, void* p) noexcept
    {
        auto* self = static_cast<Calibrator<N, Algo>*>(p);
        ++self->iterCount_;
        return self->userFn_ ? self->userFn_(n, x, grad, self->userData_) : 0.0;
    }

    template <std::size_t N, ::nlopt::algorithm Algo>
    void Calibrator<N, Algo>::addInequalityConstraint(
        NloptFunction c,
        void* data) noexcept
    {
        opt_.add_inequality_constraint(c, data, config_.tol);
    }

    template <std::size_t N, ::nlopt::algorithm Algo>
    void Calibrator<N, Algo>::setMinObjective(NloptFunction f, void* data) noexcept
    {
        iterCount_ = 0U;
        userFn_ = f;
        userData_ = data;

        // Route NLopt callback into this instance
        opt_.set_min_objective(&Calibrator<N, Algo>::ObjectiveThunk, this);
    }

    template <std::size_t N, ::nlopt::algorithm Algo>
    std::vector<double> Calibrator<N, Algo>::optimize() noexcept
    {
        // Copy initial guess to working vector
        std::vector<double> x(initGuess_.cbegin(), initGuess_.cend());
        double sse{ 0.0 };

        // Start timer
        timer_.StartStopWatch();

        // Run optimizer
        ::nlopt::result successCode = opt_.optimize(x, sse);

        // End timer
        timer_.StopStopWatch();

        // Warn if upper or lower bounds are touched
        utils::warnBoundsHit
        (
            x,
            lowerBounds_,
            upperBounds_,
            config_.paramNames
        );

        // Log calibration results 
        utils::logResults(
            x,                                       // Parameters
            config_.paramNames,                      // Parameter names
            sse,                                     // SSE
            iterCount_,                              // Iterations
            timer_.GetTime<std::milli>(),            // Elapsed [ms]
            (successCode > ::nlopt::FAILURE)           // Success flag
        );

        return x;
    }

    template <std::size_t N, ::nlopt::algorithm Algo>
    const double& Calibrator<N, Algo>::eps() const noexcept
    {
        return config_.eps;
    }

    template <std::size_t N, ::nlopt::algorithm Algo>
    double Calibrator<N, Algo>::tol() const noexcept
    {
        return config_.tol;
    }
}
