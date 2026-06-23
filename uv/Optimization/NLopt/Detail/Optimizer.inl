// SPDX-License-Identifier: Apache-2.0

#include "Base/Errors/Errors.hpp"
#include "Base/Macros/Require.hpp"
#include "Optimization/Helpers.hpp"
#include "Optimization/NLopt/Detail/MapAlgorithm.hpp"
#include "Optimization/NLopt/Detail/NLoptStatus.hpp"

#include <nlopt.hpp>
#include <string_view>

namespace uv::opt::nlopt
{
// NLopt's public C/C++ callback interface requires function pointers carrying
// state through void*. These exact signatures are an ABI boundary, not type erasure by
// choice.
template <std::size_t N, Algorithm Algo>
Optimizer<N, Algo>::Optimizer(const Config<N>& config)
    : config_(config),
      opt_(detail::toNlopt(Algo), N)
{
    opt_.set_ftol_rel(config_.ftolRel);
    opt_.set_maxeval(config_.maxEval);
    opt_.set_exceptions_enabled(false);
}

template <std::size_t N, Algorithm Algo>
Optimizer<N, Algo> Optimizer<N, Algo>::fresh() const
{
    return Optimizer<N, Algo>{config_};
}

template <std::size_t N, Algorithm Algo> void Optimizer<N, Algo>::setGuessBounds(
    std::array<double, N> initGuess,
    std::array<double, N> lowerBounds,
    std::array<double, N> upperBounds
)
{

    clampBounds(initGuess, lowerBounds, upperBounds);

    initGuess_ = initGuess;
    lowerBounds_ = lowerBounds;
    upperBounds_ = upperBounds;

    opt_.set_lower_bounds(Vector<double>(lowerBounds_.begin(), lowerBounds_.end()));
    opt_.set_upper_bounds(Vector<double>(upperBounds_.begin(), upperBounds_.end()));
}

template <std::size_t N, Algorithm Algo> double Optimizer<N, Algo>::objectiveThunk(
    unsigned n,
    const double* x,
    double* grad,
    void* p // NOSONAR -- NLopt's C callback ABI requires an opaque void pointer.
) noexcept
{
    auto* self = static_cast<Optimizer<N, Algo>*>(p);
    ++self->iterCount_;
    return self->userFn_ ? self->userFn_(n, x, grad, self->userData_) : 0.0;
}

template <std::size_t N, Algorithm Algo>
void Optimizer<N, Algo>::addInequalityConstraint(NloptFunction c, void* data) // NOSONAR
{
    opt_.add_inequality_constraint(c, data, config_.tol);
}

template <std::size_t N, Algorithm Algo> void
Optimizer<N, Algo>::addInequalityMConstraint( // NOSONAR -- Exact NLopt C callback ABI.
    std::size_t m,
    NloptMFunction c,
    void* data
)
{
    Vector<double> tol(m, config_.tol);
    opt_.add_inequality_mconstraint(c, data, tol);
}

template <std::size_t N, Algorithm Algo>
void Optimizer<N, Algo>::setMinObjective(NloptFunction f, void* data) // NOSONAR
{
    iterCount_ = 0U;
    userFn_ = f;
    userData_ = data;

    if (config_.verbose)
    {
        opt_.set_min_objective(&Optimizer<N, Algo>::objectiveThunk, this);
    }
    else
    {
        opt_.set_min_objective(userFn_, userData_);
    }
}

template <std::size_t N, Algorithm Algo> Vector<double> Optimizer<N, Algo>::optimize()
{

    Vector<double> x(initGuess_.cbegin(), initGuess_.cend());
    double sse{0.0};

    timer_.StartStopWatch();

    ::nlopt::result successCode = opt_.optimize(x, sse);

    timer_.StopStopWatch();

    if (config_.verbose)
    {
        warnBoundsHit(x, lowerBounds_, upperBounds_);

        logResults(
            x,
            config_.paramNames,
            sse,
            iterCount_,
            timer_.GetTime<std::milli>(),
            (successCode > ::nlopt::FAILURE) ||
                (successCode == ::nlopt::ROUNDOFF_LIMITED),
            detail::toString(successCode)
        );
    }

    return x;
}

template <std::size_t N, Algorithm Algo>
void Optimizer<N, Algo>::setUserValue(double v) noexcept
{
    userValue_ = v;
}

template <std::size_t N, Algorithm Algo>
const double& Optimizer<N, Algo>::eps() const noexcept
{
    return config_.eps;
}

template <std::size_t N, Algorithm Algo> double Optimizer<N, Algo>::tol() const noexcept
{
    return config_.tol;
}

template <std::size_t N, Algorithm Algo>
const double& Optimizer<N, Algo>::userValue() const
{
    if (!userValue_.has_value()) [[unlikely]]
    {
        errors::raise(errors::ErrorCode::InvalidState, "userValue_ must be set");
    }

    return *userValue_;
}
} // namespace uv::opt::nlopt
