/**
* Calibrator.hpp
* Author: Alvaro Sanchez de Carlos
*/

#include "Utils/Aux/Helpers.hpp"
#include "Utils/IO/ConsoleRedirect.hpp"
#include "Utils/IO/Log.hpp"

#include <memory>
#include <iomanip>   
#include <format>    
#include <algorithm>

namespace uv::cal::ceres
{
    template <std::size_t N, typename Policy>
    Calibrator<N, Policy>::Calibrator(const Config<N>& config) : config_(config) {}

    template <std::size_t N, typename Policy>
    void Calibrator<N, Policy>::setGuessBounds(const std::array<double, N>& initGuess,
        const std::array<double, N>& lowerBounds,
        const std::array<double, N>& upperBounds) noexcept
    {
        // Set member variables
        x_ = initGuess;
        lowerBounds_ = lowerBounds;
        upperBounds_ = upperBounds;

        // Clamp initial guess within upper and lower bounds
        utils::clamp<N>(x_, lowerBounds_, upperBounds_, config_.paramNames);

        // Set initial guess
        problem_.AddParameterBlock(x_.data(), static_cast<int>(N));

        // Apply per-parameter bounds on that same block
        for (int i = 0; i < static_cast<int>(N); ++i)
        {
            problem_.SetParameterLowerBound(x_.data(), i, lowerBounds_[i]);
            problem_.SetParameterUpperBound(x_.data(), i, upperBounds_[i]);
        }
    }

    template <std::size_t N, typename Policy>
    void Calibrator<N, Policy>::addAnalyticResidual(std::unique_ptr<::ceres::CostFunction> cf) noexcept
    {
        problem_.AddResidualBlock
        (
            cf.release(),                                        
            Policy::makeLoss(config_.lossScale).release(),       
            x_.data()
        );
    }

    template <std::size_t N, typename Policy>
    std::array<double, N> Calibrator<N, Policy>::optimize()
    {   
        // Set calibration options
        ::ceres::Solver::Options options;
        options.trust_region_strategy_type = Policy::trustRegionStrategy; 
        options.linear_solver_type         = Policy::linearSolver;         
        options.max_num_iterations         = config_.maxEval;            
        options.function_tolerance         = config_.functionTol;          
        options.parameter_tolerance        = config_.paramTol;           
        options.gradient_tolerance         = config_.gradientTol;        
        options.num_threads                = std::max(1u, std::thread::hardware_concurrency()); 
        
       
        // Capture Ceres' progress output and redirect it to the unified UV logger
        ::ceres::Solver::Summary summary;
        {
            // Enable live Ceres iteration table only if verbose mode is on
            utils::ConsoleRedirect capture;
            options.minimizer_progress_to_stdout = config_.verbose;

            // Solve the problem
            ::ceres::Solve(options, &problem_, &summary);

            // Print final solver report if verbose mode is on
            if (config_.verbose)  UV_INFO(summary.FullReport());
        }

        // Warn if upper or lower bounds are touched
        utils::warnBoundsHit
        (
            x_,
            lowerBounds_,
            upperBounds_,
            config_.paramNames
        );

        // Log calibration results 
        utils::logResults
        (
            x_,                                                     // Parameters
            config_.paramNames,                                     // Parameter names
            summary.final_cost * 2.0,                               // SSE
            summary.iterations.size(),                              // Iterations
            summary.total_time_in_seconds * 1000.0,                 // Elapsed [ms]
            (summary.termination_type == ::ceres::CONVERGENCE ||
                summary.termination_type == ::ceres::USER_SUCCESS)  // Success flag
        );

        // Return calibrated parameters
        return x_; 
    }
}

 