/**
* CalibratorNLopt.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

#include "Math/Calibration/NLopt/NLoptConfig.hpp"
#include "Utils/StopWatch/StopWatch.hpp"

#include <nlopt.hpp>
#include <vector>
#include <array>
#include <cstddef>

namespace uv
{

    template <std::size_t N, nlopt::algorithm Algo>
    class CalibratorNLopt
    {
    private:

        //--------------------------------------------------------------------------
        // Type aliases
        //--------------------------------------------------------------------------
        using NloptFunction = double (*)(unsigned, const double*, double*, void*);


        //--------------------------------------------------------------------------
        // Member variables
        //--------------------------------------------------------------------------	
        // Config & engine
        NLoptConfig<N>      config_;
        nlopt::opt        opt_;
        StopWatch           timer_;

        // Bounds & guess
        std::array<double, N> lowerBounds_;
        std::array<double, N> upperBounds_;
        std::array<double, N> initGuess_;

        // Objective state
        NloptFunction userFn_;
        void* userData_;
        unsigned iterCount_;   

        //--------------------------------------------------------------------------
        // Static variables
        //--------------------------------------------------------------------------	
        static double ObjectiveThunk(unsigned n, const double* x, double* grad, void* p) noexcept;

    public:

        //--------------------------------------------------------------------------
        // Initialization
        //--------------------------------------------------------------------------
        // Constructors
        CalibratorNLopt() = delete;
        explicit CalibratorNLopt(const NLoptConfig<N>& config);

        // Return new calibrator object with same settings
        CalibratorNLopt<N, Algo> fresh() const noexcept;

        //--------------------------------------------------------------------------
        // Calibration
        //--------------------------------------------------------------------------	
        // Set Initial Guess and Bounds
        void setGuessBounds(std::array<double, N> initGuess,
            std::array<double, N> lowerBounds,
            std::array<double, N> upperBounds) noexcept;

        // Add inequality constraints
        void addInequalityConstraint(NloptFunction c, void* data) noexcept;

        // Set objective function	
        void setMinObjective(NloptFunction f, void* data) noexcept;

        // Run calibration	
        std::vector<double> optimize() noexcept;

        //--------------------------------------------------------------------------
        // Getters
        //--------------------------------------------------------------------------
        const double& eps() const noexcept;
        double tol() const noexcept;
    };
}

#include "CalibratorNLopt.inl"