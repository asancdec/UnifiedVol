/**
* double.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

#include "Calibration/NLopt/Config.hpp"
#include "Utils/Aux/StopWatch.hpp"

#include <nlopt.hpp>
#include <array>
#include <cstddef>

namespace uv::cal::nlopt
{

    template <std::size_t N, ::nlopt::algorithm Algo>
    class Calibrator
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
        Config<N>         config_;
        ::nlopt::opt      opt_;
        utils::StopWatch  timer_;

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
        Calibrator() = delete;
        explicit Calibrator(const Config<N>& config);

        // Return new calibrator object with same settings
        Calibrator<N, Algo> fresh() const noexcept;

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
        Vector<double> optimize() noexcept;

        //--------------------------------------------------------------------------
        // Getters
        //--------------------------------------------------------------------------
        const double& eps() const noexcept;
        double tol() const noexcept;
    };
}

#include "Calibrator.inl"