/**
* CalibratorNLopt.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once
#include "Math/Calibration/CalibratorConfig.hpp"
#include <nlopt.hpp>
#include <vector>
#include <array>
#include <cstddef>

namespace uv
{
    template <::std::size_t N>
    class CalibratorNLopt
    {
    private:

        //--------------------------------------------------------------------------
        // Type aliases
        //--------------------------------------------------------------------------
        using NloptFunction = double (*)(unsigned, const double*, double*, void*);

        //--------------------------------------------------------------------------
        // Internal data structures
        //--------------------------------------------------------------------------

        struct ObjWrapCtx
        {
            NloptFunction fn{};   // User-provided objective callback
            void* user{};         // Opaque user data passed to callback
            unsigned* iter{};     // Pointer to evaluation counter
        };

        //--------------------------------------------------------------------------
        // Member variables
        //--------------------------------------------------------------------------	
        CalibratorConfig<N> config_;       // Calibration configuration
        ::nlopt::opt opt_;                   // NLopt optimizer instance
        ::nlopt::algorithm algo_;            // Nlopt algorithm used

        ::std::vector<double> lowerBounds_;  // Lower parameter bounds
        ::std::vector<double> upperBounds_;  // Upper parameter bounds
        ::std::vector<double> initGuess_;    // Initial parameter guess

        ObjWrapCtx objCtx_{};              // Wrapped objective context used internally
        unsigned iterCount_{ 0 };          // Number of objective evaluations performed

        //--------------------------------------------------------------------------
        // Static variables
        //--------------------------------------------------------------------------	
        static double ObjectiveThunk(unsigned n, const double* x, double* grad, void* p) noexcept;

    public:

        //--------------------------------------------------------------------------
        // Initialization
        //--------------------------------------------------------------------------
        CalibratorNLopt() = delete;
        explicit CalibratorNLopt(const CalibratorConfig<N>& config,
            ::nlopt::algorithm algo);

        //--------------------------------------------------------------------------
        // Cloning
        //--------------------------------------------------------------------------

        // Return new calibrator object with same settings as the current instance
        CalibratorNLopt<N> fresh() const noexcept;

        //--------------------------------------------------------------------------
        // Set Initial Guess and Bounds
        //--------------------------------------------------------------------------	
        void setGuessBounds(::std::array<double, N> initGuess,
            ::std::array<double, N> lowerBounds,
            ::std::array<double, N> upperBounds) noexcept;

        //--------------------------------------------------------------------------
        // Add inequality constraints
        //--------------------------------------------------------------------------	
        void addInequalityConstraint(NloptFunction c, void* data) noexcept;

        //--------------------------------------------------------------------------
        // Set objective function
        //--------------------------------------------------------------------------	
        void setMinObjective(NloptFunction f, void* data) noexcept;

        //--------------------------------------------------------------------------
        // Run calibration
        //--------------------------------------------------------------------------	
        ::std::vector<double> optimize() noexcept;

        //--------------------------------------------------------------------------
        // Getters
        //--------------------------------------------------------------------------
        const double& eps() const noexcept;
        double tol() const noexcept;
    };
}

#include "CalibratorNLopt.inl"