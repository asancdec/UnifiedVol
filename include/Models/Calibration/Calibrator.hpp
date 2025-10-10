/**
* Calibrator.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once
#include "Models/Calibration/Config.hpp"

#include <nlopt.hpp>

#include <vector>
#include <array>
#include <cstddef>

template <std::size_t N>
class Calibrator
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

    std::vector<double> lowerBounds_;  // Lower parameter bounds
    std::vector<double> upperBounds_;  // Upper parameter bounds
    std::vector<double> initGuess_;    // Initial parameter guess

    Config<N> config_;                 // Calibration configuration (tolerances, limits)
    nlopt::opt opt_;                   // NLopt optimizer instance

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
    Calibrator() = delete;

    explicit Calibrator(std::array<double, N> initGuess,
                        std::array<double, N> lowerBounds,
                        std::array<double, N> upperBounds,
                        const Config<N>& config,
                        nlopt::algorithm algo);

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
    std::vector<double> optimize() noexcept;
 
    //--------------------------------------------------------------------------
    // Getters
    //--------------------------------------------------------------------------
    const double& eps() const noexcept;
};

#include "Calibrator.inl"