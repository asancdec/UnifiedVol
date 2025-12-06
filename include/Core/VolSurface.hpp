/**
* VolSurface.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

#include "Core/SliceData.hpp"
#include "Core/MarketData.hpp"
#include "Utils/Types.hpp"

#include <cstddef>

namespace uv::core
{
    class VolSurface
    {
    private:

        //--------------------------------------------------------------------------
        // Member variables
        //--------------------------------------------------------------------------

        std::vector<SliceData> slices_;         // 2D volatility grid: vols[maturity][tenors slice]
        Vector<Real>    tenors_;         // Vector with different volatility surface tenors

        std::size_t numTenors_;                 // Number of tenors in years
        std::size_t numStrikes_;                // Number of strikes

    public:

        //--------------------------------------------------------------------------
        // Initialization
        //--------------------------------------------------------------------------

        VolSurface() = delete;
        explicit VolSurface(const Vector<Real>& mny,
            const Matrix<Real>& vols,
            const Vector<Real>& tenors,
            const MarketData& mktData);

        //--------------------------------------------------------------------------
        // Utilities
        //--------------------------------------------------------------------------

        // Print volatility surface on the console
        void printVol() const noexcept;

        // Print total variance surface on the console
        void printTotVar() const noexcept;

        // Print Black-Scholes Call price surface
        void printBSCall() const noexcept;

        // Return total variance matrix
        Matrix<Real> totVarMatrix() const noexcept;

        // Return the logFM matrix
        Matrix<Real> logFMMatrix() const noexcept;

        // Getters
        std::vector<SliceData>& slices() noexcept;
        const Vector<Real>& tenors() const noexcept;
        std::size_t numTenors() const noexcept;
        std::size_t numStrikes() const noexcept;
    };
}
