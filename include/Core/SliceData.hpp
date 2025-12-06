/**
* SliceData.hpp
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

#include "Core/MarketData.hpp"
#include "Utils/Types.hpp"

namespace uv::core
{
    class SliceData
    {
    private:

        //--------------------------------------------------------------------------
        // Member variables
        //--------------------------------------------------------------------------

        Real T_;                         // Slice maturity
        Real F_;                         // Forward price
        Vector<Real> mny_;          // Plain moneyness (K/S)
        Vector<Real> logFM_;        // Log-forward moneyness log(K/F)
        Vector<Real> vol_;          // Volatilities 
        Vector<Real> wT_;           // Total variance (vol²T) 
        Vector<Real> K_;            // Strikes vector
        Vector<Real> callBS_;       // Value of European call options
        MarketData mktData_;               // Market data struct
        std::size_t numStrikes_;           // Number of strikes

    public:

        //--------------------------------------------------------------------------
        // Initialization
        //--------------------------------------------------------------------------   
        SliceData() = delete;
        explicit SliceData(Real T,
            const Vector<Real>& mny,
            const Vector<Real>& vol,
            const MarketData& mktData);

        //--------------------------------------------------------------------------
        // Math functions
        //--------------------------------------------------------------------------

        // Determine minimum total variance
        Real minWT() const noexcept;

        // Determine maximum total variance
        Real maxWT() const noexcept;

        // Determine minimum log-forward moneyness
        Real minLogFM() const noexcept;

        // Determine maximum log-forward moneyness
        Real maxLogFM() const noexcept;

        // Return ATM total variance
        Real atmWT() const noexcept;

        //--------------------------------------------------------------------------
        // Getters
        //--------------------------------------------------------------------------

        Real T() const noexcept;
        Real F() const noexcept;
        std::size_t numStrikes() const noexcept;
        Real r() const noexcept;
        const Vector<Real>& mny() const noexcept;
        const Vector<Real>& logFM() const noexcept;
        const Vector<Real>& vol() const noexcept;
        const Vector<Real>& wT() const noexcept;
        const Vector<Real>& K() const noexcept;
        const Vector<Real>& callBS() const noexcept;


        //--------------------------------------------------------------------------
        // Setters
        //--------------------------------------------------------------------------
        void setWT(const Vector<Real>& wT);
        void setCallBS(const Vector<Real>& callBS);
    };
}
