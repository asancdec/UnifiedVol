/**
* MathFunctions.inl
* Author: Alvaro Sanchez de Carlos
*/

#include "Utils/Log.hpp"
#include <cmath>
#include <numbers>

namespace uv
{
    template <typename T>
    inline ::std::complex<T> log1pComplex(const ::std::complex<T>& z) noexcept
    {
        // Use series only for small |z| to avoid cancellation
        if (::std::abs(z) < static_cast<T>(1e-2))
        {
            // Start at n=1 term
            ::std::complex<T> sum = z;    

            // z^n term starts at z^1
            ::std::complex<T> zpow = z;   

            // Expand Tayor Series
            for (int n = 2; n <= 50; ++n)
            {
                // z^n
                zpow *= z;       

                // z^n / n
                ::std::complex<T> term{ zpow / T(n) };

                // Early termination when incremental term negligible
                if (::std::abs(term) <= ::std::abs(sum) * ::std::numeric_limits<T>::epsilon())
                    break;

                // Determine sign using bitwise operation
                if ((n & 1) == 0) term = -term;  

                // Accumulate terms
                sum += term;
            }

            return sum;
        }

        // Fallback: direct computation for general case
        return ::std::log(::std::complex<T>(1) + z);
    }

    template <typename T>
    inline T normalCDF(T x) noexcept
    {
        return ::std::erfc(-x / ::std::sqrt(T(2.0))) * T(0.5);
    }

    template <typename T>
    inline T normalPDF(T x) noexcept
    {
        constexpr T invSqrt2Pi{ T(1.0) / ::std::sqrt(T(2.0) * ::std::numbers::pi_v<T>)};
        return invSqrt2Pi * ::std::exp(-T(0.5) * x * x);
    }

    template <typename T>
    inline T d1BS(T t,
        T r,
        T q,
        T vol,
        T S,
        T K) noexcept
    {
        return ::std::fma(t, (r - q + vol * vol * T(0.5)), ::std::log(S / K)) / (::std::sqrt(t) * vol);
    }

    template <typename T>
    inline T blackScholes(T t,
        T r,
        T q,
        T vol,
        T S,
        T K,
        bool isCall) noexcept
    {
        T d1{ uv::d1BS(t, r, q, vol, S, K)};
        T d2{ ::std::fma(-vol, ::std::sqrt(t), d1) };

        if (isCall)
        {
            return S * ::std::exp(-q * t) * uv::normalCDF(d1) - K * ::std::exp(-r * t) * uv::normalCDF(d2);
        }
        else
        {
            return K * ::std::exp(-r * t) * uv::normalCDF(-d2) - S * ::std::exp(-q * t) * uv::normalCDF(-d1);
        }
    }

    template <typename T>
    inline T vegaBS(T d1,
        T t,
        T q,
        T S) noexcept
    {
        return S * ::std::exp(-q * t) * uv::normalPDF(d1) * ::std::sqrt(t);
    }

    template <typename T>
    inline T volgaBS(T vega,
        T d1,
        T t,
        T vol) noexcept
    {
        T d2{ ::std::fma(-vol, ::std::sqrt(t), d1) };
        return vega * (d1 * d2) / vol;
    }
}