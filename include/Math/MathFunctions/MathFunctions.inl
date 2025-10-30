/**
* MathFunctions.inl
* Author: Alvaro Sanchez de Carlos
*/

#pragma once

#include <limits>
#include <cmath>

namespace uv
{
    template <typename T>
    inline std::complex<T> log1pComplex(const std::complex<T>& z) noexcept
    {
        // Use series only for small |z| to avoid cancellation
        if (std::abs(z) < static_cast<T>(1e-2))
        {
            // Start at n=1 term
            std::complex<T> sum = z;    

            // z^n term starts at z^1
            std::complex<T> zpow = z;   

            // Expand Tayor Series
            for (int n = 2; n <= 50; ++n)
            {
                // z^n
                zpow *= z;       

                // z^n / n
                std::complex<T> term{ zpow / T(n) };

                // Early termination when incremental term negligible
                if (std::abs(term) <= std::abs(sum) * std::numeric_limits<T>::epsilon())
                    break;

                // Determine sign using bitwise operation
                if ((n & 1) == 0) term = -term;  

                // Accumulate terms
                sum += term;
            }

            return sum;
        }

        // Fallback: direct computation for general case
        return std::log(std::complex<T>(1) + z);
    }
}