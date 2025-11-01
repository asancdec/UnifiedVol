/**
* CalibratorCeres.hpp
* Author: Alvaro Sanchez de Carlos
*/

#include "Math/Calibration/Utils/CalibratorUtils.hpp"

namespace uv
{
	template <::std::size_t N>
	CalibratorCeres<N>::CalibratorCeres(const CalibratorConfig<N>& config) : config_(config) {}

    template <::std::size_t N>
    void CalibratorCeres<N>::setGuessBounds(const ::std::array<double, N>& initGuess,
        const ::std::array<double, N>& lowerBounds,
        const ::std::array<double, N>& upperBounds) noexcept
    {
        // Set member variables
        x_ = initGuess;
        lowerBounds_ = lowerBounds;
        upperBounds_ = upperBounds;

        // Clamp initial guess within upper and lower bounds
        uv::clamp<N>(x_, lowerBounds_, upperBounds_, config_.paramNames);

        // Set initial guess
        problem_.AddParameterBlock(x_.data(), static_cast<int>(N));

        // Apply per-parameter bounds on that same block
        for (int i = 0; i < static_cast<int>(N); ++i)
        {
            problem_.SetParameterLowerBound(x_.data(), i, lowerBounds_[i]);
            problem_.SetParameterUpperBound(x_.data(), i, upperBounds_[i]);
        }
    }
    
}