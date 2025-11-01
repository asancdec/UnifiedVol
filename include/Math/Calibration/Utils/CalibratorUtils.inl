/**
* CalibratorUtils.inl
* Author: Alvaro Sanchez de Carlos
*/

#include "Utils/Log.hpp"
#include <format>
#include <algorithm>

namespace uv
{
	template <::std::size_t N>
	inline void clamp(::std::array<double, N>& initGuess,
		const ::std::array<double, N>& lowerBounds,
		const ::std::array<double, N>& upperBounds,
        const ::std::array<::std::string_view, N>& paramNames) noexcept
	{
        for (::std::size_t i = 0; i < initGuess.size(); ++i)
        {
            const double before = initGuess[i];
            const double after = ::std::clamp(before, lowerBounds[i], upperBounds[i]);

            UV_WARN(after != before,
                ::std::format("Calibrator: parameter '{}' initial guess = {:.4f} "
                    "out of bounds -> clamped to {:.4f} "
                    "(lb = {:.4f}, ub = {:.4f})",
                    paramNames[i], before, after,
                    lowerBounds[i], upperBounds[i]));

            initGuess[i] = after;
        }
	}

}
