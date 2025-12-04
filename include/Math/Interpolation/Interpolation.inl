/**
* Interpolation.inl
* Author: Alvaro Sanchez de Carlos
*/

#include "Utils/Log.hpp"
#include "Errors/Errors.hpp"  

#include <cstddef>
#include <string>
#include <format>

namespace uv
{
	template <typename T>
	::std::vector<T> d1PieceWise(
		const ::std::vector<T>& x,
		const ::std::vector<T>& y,
		bool extrapolateEnd
		)
	{
		// ---------- Sanity checks ------

		// Extract x and y array sizes
		const ::std::size_t xSize{ x.size() };
		const ::std::size_t ySize{ y.size() };

		// Require interpolation vectors to be the same size
		UV_REQUIRE(
			xSize == ySize,
			ErrorCode::InvalidArgument,
			"d1PieceWise: size mismatch — x vector has " + ::std::to_string(xSize) +
			" elements but y vector has " + ::std::to_string(ySize) + " elements"
		);

		// Require interpolation vectors to have more than one element
		UV_REQUIRE(
			xSize > 1,
			ErrorCode::InvalidArgument,
			"d1PieceWise: insufficient data — need at least 2 points, but got "
			+ ::std::to_string(xSize)
		);

		// Warn if x vector is not strictly monotonically increasing
		bool isMonotonic{ true };
		std::size_t violIndex{ 0 };

		// Check for monotonicity condition violations
		for (::std::size_t i = 1; i < xSize; ++i)
		{
			// strict: x[i] must be > x[i - 1]
			if (x[i] <= x[i - 1])
			{
				isMonotonic = false;
				violIndex = i;
				break;
			}
		}
		// Print a warning message
		UV_WARN(!isMonotonic,
			::std::format(
				"Monotonicity warning: x vector is not strictly increasing "
				"(violation at index {}: x[i-1] = {:.6g}, x[i] = {:.6g})",
				violIndex, x[violIndex - 1], x[violIndex]
			)
		);

		// ---------- Piecewise derivative calculation ------

		// Initialize results vector
		// NOTE: if the end elements will not be extrapolated, trim the results array by one
		::std::vector<T> derivatives(extrapolateEnd ? xSize : xSize - 1);

		// Calculate piecewise derivatives
		for (::std::size_t i = 0; i < xSize - 1; ++i)
		{
			derivatives[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
		}

		// ---------- Linear extrapolation of the last derivative ------

		if (extrapolateEnd) derivatives.back() = derivatives[xSize - 2];

		return derivatives;
	};
}
