/**
* Interpolation.inl
* Author: Alvaro Sanchez de Carlos
*/

#include "Utils/IO/Log.hpp"
#include "Utils/Aux/Errors.hpp"      

#include <cstddef>
#include <cmath>
#include <string>
#include <format>
#include <algorithm>
#include <iostream>

namespace uv
{
	template <std::floating_point T>
	T interpolateCubicHermiteSpline(const T x,
		const std::vector<T>& xs,
		const std::vector<T>& ys,
		const std::vector<T>& dydx
	)
	{
		// ---------- Check matching dimensions ----------

		const std::size_t xsSize{ xs.size() };
		const std::size_t ysSize{ ys.size() };
		const std::size_t dSize{ dydx.size() };

		// Throw if sizes do not match
		UV_REQUIRE(
			(xsSize == ysSize) && (xsSize == dSize),
			ErrorCode::InvalidArgument,
			"interpolateCubicHermiteSpline: size mismatch — "
			"xs vector has " + std::to_string(xsSize) +
			" elements, ys vector has " + std::to_string(ysSize) +
			" elements, dydx vector has " + std::to_string(dSize) + " elements"
		);

		// ---------- Special cases  ----------

		if (xsSize == 0) return T(0);   // No elements
		if (xsSize == 1) return ys[0];  // One element

		// ---------- Check strictly monotonically increasing x values ----------

		bool isMonotonic{ true };
		std::size_t violIndex{ 0 };

		for (std::size_t i = 1; i < xsSize; ++i)
		{
			if (xs[i] <= xs[i - 1])
			{
				isMonotonic = false;
				violIndex = i;
				break;
			}
		}

		// Throw if not strictly monotonically increasing
		UV_REQUIRE(
			isMonotonic,
			ErrorCode::InvalidArgument,
			std::format(
				"interpolateCubicHermiteSpline: x vector must be strictly increasing "
				"(violation at index {}: xs[i-1] = {:.6g}, xs[i] = {:.6g})",
				violIndex, xs[violIndex - 1], xs[violIndex]
			)
		);

		// ---------- Out of bounds ----------

		const T xsMin{ xs.front() };
		const T xsMax{ xs.back() };

		// To the left
		if (x < xsMin)
		{
			UV_WARN(
				true,
				std::format(
					"interpolateCubicHermiteSpline: x = {:.4f} out of bounds -> clamped to {:.4f} "
					"(xsmin = {:.4f}, xsmax = {:.4f})",
					x, xsMin, xsMin, xsMax
				)
			);

			// Flat extrapolation
			return ys.front();
		}

		// To the right
		if (x > xsMax)
		{
			UV_WARN(
				true,
				std::format(
					"interpolateCubicHermiteSpline: x = {:.4f} out of bounds -> clamped to {:.4f} "
					"(xsmin = {:.4f}, xsmax = {:.4f})",
					x, xsMax, xsMin, xsMax
				)
			);

			// Flat extrapolation
			return ys.back();
		}

		// ---------- Calculate step sizes and secant slopes ----------

		const std::size_t numSteps{ xsSize - 1 };
		std::vector<T> h(numSteps);                   // Step sizes
		std::vector<T> S(numSteps);			          // Secant slopes

		for (std::size_t i = 0; i < numSteps; ++i)
		{
			const T hi{ xs[i + 1] - xs[i] };
			h[i] = hi;
			S[i] = (ys[i + 1] - ys[i]) / hi;
		}

		// ---------- Second and third order coefficients ----------

		std::vector<T> c2s(numSteps);      // Second order
		std::vector<T> c3s(numSteps);	   // Third order

		for (std::size_t i = 0; i < numSteps; ++i)
		{
			const T c1{ dydx[i]};          // Tangent slope
			const T m{ S[i] };		       // Secant slope
			const T invH{ T(1) / h[i] };
			const T common{ c1 + dydx[i + 1] - T(2) * m };
			
			c2s[i] = (m - c1 - common) * invH;
		    c3s[i] = (common * invH * invH);
		}

		// ---------- Search for index ----------

		// Distance from first element to iterator at upper bound
		auto it = std::upper_bound(xs.begin(), xs.end(), x);
		std::size_t idx{ static_cast<std::size_t>(std::distance(xs.begin(), it)) };

		// NOTE: indices start from 0
		if (idx >= xsSize) idx = numSteps - 1; 
		else --idx;						

		// ---------- Calculate interpolated value ----------

		const T dx{ x - xs[idx] };

		return ys[idx]
			+ dydx[idx] * dx
			+ c2s[idx] * dx * dx
			+ c3s[idx] * dx * dx * dx;
	}

	template <std::floating_point T>
	std::vector<T> pchipDerivatives(const std::vector<T>& xs,
		const std::vector<T>& ys)
	{
		// ---------- Check matching dimensions ----------
		
		const std::size_t xsSize{ xs.size() };
		const std::size_t ysSize{ ys.size() };

		// Throw if sizes do not match
		UV_REQUIRE(
			xsSize == ysSize,
			ErrorCode::InvalidArgument,
			"pchipDerivatives: size mismatch — xs vector has " + std::to_string(xsSize) +
			" elements but ys vector has " + std::to_string(ysSize) + " elements"
		);

		// ---------- Check more than one element ----------

		// Throw if only one element
		UV_REQUIRE(
			xsSize > 1,
			ErrorCode::InvalidArgument,
			"pchipDerivatives: insufficient data — need at least 2 points, but got "
			+ std::to_string(xsSize)
		);

		// ---------- Check strictly monotonically increasing x values ----------

		bool isMonotonic{ true };
		std::size_t violIndex{ 0 };

		for (std::size_t i = 1; i < xsSize; ++i)
		{
			if (xs[i] <= xs[i - 1])
			{
				isMonotonic = false;
				violIndex = i;
				break;
			}
		}

		// Throw if not strictly monotonically increasing
		UV_REQUIRE(
			isMonotonic,
			ErrorCode::InvalidArgument,
			std::format(
				"pchipDerivatives: xs vector must be strictly increasing "
				"(violation at index {}: xs[i-1] = {:.6g}, xs[i] = {:.6g})",
				violIndex, xs[violIndex - 1], xs[violIndex]
			)
		);

		// ---------- Special case: two elements  ----------

		if (xsSize == 2)
		{
			// Compute and return single secant slope
			const T S{ (ys[1] - ys[0]) / (xs[1] - xs[0]) };
			return std::vector<T>{ S, S };
		}

		// ---------- Calculate step sizes and secant slopes ----------
	
		const std::size_t numSteps{xsSize - 1};
		std::vector<T> h(numSteps);                   // Step sizes
		std::vector<T> S(numSteps);			          // Secant slopes

		for (std::size_t i = 0; i < numSteps; ++i)
		{	
			const T hi{ xs[i + 1] - xs[i] };
			h[i] = hi;
			S[i] = (ys[i + 1] - ys[i]) / hi;
		}

		// ---------- Calculate middle derivatives ----------

		std::vector<T> dydx(xsSize);          // All initialized to zero

		for (std::size_t i = 1; i < numSteps; ++i)
		{	
			const T S1{ S[i - 1] };
			const T S2{ S[i] };
			const T S1timesS2{ S1 * S2 };

			// Check slope signs
			// NOTE: else, the derivative is set to zero
			if (S1timesS2 > T(0))
			{
				const T h1{ h[i - 1] };
				const T h2{ h[i] };

				// Calculate weight
				const T weight{ (h1 + T(2) * h2) / (T(3) * (h1 + h2)) };

				// Calculate tangent slope (derivative)
				dydx[i] = S1timesS2 / (weight * S2 + (T(1) - weight) * S1);
			}
		}

		// ---------- Calculate edge derivatives ----------

		// Left endpoint
		dydx[0] = pchipEndpointSlope<T>(
			h[0],      // h1 = xs[1] - xs[0]
			h[1],      // h2 = xs[2] - xs[1]
			S[0],      // S1 = (ys[1] - ys[0]) / h1
			S[1]       // S2 = (ys[2] - ys[1]) / h2
		);

		// Right endpoint
		dydx.back() = pchipEndpointSlope<T>(
			h[numSteps - 1],   // h1 = xs[n-1] - xs[n-2]
			h[numSteps - 2],   // h2 = xs[n-2] - xs[n-3]
			S[numSteps - 1],   // S1 = (ys[n-1] - ys[n-2]) / h1
			S[numSteps - 2]    // S2 = (ys[n-2] - ys[n-3]) / h2
		);

		return dydx;
	}

	template <std::floating_point T>
	T pchipEndpointSlope(const T h1,
		const T h2,
		const T S1,
		const T S2
	) noexcept
	{
		// ---------- Derivative ----------

		T d{ ((T(2) * h1 + h2) * S1 - h1 * S2) / (h1 + h2) };

		// ---------- Shape-preserving conditions ----------
		// (1) If d points against the local trend, clamp to zero
		if (std::signbit(d) != std::signbit(S1))
		{
			d = T(0);
		}
		// (2) If adjacent secant slopes disagree AND magnitude too large, clamp
		else if ((std::signbit(S1) != std::signbit(S2)) && std::abs(d) > T(3) * std::abs(S1))
		{
			d = T(3) * S1;
		}
		return d;
	}
}
