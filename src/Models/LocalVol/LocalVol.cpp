/*
* LocalVol.cpp
* Author: Alvaro Sanchez de Carlos
*/

#include "Models/LocalVol/LocalVol.hpp"
#include "Errors/Errors.hpp"
#include "Math/Interpolation/Interpolation.hpp"
#include "Utils/Data/Aux/Aux.hpp"
#include "Models/SVI/SVI.hpp"
#include "Core/SliceData.hpp"

#include <cstddef>
#include <string>
#include <cmath>
#include <iostream>

using namespace uv;

namespace uv::local_vol
{
	VolSurface build(const VolSurface& volSurface, const std::vector<SVISlice>& sviSlices)
	{
		// ---------- Sanity checks ------
		
		// Extract the number of tenors
		const std::size_t numTenors{ volSurface.numMaturities() };

		// Extract the number of strikes
		const std::size_t numStrikes{ volSurface.numStrikes() };

		// Ensure matching dimensions
		UV_REQUIRE(
			numTenors == sviSlices.size(),
			ErrorCode::InvalidArgument,
			"LocalVol: number of maturities (" + std::to_string(numTenors) +
			") does not match number of SVI slices (" +
			std::to_string(sviSlices.size()) + ")"
		);

		// Copy the volatility surface object
		VolSurface lvVolSurf{ volSurface };

		// Extract volatility surface tenors
		const std::vector<double>& tenors{ lvVolSurf.maturities() };

		// Ensure volatility surface has the same maturities as the SVI parameters
		for (std::size_t i = 0; i < numTenors; ++i)
		{
			const double volMat{ tenors[i] };
			const double sviMat{ sviSlices[i].T };

			UV_REQUIRE(
				std::abs(volMat - sviMat) < 1e-12,
				ErrorCode::InvalidArgument,
				"LocalVol::build: maturity mismatch at index " + std::to_string(i) +
				" — volSurface maturity = " + std::to_string(volMat) +
				", SVI maturity = " + std::to_string(sviMat)
			);
		}

		// ---------- Calculate local total variance ------

		// Calculate dw/dT matrix
		const std::vector<std::vector<double>> dwdT{detail::dwdT(tenors, lvVolSurf.totVarMatrix())};

		// Extract volatility surface slices (will be modified in place)
		std::vector<SliceData>& volSlices{ lvVolSurf.slices()};

		// Calculate and set local total volatility variance
		for (std::size_t i = 0; i < numTenors; ++i)
		{
			// Temporary vector to store local total variance calculations
			std::vector<double> localTotVar(numStrikes);

			// Extract volatility surface slice
			SliceData& volSlice{volSlices[i]};

			// Extract logFM
			const std::vector<double>& logFM {volSlice.logFM()};

			// Extract SVI slice data
			const SVISlice& sviSlice{ sviSlices[i]};

			// Extract parameters
			const double a{ sviSlice.a };
			const double b{ sviSlice.b };
			const double rho{ sviSlice.rho };
			const double m{ sviSlice.m };
			const double sigma{ sviSlice.sigma };
			const double T{ sviSlice.T };

			for (std::size_t j = 0; j < numStrikes; ++j)
			{
				localTotVar[j] = dwdT[i][j] * T / svi::gk(a, b, rho, m, sigma, logFM[j]);

				//std::cout << dwdT[i][j] << "\n";
			}

			// Set local total variance
			volSlice.setWT(localTotVar);
		}

		return lvVolSurf;
	}

} // namespace uv::local_vol

namespace uv::local_vol::detail
{	
	using namespace uv;

	std::vector<std::vector<double>> dwdT(const std::vector<double>& tenors, 
		const std::vector<std::vector<double>>& totVarMatrix) noexcept
	{
		// Extract the number of tenors
		const std::size_t numTenors{ tenors.size()};

		// Extract the number of strikes
		const std::size_t numStrikes{ totVarMatrix[0].size()};

		// Generate dw/dT matrix to store results
		std::vector<std::vector<double>> dwdT(numStrikes, std::vector<double>(numTenors));

		// Transpose total variance matrix 
		// Rows of tenors and columns of strikes
		const std::vector<std::vector<double>> totVarT{ transposeMatrix(totVarMatrix) };

		// Calculate dw/dT
		for (std::size_t i = 0; i < numStrikes; ++i)
		{
			// Calculates piecewise derivatives
			dwdT[i] = d1PieceWise<double>(tenors, totVarT[i], /*extrapolateEnd=*/true);
		}

		// Transpose the dw/dT matrix
		// Rows of strikes and columns of tenors
		return transposeMatrix(dwdT);
	}

} // namespace uv::local_vol