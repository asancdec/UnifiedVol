/*
* LocalVol.cpp
* Author: Alvaro Sanchez de Carlos
*/


#include "Models/LocalVol/LocalVol.hpp"
#include "Utils/Aux/Errors.hpp"    "
#include "Math/Interpolation/Interpolation.hpp"
#include "Utils/Aux/Helpers.hpp"
#include "Models/SVI/SVI.hpp"
#include "Core/SliceData.hpp"

#include <cstddef>
#include <string>
#include <cmath>
#include <iostream>

using namespace uv;

namespace uv::localvol
{
	VolSurface buildSurface(const VolSurface& volSurface,
		const Vector<models::svi::SVISlice>& sviSlices)
	{

		// ---------- Check matching dimensions ----------
		
		const std::size_t numTenors{ volSurface.numTenors() };
		const std::size_t numStrikes{ volSurface.numStrikes() };
		const std::size_t numSlices{ sviSlices.size() };

		// Throw if dimensions do not match
		UV_REQUIRE
		(
			numTenors == numSlices,
			ErrorCode::InvalidArgument,
			"build: number of tenors (" + std::to_string(numTenors) +
			") does not match number of SVI slices (" +
			std::to_string(numSlices) + ")"
		);

		// ---------- Check matching tenors ----------

		const Vector<double>& tenors{ volSurface.tenors() };

		for (std::size_t i = 0; i < numTenors; ++i)
		{
			const double volTenor{ tenors[i] };
			const double sviTenor{ sviSlices[i].T };

			// Throw 
			UV_REQUIRE
			(
				std::abs(volTenor - sviTenor) < 1e-12,	// tolerance
				ErrorCode::InvalidArgument, 
				"localvol::build: tenor mismatch at index " 
				+ std::to_string(i) +
				" — volSurface tenor = " + std::to_string(volTenor) +
				", SVI tenor = " + std::to_string(sviTenor)
			);
		}

		// ---------- Calculate dw/dT matrix ------

		const Matrix<double> dwdT
		{
			detail::dwdT
			(
				tenors, 
				volSurface.totVarMatrix()
			)
		};


		// Copy the volatility surface object
		VolSurface lvVolSurf{ volSurface };


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
			const models::svi::SVISlice& sviSlice{ sviSlices[i]};

			// Extract parameters
			const double a{ sviSlice.a };
			const double b{ sviSlice.b };
			const double rho{ sviSlice.rho };
			const double m{ sviSlice.m };
			const double sigma{ sviSlice.sigma };
			const double T{ sviSlice.T };

			for (std::size_t j = 0; j < numStrikes; ++j)
			{
				localTotVar[j] = dwdT[i][j] * T / models::svi::gk(a, b, rho, m, sigma, logFM[j]);
			}

			// Set local total variance
			volSlice.setWT(localTotVar);
		}

		return lvVolSurf;
	}

} // namespace uv::localvol

namespace uv::localvol::detail
{	
	using namespace uv;

	Matrix<double> dwdT(const std::vector<double>& tenors, 
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
		const std::vector<std::vector<double>> totVarT{ utils::transposeMatrix(totVarMatrix) };

		// Calculate dw/dT
		for (std::size_t i = 0; i < numStrikes; ++i)
		{
			// Calculates piecewise derivatives
			dwdT[i] = pchipDerivatives(tenors, totVarT[i]);
		}

		// Transpose the dw/dT matrix
		// Rows of strikes and columns of tenors
		return utils::transposeMatrix(dwdT);
	}

} // namespace uv::localvol