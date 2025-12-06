/*
* LocalVol.cpp
* Author: Alvaro Sanchez de Carlos
*/


#include "Models/LocalVol/LocalVol.hpp"
#include "Utils/Aux/Errors.hpp"
#include "Math/Interpolation.hpp"
#include "Utils/Aux/Helpers.hpp"
#include "Models/SVI/SVI.hpp"
#include "Core/SliceData.hpp"

#include <cstddef>
#include <string>
#include <cmath>
#include <iostream>

using namespace uv;

namespace uv::models::localvol
{
	core::VolSurface buildSurface(const core::VolSurface& volSurface,
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

		const Vector<Real>& tenors{ volSurface.tenors() };

		for (std::size_t i = 0; i < numTenors; ++i)
		{
			const Real volTenor{ tenors[i] };
			const Real sviTenor{ sviSlices[i].T };

			// Throw if suface and SVI tenors do not match
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

		// ---------- Calculate local total variance matrix ----------

		const Matrix<Real> localtoVar
		{ 
			detail::localTotVar
			(
				tenors,                          
				volSurface.logFMMatrix(),    // Log(F/K) strikes
				volSurface.totVarMatrix(),   // Total variance
				sviSlices
			)
		};

		// ---------- Build the volatility surface object ----------

		// NOTE: Copy the input volatility surface
		core::VolSurface lvVolSurf{ volSurface };   
		std::vector<core::SliceData>& volSlices{ lvVolSurf.slices()};

		for (std::size_t i = 0; i < numTenors; ++i) volSlices[i].setWT(localtoVar[i]);

		return lvVolSurf;
	}

} // namespace uv::models::localvol

namespace uv::models::localvol::detail
{	
	Matrix<Real> localTotVar(const Vector<Real>& tenors,
		const Matrix<Real>& logFM,
		const Matrix<Real>& totVar,
		const Vector<svi::SVISlice>& sviSlices)
	{
		// ---------- Initialize data ----------

		const std::size_t numTenors{ tenors.size() };
		const std::size_t numStrikes{ totVar[0].size() };
		Matrix<Real> dwdT(numStrikes, Vector<Real>(numTenors));
		Matrix<Real> results(numTenors, Vector<Real>(numStrikes));

		// ---------- Calculate dw/dT matrix ----------

		const Matrix<Real> totVarT{ utils::transposeMatrix(totVar) };

		for (std::size_t i = 0; i < numStrikes; ++i)
		{
			dwdT[i] = math::pchipDerivatives(tenors, totVarT[i]);
		}

		// Transpose back into original dimensions
		dwdT = utils::transposeMatrix(dwdT);

		// ---------- Local total variance from SVI ----------

		for (std::size_t i = 0; i < numTenors; ++i)
		{
			// ---------- SVI data ----------

			const svi::SVISlice& sviSlice{ sviSlices[i] };
			const Real a{ sviSlice.a };
			const Real b{ sviSlice.b };
			const Real rho{ sviSlice.rho };
			const Real m{ sviSlice.m };
			const Real sigma{ sviSlice.sigma };
			const Real T{ sviSlice.T };

			// ---------- Gatheral-style local total variance ----------

			for (std::size_t j = 0; j < numStrikes; ++j)
			{
				const Real gk{ svi::gk(a, b, rho, m, sigma, logFM[i][j]) };
				results[i][j] = dwdT[i][j] * T / gk;
			}
		}

		return results;
	}

} // uv::models::localvol::detail