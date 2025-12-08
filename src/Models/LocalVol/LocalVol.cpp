// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Functions.cpp
 * Author:      Álvaro Sánchez de Carlos
 * Created:     2025-12-08
 *
 * Description:
 *   [Brief description of what this file declares or implements.]
 *
 * Copyright (c) 2025 Álvaro Sánchez de Carlos
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at:
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under this License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the LICENSE for the specific language governing permissions and
 * limitations under this License.
 */



#include "Models/LocalVol/Functions.hpp"
#include "Utils/Aux/Errors.hpp"
#include "Math/Interpolation.hpp"
#include "Math/Functions.hpp"
#include "Core/Functions.hpp"
#include "Models/SVI/Functions.hpp"
#include "Core/SliceData.hpp"
#include "Models/LocalVol/Pricer.hpp"

#include <cstddef>
#include <string>
#include <cmath>
#include <iostream>

using namespace uv;

namespace uv::models::localvol
{
	core::VolSurface buildSurface(const core::VolSurface& volSurface,
		const Vector<models::svi::Params>& paramsSVI)
	{

		// ---------- Check matching dimensions ----------
		
		const std::size_t numTenors{ volSurface.numTenors() };
		const std::size_t numStrikes{ volSurface.numStrikes() };
		const std::size_t numSlices{ paramsSVI.size() };

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
			const Real sviTenor{ paramsSVI[i].T };

			// Throw if surface and SVI tenors do not match
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
				volSurface.logKFMatrix(),    // Log(F/K) strikes
				volSurface.totVarMatrix(),   // Total variance
				paramsSVI
			)
		};

		// ---------- Build the volatility surface object ----------

		// NOTE: Copy the input volatility surface
		core::VolSurface lvVolSurf{ volSurface };   
		std::vector<core::SliceData>& volSlices{ lvVolSurf.slices()};

		for (std::size_t i = 0; i < numTenors; ++i) volSlices[i].setWT(localtoVar[i]);

		return lvVolSurf;
	}

	Vector<Real> price(const core::VolSurface& localVolSurface,
		const core::MarketData& marketData,
		const std::size_t NT,
		const std::size_t NF,
		const unsigned int X)
	{
		// ---------- Initialization ----------

		Vector<Real> surfaceTenors{ localVolSurface.tenors() };
		Vector<Real> sufaceStrikes{ localVolSurface.strikes()};

		Pricer pricer
		(
			surfaceTenors,
			sufaceStrikes,
			marketData,
			NT,
			NF,
			X
		);

		// ---------- Pricing ----------

		return pricer.price
		(
			localVolSurface.varMatrix(),
			surfaceTenors,
			sufaceStrikes
		);
	}

} // namespace uv::models::localvol

namespace uv::models::localvol::detail
{	
	Matrix<Real> localTotVar(const Vector<Real>& tenors,
		const Matrix<Real>& logKF,
		const Matrix<Real>& totVar,
		const Vector<svi::Params>& paramsSVI)
	{
		// ---------- Initialize data ----------

		const std::size_t numTenors{ tenors.size() };
		const std::size_t numStrikes{ totVar[0].size() };
		Matrix<Real> dwdT(numStrikes, Vector<Real>(numTenors));
		Matrix<Real> results(numTenors, Vector<Real>(numStrikes));

		// ---------- Calculate dw/dT matrix ----------

		const Matrix<Real> totVarT{ uv::core::transposeMatrix(totVar) };

		for (std::size_t i = 0; i < numStrikes; ++i)
		{
			dwdT[i] = math::interp::pchipDerivatives(tenors, totVarT[i]);
		}

		// Transpose back into original dimensions
		dwdT = uv::core::transposeMatrix(dwdT);

		// ---------- Local total variance from SVI ----------

		for (std::size_t i = 0; i < numTenors; ++i)
		{
			// ---------- SVI data ----------

			const svi::Params& sviSlice{ paramsSVI[i] };
			const Real a{ sviSlice.a };
			const Real b{ sviSlice.b };
			const Real rho{ sviSlice.rho };
			const Real m{ sviSlice.m };
			const Real sigma{ sviSlice.sigma };
			const Real T{ sviSlice.T };

			// ---------- Gatheral-style local total variance ----------

			for (std::size_t j = 0; j < numStrikes; ++j)
			{
				const Real gk{ svi::gk(a, b, rho, m, sigma, logKF[i][j]) };
				results[i][j] = dwdT[i][j] * T / gk;
			}
		}

		return results;
	}

} // uv::models::localvol::detail