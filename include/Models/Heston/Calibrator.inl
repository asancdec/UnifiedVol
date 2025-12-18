// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Calibrator.inl
 * Author:      Alvaro Sanchez de Carlos
 * Created:     2025-12-08
 *
 * Description:
 *   [Brief description of what this file declares or implements.]
 *
 * Copyright (c) 2025 Alvaro Sanchez de Carlos
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

#include "Core/Functions.hpp"
#include "Core/Matrix/Functions.hpp"

#include <ceres/ceres.h>

#include <memory>
#include <algorithm>
#include <span>

namespace uv::models::heston::calibrator
{
	template <std::size_t N, typename Policy>
	Params calibrate(const Vector<Real>& tenors,
		const Vector<Real>& strikes,
		const Vector<Real>& forwards,
		const Vector<Real>& rates,
		const core::Matrix<Real>& callM,
		Pricer<N>& pricer,
		opt::ceres::Optimizer<5, Policy>& optimizer)
	{
		// ---------- Validate input data  ----------

		detail::validateInputs(tenors, strikes, forwards, rates, callM);

		// ---------- Convert data ----------

		// Ceres API accepts doubles only
		const Vector<double> tenorsD{ core::convertVector<double>(tenors) };
		const Vector<double> strikesD{ core::convertVector<double>(strikes) };
		const Vector<double> forwardsD{ core::convertVector<double>(forwards) };
		const Vector<double> ratesD{ core::convertVector<double>(rates) };
		const core::Matrix<double> callMD{ core::convertMatrix<double>(callM) };

		const std::size_t numTenors{ tenorsD.size() };
		const std::size_t numStrikes{ strikesD.size() };

		// ---------- Set initial guess and bounds ----------

		optimizer.setGuessBounds
		(
			detail::initGuess(),
			detail::lowerBounds(),
			detail::upperBounds()
		);

		// ---------- Add residual blocks ----------

		for (std::size_t i = 0; i < numTenors; ++i)
		{
			// Extract data
			std::span<const double> callMDRow{callMD[i]};
			const double T{ tenorsD[i] };
			const double F{ forwardsD[i]};
			const double r{ ratesD[i] };

			for (std::size_t j = 0; j < numStrikes; ++j)
			{
				optimizer.addAnalyticResidual
				(
					std::make_unique<detail::PriceResidualJac<N>>
					(
						T, F, r, strikesD[j], callMDRow[j], pricer
					)
				);
			}
		}

		// ---------- Run calibration ----------

		const std::array<double, 5> params{ optimizer.optimize() };
		return Params
		{
			Real(params[0]),  // kappa
			Real(params[1]),  // theta
			Real(params[2]),  // sigma
			Real(params[3]),  // rho
			Real(params[4])   // v0 
		};
	}

	template <std::size_t N>
	core::VolSurface buildSurface(const core::VolSurface& volSurface,
		const Pricer<N>& pricer)
	{
		// ---------- Extract data ----------
		
		const std::size_t numTenors{ volSurface.numTenors() };
		const std::size_t numStrikes{ volSurface.numStrikes() };
		const Vector<Real>& tenors{ volSurface.tenors() };
		const Vector<Real>& strikes{ volSurface.strikes() };
		const Vector<Real>& forwards{ volSurface.forwards() };
		const Vector<Real>& rates{ volSurface.rates() };

		// ---------- Copy and set surface ----------

		core::VolSurface hestonVolSurface{ volSurface };
		hestonVolSurface.setCallPrices
		(
			core::applyIndexed<Real>
			(
				numTenors,
				numStrikes,
				[&](std::size_t i, std::size_t j)
				{
					return pricer.callPrice
					(
						tenors[i],
						forwards[i],
						rates[i],
						strikes[j]
					);
				}
			)
		);

		return hestonVolSurface;
	}

} // namespace uv::models::heston::calibrator:


namespace uv::models::heston::calibrator::detail
{   
	template <std::size_t N>
	struct PriceResidualJac final
		: public ceres::SizedCostFunction<1, 5> 
	{
		const double T_, F_, r_, K_, callPriceMkt_;
		const Pricer<N>* pricer_; 

		PriceResidualJac(double T, double F, double r, double K, double callPriceMkt,
			const Pricer<N>& pricer) noexcept
			: T_(T), F_(F), r_(r), K_(K), callPriceMkt_(callPriceMkt), pricer_(&pricer) {}

		bool Evaluate(double const* const* parameters,
			double* residuals,
			double** jacobians) const override
		{
			const double* p = parameters[0];
			const auto pg = pricer_->callPriceWithGradient
			(
				p[0], p[1], p[2], p[3], p[4], T_, F_, r_, K_
			);

			residuals[0] = double(pg[0] - callPriceMkt_);

			if ((jacobians != nullptr) && (jacobians[0] != nullptr)) 
			{
				double* J = jacobians[0];
				J[0] = double(pg[1]);
				J[1] = double(pg[2]);
				J[2] = double(pg[3]);
				J[3] = double(pg[4]);
				J[4] = double(pg[5]);
			}
			return true;
		}
	};
} // namespace uv::models::heston::calibrator::detail