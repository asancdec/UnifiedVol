// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Calibrator.inl
 * Author:      Alvaro Sanchez de Carlos
 * Created:     2025-12-08
 *
 * Description:
 *   Heston calibration routines.
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
	template <std::floating_point T, std::size_t N, typename Policy>
	Params<T> calibrate
	(
		const Vector<T>& tenors,
		const Vector<T>& strikes,
		const Vector<T>& forwards,
		const Vector<T>& rates,
		const core::Matrix<T>& callM,
		Pricer<T, N>& pricer,
		opt::ceres::Optimizer<5, Policy>& optimizer,
		const opt::WeightATM<double>& weightATM
	)
	{
		// ---------- Validate  ----------

		detail::validateInputs<T>(tenors, strikes, forwards, rates, callM);

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

		// Allocate
		Vector<double> bufferWeights(numStrikes);
		Vector<double> logKF(numStrikes);

		for (std::size_t i = 0; i < numTenors; ++i)
		{
			// Extract data

			std::span<const double> callMDRow{callMD[i]};
			const double t{ tenorsD[i] };
			const double F{ forwardsD[i]};
			const double r{ ratesD[i] };

			// logKF strikes
			for (std::size_t j{ 0 }; j < numStrikes; ++j)
				logKF[j] = std::log(strikesD[j] / F);

			// ATM weights
			opt::weightsATM<double>
			(
				logKF,
				weightATM,
				bufferWeights
			);

			for (std::size_t j = 0; j < numStrikes; ++j)
			{
				optimizer.addAnalyticResidual
				(
					std::make_unique<detail::PriceResidualJac<T, N>>
					(
						t, 
						F, 
						r, 
						strikesD[j], 
						callMDRow[j], 
						bufferWeights[j],
						pricer
					)
				);
			}
		}

		// ---------- Run calibration ----------

		const std::array<double, 5> params{ optimizer.optimize() };
		return Params<T>
		{
			T(params[0]),  // kappa
			T(params[1]),  // theta
			T(params[2]),  // sigma
			T(params[3]),  // rho
			T(params[4])   // v0 
		};
	}

	template <std::floating_point T, std::size_t N>
	core::VolSurface<T> buildSurface
	(
		const core::VolSurface<T>& volSurface,
		const Pricer<T, N>& pricer
	)
	{
		// ---------- Extract ----------
		
		const std::size_t numTenors{ volSurface.numTenors() };
		const std::size_t numStrikes{ volSurface.numStrikes() };
		const Vector<T>& tenors{ volSurface.tenors() };
		const Vector<T>& strikes{ volSurface.strikes() };
		const Vector<T>& forwards{ volSurface.forwards() };
		const Vector<T>& rates{ volSurface.rates() };

		// ---------- Copy and set surface ----------

		core::VolSurface<T> hestonVolSurface{ volSurface };
		hestonVolSurface.setCallPrices
		(
			core::generateIndexed<T>
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

} // namespace uv::models::heston::calibrator


namespace uv::models::heston::calibrator::detail
{   
	template <std::floating_point T>
	void validateInputs
	(
		const Vector<T>& tenors,
		const Vector<T>& strikes,
		const Vector<T>& forwards,
		const Vector<T>& rates,
		const core::Matrix<T>& callM
	)
	{
		UV_REQUIRE(!tenors.empty(),
			ErrorCode::InvalidArgument,
			"validateInputs: tenors is empty");

		UV_REQUIRE(!strikes.empty(),
			ErrorCode::InvalidArgument,
			"validateInputs: strikes is empty");

		UV_REQUIRE(!forwards.empty(),
			ErrorCode::InvalidArgument,
			"validateInputs: forwards is empty");

		UV_REQUIRE(!rates.empty(),
			ErrorCode::InvalidArgument,
			"validateInputs: rates is empty");

		UV_REQUIRE(!callM.empty(),
			ErrorCode::InvalidArgument,
			"validateInputs: callM is empty");

		const std::size_t numTenors{ tenors.size() };
		const std::size_t numStrikes{ strikes.size() };

		UV_REQUIRE(forwards.size() == numTenors,
			ErrorCode::InvalidArgument,
			"validateInputs: forwards size must equal tenors size");

		UV_REQUIRE(rates.size() == numTenors,
			ErrorCode::InvalidArgument,
			"validateInputs: rates size must equal tenors size");

		UV_REQUIRE(callM.rows() == numTenors,
			ErrorCode::InvalidArgument,
			"validateInputs: callM rows must equal number of tenors");

		UV_REQUIRE(callM.cols() == numStrikes,
			ErrorCode::InvalidArgument,
			"validateInputs: callM columns must equal strikes size");

		// Tenors strictly increasing + positive
		for (std::size_t i = 0; i < numTenors; ++i)
		{
			UV_REQUIRE(tenors[i] > Real(0),
				ErrorCode::InvalidArgument,
				"validateInputs: tenors must be > 0");

			if (i > 0)
			{
				UV_REQUIRE(tenors[i] > tenors[i - 1],
					ErrorCode::InvalidArgument,
					"validateInputs: tenors must be strictly increasing");
			}
		}

		// Strikes strictly increasing + positive
		for (std::size_t j = 0; j < numStrikes; ++j)
		{
			UV_REQUIRE(strikes[j] > Real(0),
				ErrorCode::InvalidArgument,
				"validateInputs: strikes must be > 0");

			if (j > 0)
			{
				UV_REQUIRE(strikes[j] > strikes[j - 1],
					ErrorCode::InvalidArgument,
					"validateInputs: strikes must be strictly increasing");
			}
		}

		// Forwards positive
		for (std::size_t i = 0; i < numTenors; ++i)
		{
			UV_REQUIRE(forwards[i] > Real(0),
				ErrorCode::InvalidArgument,
				"validateInputs: forwards must be > 0");
		}
	}

	template <std::floating_point T, std::size_t N>
	struct PriceResidualJac final
		: public ceres::SizedCostFunction<1, 5> 
	{
		const double T_, F_, r_, K_, callPriceMkt_, w_;
		const Pricer<T, N>* pricer_; 

		PriceResidualJac
		(
			double t, 
			double F, 
			double r, 
			double K, 
			double callPriceMkt,
			double w,
			const Pricer<T, N>& pricer
		) noexcept
			: 
			T_(t), 
			F_(F), 
			r_(r),
			K_(K), 
			w_(w),
			callPriceMkt_(callPriceMkt),
			pricer_(&pricer) 
		{}

		bool Evaluate
		(
			double const* const* parameters,
			double* residuals,
			double** jacobians
		) const override
		{
			const double* p = parameters[0];
			const auto pg = pricer_->callPriceWithGradient
			(
				p[0], 
				p[1], 
				p[2], 
				p[3], 
				p[4], 
				T_, 
				F_, 
				r_,
				K_
			);

			residuals[0] = double((pg[0] - callPriceMkt_) * w_) ;

			if ((jacobians != nullptr) && (jacobians[0] != nullptr)) 
			{
				double* J = jacobians[0];
				J[0] = double(pg[1] * w_);
				J[1] = double(pg[2] * w_);
				J[2] = double(pg[3] * w_);
				J[3] = double(pg[4] * w_);
				J[4] = double(pg[5] * w_);
			}
			return true;
		}
	};
} // namespace uv::models::heston::calibrator::detail