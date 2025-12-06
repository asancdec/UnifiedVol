/**
* Calibrator.inl
* Author: Alvaro Sanchez de Carlos
*/

#include <ceres/ceres.h>

#include <memory>
#include <algorithm>

namespace uv::models::heston::calibrator
{
	template <std::size_t N, typename Policy>
	core::VolSurface calibrate(const core::VolSurface& mktVolSurf,
		Pricer<N>& pricer,
		cal::ceres::Calibrator<5, Policy>& calibrator)
	{
		// Copy market volatility surface
		core::VolSurface hestonVolSurf{ mktVolSurf };

		// Set initial guess, lower and upper bounds
		calibrator.setGuessBounds
		(
			detail::initGuess(),
			detail::lowerBounds(),
			detail::upperBounds()
		);

		// Add residual blocks
		for (const auto& slice : hestonVolSurf.slices())
		{
			// Extract slice parameters
			const double T{ double(slice.T()) };
			const double F{ double(slice.F()) };
			const double r{ double(slice.r() )};
			Vector<double> Kcopy{ slice.K().cbegin(), slice.K().cend() };
			const Vector<double>& K{ Kcopy };
			Vector<double> callPriceMktCopy{ slice.callBS().cbegin(), slice.callBS().cend() };
			const Vector<double>& callPriceMkt{ callPriceMktCopy };

			for (std::size_t i = 0; i < K.size(); ++i)
			{
				calibrator.addAnalyticResidual
				(
					std::make_unique<detail::PriceResidualJac<N>>
					(
						T, F, r, K[i], callPriceMkt[i], pricer
					)
				);
			}
		}

		// Calibrate Pricer parameters
		pricer.setParams(calibrator.optimize());

		// Fill in the calibrated volatility surface
		for (auto& slice : hestonVolSurf.slices())
		{
			// Extract slice parameters
			const Real T{ slice.T() };
			const Real F{ slice.F() };
			const Real r{ slice.r() };
			const Vector<Real>& K{ slice.K() };
			Vector<Real> modelCall(K.size());

			// Fill modelCall with call prices for each strike
			std::transform
			(
				K.begin(),
				K.end(),
				modelCall.begin(),
				[T, F, r, &pricer](Real Ki) noexcept { return pricer.callPrice(T, F, r, Ki); }
			);

			// Store the calculated prices in the volatility surface object
			slice.setCallBS(std::move(modelCall));
		}
		return hestonVolSurf;  // RVO C++20
	}
} // namespace uv::models::heston_calibrator

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
			// returns [P, dP/dκ, dP/dθ, dP/dσ, dP/dρ, dP/dv0]
			const double* p = parameters[0];
			const auto pg = pricer_->callPriceWithGradient
			(
				p[0], p[1], p[2], p[3], p[4], T_, F_, r_, K_
			);

			residuals[0] = double(pg[0] - callPriceMkt_);

			if (jacobians && jacobians[0]) 
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