/**
* HestonCalibrator.inl
* Author: Alvaro Sanchez de Carlos
*/

#include <ceres/ceres.h>

#include <memory>
#include <algorithm>

namespace uv
{
	template <::std::size_t N>
	struct HestonCalibrator::PriceResidualND
	{
		// Member variables
		const double T_, F_, r_, K_, callPriceMkt_;
		const HestonPricer<N>* pricer_;

		// Constructor
		explicit PriceResidualND(double T,
			double F,
			double r,
			double K,
			double callPriceMkt,
			const HestonPricer<N>& pricer) noexcept :
			T_(T), F_(F), r_(r), K_(K), callPriceMkt_(callPriceMkt), pricer_(&pricer) {}

		// Overload () operator with Ceres functor signature
		bool operator()(const double* const p, double* residual) const noexcept
		{
			const double callPriceModel = pricer_->callPrice
			(
				p[0], p[1], p[2], p[3], p[4],   // kappa, theta, sigma, rho, v0
				T_, F_, r_, K_
			);

			// Define the residual function
			residual[0] = callPriceModel - callPriceMkt_;
			return true;
		}
	};

	template <::std::size_t N, typename Policy>
	VolSurface HestonCalibrator::calibrate(const VolSurface& mktVolSurf,
		HestonPricer<N>& pricer,
		CalibratorCeres<5, Policy>& calibrator)
	{
		// Copy market volatility surface
		VolSurface hestonVolSurf{ mktVolSurf };

		// Set initial guess, lower and upper bounds
		calibrator.setGuessBounds
		(
			HestonCalibrator::initGuess(),
			HestonCalibrator::lowerBounds(),
			HestonCalibrator::upperBounds()
		);

		// Add residual blocks
		for (const auto& slice : hestonVolSurf.slices())
		{
			// Extract slice parameters
			const double T{ slice.T() };
			const double F{ slice.F() };
			const double r{ slice.r() };
			const ::std::vector<double>& K{ slice.K() };
			const ::std::vector<double>& callPriceMkt(slice.callBS());

			// Add numerical residue blocks
			for (std::size_t i = 0; i < K.size(); ++i)
			{	
				calibrator.template addNumericResidual<PriceResidualND<N>>
					(PriceResidualND<N>(T, F, r, K[i], callPriceMkt[i], pricer));
			}
		}

		// Calibrate Heston Pricer parameters
		pricer.setHestonParams(calibrator.optimize());

		// Fill in the calibrated volatility surface
		for (auto& slice : hestonVolSurf.slices())
		{
			// Extract slice parameters
			const double T{ slice.T() };
			const double F{ slice.F() };
			const double r{ slice.r() };
			const ::std::vector<double>& K{ slice.K() };
			std::vector<double> modelCall(K.size());

			// Fill modelCall with call prices for each strike
			std::transform
			(
				K.begin(),
				K.end(),
				modelCall.begin(),
				[T, F, r, &pricer](double Ki) noexcept { return pricer.callPrice(T, F, r, Ki); }
			);

			// Store the calculated prices in the volatility surface object
			slice.setCallBS(std::move(modelCall));
		}

		return hestonVolSurf;  // RVO C++20
	}
}