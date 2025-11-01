/**
* HestonCalibrator.cpp
* Author: Alvaro Sanchez de Carlos
*/

#include "Models/Heston/HestonCalibrator.hpp"
#include <ceres/ceres.h>
#include <array>
#include <format>
#include <sstream>
#include <iomanip> 

namespace uv
{

	struct HestonCalibrator::PriceResidualND
	{	
		// Member variables
		const double T_, F_, r_, K_, callPriceMkt_;
		::std::shared_ptr<const HestonPricer> pricer_;

		// Constructor
		explicit PriceResidualND(double T, 
			double F, 
			double r, 
			double K,
			double callPriceMkt,
			::std::shared_ptr<const HestonPricer> pricer) noexcept:
		T_(T), F_(F), r_(r), K_(K),	callPriceMkt_(callPriceMkt), pricer_(std::move(pricer)) {}

		// Overload () operator with Ceres functor signature
		bool operator()(const double* const p, double* residual) const noexcept
		{
			const double callPriceModel
			{ 
				pricer_->callPrice
				(
					p[0],	// kappa
					p[1],	// theta
					p[2],	// sigma
					p[3],	// rho
					p[4],   // v0
					T_,		// T
					F_,		// F
					r_,		// r
					K_		// K
				)
			};

			residual[0] = callPriceModel - callPriceMkt_;
			return true;
		}
	};

	void HestonCalibrator::calibrate(const VolSurface& mktVolSurf,
		::std::shared_ptr<const HestonPricer> pricer,
		CalibratorCeres<5>& calibrator)
	{
		// Copy market volatility surface
		VolSurface hestonVolSurf{ mktVolSurf };

		// Set initial guess, lower and upper bounds
		calibrator.setGuessBounds
		(
			::std::move(initGuess()),
			::std::move(lowerBounds()),
			::std::move(upperBounds())
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

			// Iterate through volatilities vector
			for (std::size_t i = 0; i < K.size(); ++i)
			{
				// Build the cost function using central finite difference scheme
				auto* cost = new ::ceres::NumericDiffCostFunction<PriceResidualND, ::ceres::CENTRAL, 1, 5>
					(
						new PriceResidualND
						(
							T,
							F,
							r,
							K[i],
							callPriceMkt[i],
							pricer
						)
					);

				// Huber Loss function against outliers: experiment
				::ceres::LossFunction* loss = new ::ceres::HuberLoss(1.0);

				calibrator.problem_.AddResidualBlock(cost, loss, calibrator.x_.data());
			}
		}

		::ceres::Solver::Options options;
		options.trust_region_strategy_type = ::ceres::LEVENBERG_MARQUARDT;
		options.linear_solver_type = ::ceres::DENSE_QR;
		options.max_num_iterations = 10;           // tune
		options.function_tolerance = 1e-12;        // tune
		options.parameter_tolerance = 1e-12;       // tune
		options.gradient_tolerance = 0.0;          // usually leave 0 with LM
		options.num_threads = std::max(1u, std::thread::hardware_concurrency());

		::ceres::Solver::Summary summary;
		::ceres::Solve(options, &calibrator.problem_, &summary);


		// Print full report
		UV_INFO(summary.FullReport());

		// --- Print final parameters (manual join) ---
		std::ostringstream oss;
		oss << ::std::fixed << std::setprecision(6);  // optional formatting

		for (::std::size_t i = 0; i < calibrator.x_.size(); ++i)
		{
			oss << calibrator.x_[i];
			if (i + 1 < calibrator.x_.size())
				oss << ", ";
		}

		UV_INFO(::std::format("Final params: [{}]", oss.str()));

	}

	::std::array<double, 5> HestonCalibrator::initGuess() noexcept
	{
		return
		{
			2.5,     // kappa
			0.09,    // theta
			0.60,    // sigma
			-0.50,   // rho
			0.09     // vo
		};
	}

	::std::array<double, 5> HestonCalibrator::lowerBounds() noexcept
	{
		return
		{
			0.001,   // kappa
			0.001,   // theta
			0.001,   // sigma
			-0.999,  // rho
			0.001    // vo
		};

	}

	::std::array<double, 5> HestonCalibrator::upperBounds() noexcept
	{
		return
		{
			10.0,   // kappa
			0.5,    // theta
			10.0,   // sigma
			0.999,  // rho
			0.5     // vo
		};
	}
}

