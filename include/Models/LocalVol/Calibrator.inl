// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Calibrator.inl
 * Author:      Alvaro Sanchez de Carlos
 * Created:     2026-01-22
 *
 * Description:
 *   Local volatility surface calibration namespace
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

#include "Utils/Aux/Errors.hpp"
#include "Models/LocalVol/VarianceView.hpp"
#include "Math/Interpolation/Interpolator.hpp"
#include "Math/Interpolation/Policies.hpp"

#include <format>

namespace uv::models::localvol::calibrator
{

	template
	<
		std::floating_point T,
		std::size_t NT,
		std::size_t NX,
        class Interpolator
	>
	Surface<T> calibrate
	(
		const core::Matrix<T>& callM,
		const Vector<T>& tenors,
		const core::Matrix<T>& logKF,
        const core::Matrix<T>& totVar, 
        Pricer<T, NT, NX, Interpolator>& pricer
	)
	{
		// ---------- Validate ----------
		
        detail::validate
        (
            callM, 
            tenors,
            logKF,
            totVar
        );

		// ---------- Allocate ----------

        // Size
        std::size_t nTenors{ logKF.rows() };
		std::size_t nStrikes{ logKF.cols() };

        // Buffers
        Vector<T> localVar(nStrikes);
        Vector<T> dydx(nStrikes);

        // Parameters
        Vector<T> params;

        // Derivatives policy
        // NOTE: the architecture ensures the derivatives calculation policy
        // matches that of the the interpolator used in the Pricer instance
        typename Pricer<T, NT, NX, Interpolator>::derivatives_type derivPolicy{};

        // ---------- Initial guess ----------

        detail::coldGuess<T>
        (
            tenors.front(),
            totVar[0],
            localVar
        );

        // ---------- Extract ----------

        std::span<const T> logKFSlice{logKF[0]};
        T tenor{tenors[0]};

        // Calculate derivatives at knots
        derivPolicy
        (
            logKFSlice,
            localVar,
            dydx
        );

        // Pack into VarianceView struct
        VarianceView<T> localVarView
        {
            logKFSlice,
            localVar,
            dydx
        };

        // Price
        Vector<T> prices 
        {
            pricer.priceNormalized
            (
                tenor,
                logKFSlice,
                localVarView
            )
        };

        //for (std::size_t i = 0; i < results.size(); ++i)
        //{
        //    std::cout << results[i] * 485.77548 << ' ';
        //}

		return Surface<T>
		{
			tenors,
			logKF,
		    logKF
		};
	}

} // namespace uv::models::localvol

namespace uv::models::localvol::calibrator::detail
{
	template<std::floating_point T>
    void validate
    (
        const core::Matrix<T>& callM,
        const Vector<T>& tenors,
        const core::Matrix<T>& logKF,
        const core::Matrix<T>& totVar
    )
    {
        // ---------- Size ----------

        const std::size_t nRows{ callM.rows() };
        const std::size_t nCols{ callM.cols() };

        UV_REQUIRE
        (
            nRows > 0 && nCols > 0,
            ErrorCode::InvalidArgument,
            std::format
            (
                "validate: callM must be non-empty, got {}x{}",
                nRows, nCols
            )
        );

        UV_REQUIRE
        (
            nRows == logKF.rows() && nCols == logKF.cols(),
            ErrorCode::InvalidArgument,
            std::format
            (
                "validate: callM/logKF size mismatch — callM is {}x{}, logKF is {}x{}",
                nRows, nCols, logKF.rows(), logKF.cols()
            )
        );

        UV_REQUIRE
        (
            nRows == totVar.rows() && nCols == totVar.cols(),
            ErrorCode::InvalidArgument,
            std::format
            (
                "validate: callM/totVar size mismatch — callM is {}x{}, totVar is {}x{}",
                nRows, nCols, totVar.rows(), totVar.cols()
            )
        );

        UV_REQUIRE
        (
            tenors.size() == nRows,
            ErrorCode::InvalidArgument,
            std::format
            (
                "validate: row mismatch — tenors={}, rows={}",
                tenors.size(), nRows
            )
        );

        // ---------- Values ----------

        // Tenors
        for (std::size_t i = 0; i < nRows; ++i)
        {
            UV_REQUIRE
            (
                std::isfinite(tenors[i]),
                ErrorCode::InvalidArgument,
                std::format
                (
                    "validate: tenors[{}]={} is not finite",
                    i, tenors[i]
                )
            );

            UV_REQUIRE
            (
                tenors[i] > 0.0,
                ErrorCode::InvalidArgument,
                std::format
                (
                    "validate: tenors[{}]={} must be > 0",
                    i, tenors[i]
                )
            );

            if (i > 0)
            {
                UV_REQUIRE
                (
                    tenors[i] > tenors[i - 1],
                    ErrorCode::InvalidArgument,
                    std::format
                    (
                        "validate: tenors must be strictly increasing, but tenors[{}]={} <= tenors[{}]={}",
                        i, tenors[i], i - 1, tenors[i - 1]
                    )
                );
            }
        }

        // logKF, callM, totVar matrices
        for (std::size_t t = 0; t < nRows; ++t)
        {
            // ---------- Values ----------

            for (std::size_t k = 0; k < nCols; ++k)
            {
                const T c{ callM[t][k] };
                const T xk{ logKF[t][k] };
                const T w{ totVar[t][k] };

                UV_REQUIRE
                (
                    std::isfinite(c),
                    ErrorCode::InvalidArgument,
                    std::format
                    (
                        "validate: callM({}, {})={} is not finite",
                        t, k, c
                    )
                );

                UV_REQUIRE
                (
                    c >= 0.0,
                    ErrorCode::InvalidArgument,
                    std::format
                    (
                        "validate: callM({}, {})={} must be >= 0",
                        t, k, c
                    )
                );

                UV_REQUIRE
                (
                    std::isfinite(xk),
                    ErrorCode::InvalidArgument,
                    std::format
                    (
                        "validate: logKF({}, {})={} is not finite",
                        t, k, xk
                    )
                );

                UV_REQUIRE
                (
                    std::isfinite(w),
                    ErrorCode::InvalidArgument,
                    std::format
                    (
                        "validate: totVar({}, {})={} is not finite",
                        t, k, w
                    )
                );

                UV_REQUIRE
                (
                    w >= 0.0,
                    ErrorCode::InvalidArgument,
                    std::format
                    (
                        "validate: totVar({}, {})={} must be >= 0",
                        t, k, w
                    )
                );
            }

            // ---------- Monotonicity ----------

            for (std::size_t k = 1; k < nCols; ++k)
            {
                const T prev{ logKF[t][k - 1] };
                const T curr{ logKF[t][k] };

                UV_REQUIRE
                (
                    curr > prev,
                    ErrorCode::InvalidArgument,
                    std::format
                    (
                        "validate: logKF row {} must be strictly increasing, \
                         but logKF({}, {})={} <= logKF({}, {})={}",
                        t, t, k, curr, t, k - 1, prev
                    )
                );
            }
        }

        // ---------- No-arbitrage ----------

         // Total variance
         for (std::size_t k = 0; k < nCols; ++k)
         {
             for (std::size_t t = 1; t < nRows; ++t)
             {
                 const T prev{ totVar[t - 1][k] };
                 const T curr{ totVar[t][k] };
        
                 UV_REQUIRE
                 (
                     curr >= prev,
                     ErrorCode::InvalidArgument,
                     std::format
                     (
                         "validate: totVar must be non-decreasing in tenor at strike col {}, \
                         but totVar({}, {})={} < totVar({}, {})={}",
                         k, t, k, curr, t - 1, k, prev
                     )
                 );
             }
         }
    }
    
    template<std::floating_point T>
    void coldGuess
    (
        T tenor,
        std::span<const T> totVar,
        std::span<T> localVar
    ) noexcept
    {
        T invTenor{1.0/ tenor};

        for (std::size_t i{ 0 }; i < localVar.size(); ++i)
        {
            localVar[i] = totVar[i] * invTenor;
        }
    }




} // namespace uv::models::localvol::calibrator::detail