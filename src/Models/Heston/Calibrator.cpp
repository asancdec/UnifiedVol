// SPDX-License-Identifier: Apache-2.0
/*
 * File:        HestonCalibrator.cpp
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


#include "Models/Heston/Calibrator.hpp"
#include "Utils/Aux/Errors.hpp"
#include "Utils/Types.hpp"

#include <array>
#include <cmath>
#include <cstddef>

namespace uv::models::heston::calibrator::detail
{	
    void validateInputs(const Vector<Real>& tenors,
        const Vector<Real>& strikes,
        const Vector<Real>& forwards,
        const Vector<Real>& rates,
        const core::Matrix<Real>& callM)
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

	std::array<double, 5> initGuess() noexcept
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

	std::array<double, 5> lowerBounds() noexcept
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

	std::array<double, 5> upperBounds() noexcept
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
} // namespace uv::models::heston_calibrator::detail

