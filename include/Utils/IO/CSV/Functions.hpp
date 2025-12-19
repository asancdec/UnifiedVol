// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Functions.hpp
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


#pragma once

#include "Core/VolSurface.hpp"

#include <string>
#include <concepts>

namespace uv::utils
{
	/**
	 * @brief Reads a CSV file containing a volatility surface and returns a VolSurface object.
	 *
	 * The CSV file is expected to have the following structure:
	 * - The first row contains moneyness (K/S).
	 * - The first column of each subsequent row contains tenors.
	 * - The remaining cells contain implied volatilities for each strike-maturity pair.
	 *
	 * @param filename Path to the CSV file.
	 * @return VolSurface Object containing the strikes, tenors, and implied volatilities.
	 */
	template <std::floating_point T>
	core::VolSurface<T> readVolSurface(const std::string& filename, const core::MarketData<T>& mktData);

} // namespace uv::utils

#include "Functions.inl"