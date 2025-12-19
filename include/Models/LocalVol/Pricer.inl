// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Pricer.inl
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


#include "Core/Matrix/Functions.hpp"
#include "Math/PDE/Functions.hpp"
#include "Utils/IO/Functions.hpp"

#include <cmath>
#include <iostream>

namespace uv::models::localvol
{
    template <std::floating_point T>
    Pricer<T>::Pricer(std::span<const T> r,
        std::span<const T> tGrid,
        std::span<const T> xGrid)
        :
        r_(r.begin(), r.end()),
        tGrid_(tGrid.begin(), tGrid.end()),
        xGrid_(xGrid.begin(), xGrid.end()),
        nt_(tGrid_.size()),
        nx_(xGrid_.size()),
        nxMid_(nx_ - 1),
        dt_(tGrid_[1] - tGrid_[0]),
        dx_(xGrid_[1] - xGrid_[0]),
        initGuess_(math::pde::fokkerPlankInit<T>(0, xGrid_)),
        xMidGrid_(math::pde::createXMidGrid<T>(xGrid_))
    {
        // TODO: Validate inputs
        // TODO: Set grid of strikes and tenors to price
    }

    template <std::floating_point T>
    T Pricer<T>::price(T localVar) const
    {
        // TODO: Assume non constant localVar
        // TODO: Return a price grid, not just a T

        core::Matrix<T> B(nt_, nxMid_, -0.5 * localVar);
        core::Matrix<T> C{ B };  // TODO:: derivative of vol

        core::Matrix<T> weight{ math::pde::changCooperWeights<T>(B, C, dx_) };

        weight.print();


        return T(3);
    }

} // namespace uv::models::localvol