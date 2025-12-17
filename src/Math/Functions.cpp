// SPDX-License-Identifier: Apache-2.0
/*
 * File:        CSVRead.hpp
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

#include "Math/Functions.hpp"
#include "Utils/Aux/Errors.hpp"
#include "Utils/Types.hpp"

#include <cstddef>
#include <span>

namespace uv::math
{
    core::Matrix<Real> blackScholes(const Vector<Real>& t,
        const Vector<Real>& r,
        const Vector<Real>& q,
        const core::Matrix<Real>& vol,
        Real S,
        const Vector<Real>& K,
        bool isCall
    )
    {
        // ---------- Validate inputs ----------

        const std::size_t Nt{ t.size() };
        const std::size_t Nk{ K.size() };

        UV_REQUIRE(
            r.size() == Nt &&
            q.size() == Nt &&
            vol.rows() == Nt &&
            Nk == vol.cols(),
            ErrorCode::InvalidArgument,
            "blackScholes(matrix): input size mismatch"
        );

        // ---------- Calculate Black-Scholes prices ----------

        core::Matrix<Real> out(Nt, Nk);

        for (std::size_t i = 0; i < Nt; ++i)
        {
            std::span<Real> outRow{ out[i] };
            std::span<const Real> volRow{ vol[i] };
            const Real ti{ t[i] };
            const Real ri{ r[i] };
            const Real qi{ q[i] };


            for (std::size_t j = 0; j < Nk; ++j)
            {
                outRow[j] = math::blackScholes(
                    ti, ri, qi, volRow[j],
                    S, K[j], isCall
                );
            }
        }

        return out;
    }

} // nampespace uv::math