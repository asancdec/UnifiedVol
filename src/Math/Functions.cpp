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

#include <iostream>

namespace uv::math
{
    Vector<Real> blackScholes(Real t,
        Real r,
        Real q,
        const Vector<Real>& vol,
        Real S,
        const Vector<Real>& K,
        bool isCall) noexcept
    {

        // ---------- Validate inputs ----------

        const std::size_t N{ vol.size() };

        UV_REQUIRE(
            K.size() == N,
            ErrorCode::InvalidArgument,
            "blackScholes: vol / strike size mismatch"
        );

        // ---------- Calculate Black-Scholes prices ----------

        Vector<Real> out(N);

        for (std::size_t i = 0; i < N; ++i)
        {
            out[i] = uv::math::blackScholes(
                t,
                r,
                q,
                vol[i],
                S,
                K[i],
                isCall
            );
        }

        return out;
    }

    Matrix<Real> blackScholes(const Vector<Real>& t,
        const Vector<Real>& r,
        const Vector<Real>& q,
        const Matrix<Real>& vol,
        Real S,
        const Vector<Real>& K,
        bool isCall
    ) noexcept
    {
        std::cout
            << "Nt=" << t.size()
            << " K=" << K.size()
            << " vol0=" << (vol.empty() ? 0ULL : vol.front().size())
            << " volLast=" << (vol.empty() ? 0ULL : vol.back().size())
            << '\n';

        // ---------- Validate inputs ----------

        const std::size_t Nt{ t.size() };

        UV_REQUIRE(
            r.size() == Nt &&
            q.size() == Nt &&
            vol.size() == Nt,
            ErrorCode::InvalidArgument,
            "blackScholes(matrix): input size mismatch"
        );

        // ---------- Calculate Black-Scholes prices ----------

        Matrix<Real> out(Nt);

        for (std::size_t i = 0; i < Nt; ++i)
        {
            out[i] = uv::math::blackScholes(
                t[i],
                r[i],
                q[i],
                vol[i],   // Vector<T>
                S,
                K,
                isCall
            );
        }

        return out;
    }

} // nampespace uv::math
