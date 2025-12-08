// SPDX-License-Identifier: Apache-2.0
/*
 * File:        VolSurface.cpp
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


#include "Core/VolSurface.hpp"
#include "Utils/Aux/Errors.hpp"
#include "Utils/IO/Log.hpp"

#include <iomanip>
#include <sstream>
#include <algorithm>
#include <string>

namespace uv::core
{
    VolSurface::VolSurface(const Vector<Real>& mny,
        const Matrix<Real>& vols,
        const Vector<Real>& tenors,
        const MarketData& mktData)
        : tenors_(tenors),
        numTenors_(tenors.size())
    {   
        // ---------- Sanity checks ----------
        
        // Number of volatility slices (one per maturity)
        const std::size_t dim0{ vols.size()};

        // Validate that the volatility matrix provides one slice per maturity
        UV_REQUIRE(
            dim0 == numTenors_,
            ErrorCode::InvalidArgument,
            "VolSurface: number of vol slices (" + std::to_string(dim0) +
            ") does not match number of tenors (" + std::to_string(numTenors_) + ")"
        );

        // Number of strikes in the strike grid (must be consistent across slices)
        const std::size_t dim1{ mny.size() };

        // Validate that each volatility slice has the same number of strikes
        for (std::size_t i = 0; i < dim0; ++i)
        {
            const std::size_t t{ vols[i].size() };
            UV_REQUIRE(
                t == dim1,
                ErrorCode::InvalidArgument,
                "VolSurface: inconsistent strike dimension — expected " +
                std::to_string(t) + " strikes, but slice " + std::to_string(i) +
                " has " + std::to_string(dim1)
            );
        }

        // ---------- Initialize member variables ----------

        slices_.reserve(numTenors_);
        for (size_t i = 0; i < numTenors_; ++i)
        {
            slices_.emplace_back
            (
                SliceData(tenors_[i],
                    mny,
                    vols[i],
                    mktData)
            );
        }

        numStrikes_ = dim1;
        strikes_.resize(numStrikes_);
        strikes_ = slices_[0].K();
    }

    void VolSurface::printVol() const noexcept
    {
        std::ostringstream oss;
        oss << '\n';
        oss << "T\\%S\t";

        // Header row (moneyness)
        for (const auto& m : slices_[0].mny())
            oss << std::fixed << std::setprecision(2) << m << '\t';
        oss << '\n';

        // Each tenor row
        for (size_t i = 0; i < numTenors_; ++i)
        {
            oss << std::fixed << std::setprecision(2) << tenors_[i] << '\t';

            for (const auto& v : slices_[i].vol())
                oss << std::fixed << std::setprecision(5) << v << '\t';

            oss << '\n';
        }

        UV_INFO(oss.str());
    }

    void VolSurface::printTotVar() const noexcept
    {
        std::ostringstream oss;
        oss << '\n';
        oss << "T\\k\t";

        // Header row (moneyness)
        for (const auto& m : slices_[0].mny())
            oss << std::fixed << std::setprecision(2) << m << '\t';
        oss << '\n';

        // Each tenor row
        for (size_t i = 0; i < numTenors_; ++i)
        {
            oss << std::fixed << std::setprecision(2) << tenors_[i] << '\t';

            for (const auto& v : slices_[i].wT())

                oss << std::fixed << std::setprecision(5) << v << '\t';

            oss << '\n';
        }

        UV_INFO(oss.str());
    }

    void VolSurface::printBSCall() const noexcept
    {
        std::ostringstream oss;
        oss << '\n';
        oss << "T\\k\t";

        // Header row (moneyness)
        for (const auto& m : slices_[0].mny())
            oss << std::fixed << std::setprecision(2) << m << '\t';
        oss << '\n';

        // Each tenor row
        for (size_t i = 0; i < numTenors_; ++i)
        {
            oss << std::fixed << std::setprecision(2) << tenors_[i] << '\t';

            for (const auto& v : slices_[i].callBS())
                oss << std::fixed << std::setprecision(7) << v << '\t';

            oss << '\n';
        }
        UV_INFO(oss.str());
    }

    Matrix<Real> VolSurface::totVarMatrix() const noexcept
    {
        // Allocate matrix:
        // outer dimension = tenors (rows)
        // inner dimension = strikes (columns)        
        Matrix<Real> matrix(
            numTenors_, Vector<Real>(numStrikes_)
        );

        // Populate total variance slices
        for (std::size_t i = 0; i < numTenors_; ++i)
        {
            // Extract total variance slice at tenor i
            matrix[i] = slices_[i].wT();
        }

        return matrix;
    }

    Matrix<Real> VolSurface::volMatrix() const noexcept
    {
        // Allocate matrix:
        // outer dimension = tenors (rows)
        // inner dimension = strikes (columns)        
        Matrix<Real> matrix(
            numTenors_, Vector<Real>(numStrikes_)
        );

        // Populate total variance slices
        for (std::size_t i = 0; i < numTenors_; ++i)
        {
            // Extract total variance slice at tenor i
            matrix[i] = slices_[i].vol();
        }

        return matrix;
    }

    Matrix<Real> VolSurface::logFMMatrix() const noexcept
    {
        // Allocate matrix:
        // outer dimension = tenors (rows)
        // inner dimension = strikes (columns)        
        Matrix<Real> matrix(
            numTenors_, Vector<Real>(numStrikes_)
        );

        // Populate total variance slices
        for (std::size_t i = 0; i < numTenors_; ++i)
        {
            // Extract total variance slice at tenor i
            matrix[i] = slices_[i].logFM();
        }

        return matrix;
    }

    std::vector<SliceData>& VolSurface::slices() noexcept
    {
        return slices_;
    }

    const Vector<Real>& VolSurface::tenors() const noexcept
    {
        return tenors_;
    }

    const Vector<Real>& VolSurface::strikes() const noexcept
    {
        return strikes_;
    }

    std::size_t VolSurface::numTenors() const noexcept
    {
        return numTenors_;
    }

    std::size_t VolSurface::numStrikes() const noexcept
    {
        return numStrikes_;
    }

}
