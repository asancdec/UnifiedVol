// SPDX-License-Identifier: Apache-2.0
/*
 * File:        VolSurface.hpp
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

#include "Core/MarketData.hpp"
#include "Core/Matrix/Matrix.hpp"

#include "concepts"
#include <cstddef>

namespace uv::core
{
    /**
     * @brief Volatility surface container on a fixed (tenor × moneyness) grid.
     *
     * Holds an implied volatility surface with rows indexed by tenor and columns
     * indexed by moneyness (and the corresponding strike grid). The class also
     * stores commonly used derived quantities: forwards, log(F/K), total variance,
     * and Black–Scholes call prices.
     *
     * Indexing convention:
     * - Row i: tenor tenors_[i]
     * - Col j: moneyness mny_[j] and strike strikes_[j]
     *
     * Rates and dividends are currently treated as flat term structures, stored
     * as one value per tenor.
     */
    template <std::floating_point T>
    class VolSurface
    {
    private:

        //--------------------------------------------------------------------------
        // Member variables
        //--------------------------------------------------------------------------

        Vector<T> tenors_;        ///< Tenors (years), one per surface row.
        std::size_t numTenors_;   ///< Number of tenors (rows).

        Vector<T> mny_;           ///< Moneyness grid (%S), one per surface column.
        std::size_t numStrikes_;  ///< Number of strikes/moneyness points (cols).

        Matrix<T> volMatrix_;     ///< Implied vols [tenor][strike].

        T S_;                     ///< Spot price.
        Vector<T> rates_;         ///< Risk-free rates per tenor (flat per row).
        Vector<T> dividends_;     ///< Dividend yields per tenor (flat per row).

        Vector<T> strikes_;       ///< Strike grid derived from mny_ and S_.
        Vector<T> forwards_;      ///< Forward prices per tenor.

        Matrix<T> callPrices_;    ///< Black–Scholes call prices [tenor][strike].
        Matrix<T> logKFMatrix_;   ///< log(F/K) values [tenor][strike] (or header row).
        Matrix<T> totVarMatrix_;  ///< Total variance w = vol^2 * T [tenor][strike].

        //--------------------------------------------------------------------------
        // Internal helpers
        //--------------------------------------------------------------------------

        /**
         * @brief Compute forwards_ from S_, rates_, dividends_, and tenors_.
         *
         * F_i = S * exp((r_i - q_i) * T_i)
         */
        void setForwards_() noexcept;

        /**
         * @brief Compute call prices using Black-Scholes formula
         */
        void setCallPrices_();

        /**
         * @brief Compute logKFMatrix_ using forwards_ and strikes_.
         *
         * logKFMatrix_[i][j] = log(forwards_[i] / strikes_[j])
         */
        void setLogKFMatrix_();

        /**
         * @brief Compute totVarMatrix_ from an implied vol matrix.
         *
         * totVarMatrix_[i][j] = volMatrix[i][j]^2 * tenors_[i]
         */
        void setTotVar_(const Matrix<T>& volMatrix);

        /**
         * @brief Compute volMatrix_ from a total variance matrix.
         *
         * volMatrix_[i][j] = sqrt(totVarMatrix[i][j] / tenors_[i])
         */
        void setVolFromVar_(const Matrix<T>& totVarMatrix);

        /**
         * @brief Compute volMatrix_ from Call Prices solving for
         * implied volatility using the Black-Scholes formula.
         */
        void setVolFromPrices_(const Matrix<T>& callPrices);

    public:
        //--------------------------------------------------------------------------
        // Initialization
        //--------------------------------------------------------------------------

        VolSurface() = delete;

        /**
         * @brief Construct a surface from an implied vol matrix.
         *
         * Initializes the tenor/moneyness grid, spot and term-structure data, then
         * builds derived quantities (strikes, forwards, log(F/K), total variance,
         * and Black–Scholes call prices).
         *
         * @param tenors    Tenors in years (rows).
         * @param mny       Moneyness grid (columns).
         * @param volMatrix Implied vols [tenor][strike].
         * @param mktData   Market data (spot, rates, dividends).
         */
        explicit VolSurface(Vector<T> tenors,
            Vector<T> mny,
            Matrix<T> volMatrix,
            const MarketData<T>& mktData);

        //--------------------------------------------------------------------------
        // Setters
        //--------------------------------------------------------------------------

        /**
         * @brief Set total variance matrix and update implied vols.
         *
         * Stores the given total variance surface and recomputes volMatrix_ from it.
         */
        void setTotVar(const Matrix<T>& totalVarMatrix);

        /**
         * @brief Set call price matrix and update implied vols.
         *
         * Stores the given call price surface and recomputes volMatrix_ by inverting
         * Black–Scholes at each grid point.
         */
        void setCallPrices(const Matrix<T>& callPrices);

        //--------------------------------------------------------------------------
        // Getters
        //--------------------------------------------------------------------------

        const Vector<T>& tenors() const noexcept;
        std::size_t numTenors() const noexcept;
        const Vector<T>& mny() const noexcept;
        std::size_t numStrikes() const noexcept;
        const Matrix<T>& volMatrix() const noexcept;
        T S() const noexcept;
        const Vector<T>& rates() const noexcept;
        const Vector<T>& dividends() const noexcept;
        const Vector<T>& strikes() const noexcept;
        const Vector<T>& forwards() const noexcept;
        const Matrix<T>& callPrices() const noexcept;
        const Matrix<T>& logKFMatrix() const noexcept;
        const Matrix<T>& totVarMatrix() const noexcept;

        //--------------------------------------------------------------------------
        // Printing
        //--------------------------------------------------------------------------

        /**
         * @brief Print implied vol surface to the console/logger.
         *
         * @param valuePrec Decimal precision for matrix values.
         * @param mnyFlag   If true, header is moneyness; otherwise header is log(F/K).
         */
        void printVol(unsigned int valuePrec = 5,
            bool mnyFlag = true) const noexcept;

        /**
         * @brief Print total variance surface to the console/logger.
         *
         * @param valuePrec Decimal precision for matrix values.
         * @param mnyFlag   If true, header is moneyness; otherwise header is log(F/K).
         */
        void printTotVar(unsigned int valuePrec = 5,
            bool mnyFlag = true) const noexcept;

        /**
         * @brief Print Black–Scholes call price surface to the console/logger.
         *
         * @param valuePrec Decimal precision for matrix values.
         * @param mnyFlag   If true, header is moneyness; otherwise header is log(F/K).
         */
        void printBSCall(unsigned int valuePrec = 5,
            bool mnyFlag = true) const noexcept;
    };
}

#include "VolSurface.inl"