// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Pricer.hpp
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

#include "Utils/Types.hpp"
#include "Models/Heston/Params.hpp"
#include "Models/Heston/Config.hpp"
#include "Math/Quadratures/TanHSinH.hpp"
#include "Core/VolSurface.hpp"

#include <memory>
#include <array>
#include <optional>    
#include <cstddef>  
#include <tuple>


namespace uv::models::heston
{     
    /**
     * @brief Heston call option pricer using a damped Fourier / contour-shift representation.
     *
     * Uses a Tanh-Sinh quadrature backend (provided externally) to evaluate the
     * Andersen–Lake contour-deformation integral formulation. @cite AndersenLake2018
     *
     * @tparam N Quadrature order / compile-time setting used by TanHSinH.
     */
    template <std::size_t N>
    class Pricer
    {
    private:

        //--------------------------------------------------------------------------
        // Internal struct to save intermediate calculations
        //--------------------------------------------------------------------------
        /// @cond INTERNAL
        /**
         * @brief Cached characteristic-function intermediates for fast gradient evaluation.
         *
         * @details This is an implementation detail of callPriceWithGradient(...).
         */
        struct CharFunCache
        {
            Complex<Real> psi;                 ///< psi(u) := exp(A + v0 * B)
            Complex<Real> A;                   ///< A term in CF exponent
            Complex<Real> B;                   ///< B term in CF exponent
            Complex<Real> beta;                ///< beta := kappa − sigma*rho*(i*u)
            Complex<Real> D;                   ///< D := sqrt(beta^2 + sigma^2*u(u+i))
            Complex<Real> DT;                  ///< DT := D*T
            Complex<Real> betaPlusD;           ///< beta + D
            Complex<Real> betaMinusD;          ///< beta − D (stable branch)
            Complex<Real> ui;                  ///< ui := u*i
            Complex<Real> kFac;                ///< kFac := (kappa*theta)/sigma^2
            Real invSigma2;                    ///< 1 / sigma^2
            Real kappaTheta;                   ///< kappa * theta
            Real sigma2;                       ///< sigma^2

            Complex<Real> uu;                  ///< uu := u(u+i)
            Complex<Real> eDT;                 ///< eDT := exp(-DT)
            Complex<Real> g;                   ///< g := (beta−D)/(beta+D)
            Complex<Real> Q;                   ///< Q := 1 − g*eDT
            Complex<Real> invQ;                ///< 1 / Q
            Complex<Real> invQ2;               ///< 1 / Q^2
            Complex<Real> R;                   ///< R := 1 − g
            Complex<Real> S;                   ///< S := (beta−D)*T − 2*log(Q/R)
            Complex<Real> fracB;               ///< fracB := (1 − eDT) / Q
            Complex<Real> denomG;              ///< denomG := (beta + D)^2
            Complex<Real> betaMinusDinvSigma2; ///< (betaMinusD)/sigma^2
        };
        /// @endcond

        //--------------------------------------------------------------------------
        // Private member variables
        //--------------------------------------------------------------------------
        std::optional<Params> params_;                     ///< Stored model parameters (optional).
        std::shared_ptr<const math::TanHSinH<N>> quad_;    ///< Quadrature engine (shared).
        const Config config_;                              ///< Pricer configuration (alphas, eps, etc).

        //--------------------------------------------------------------------------
        // Pricing helpers
        //--------------------------------------------------------------------------
        /**
         * @brief Residues arising from the contour shift.
         * @param alpha Damping parameter.
         * @param F Forward price.
         * @param K Strike.
         * @return Residue contribution (depends on alpha region).
         */
        static Real getResidues(Real alpha, Real F, Real K) noexcept;

        /**
         * @brief Select ITM/OTM damping parameter as a function of log-moneyness.
         * @param w log(F/K).
         * @return alpha used in contour integral.
         */
        Real getAlpha(Real w) const noexcept;

        /**
         * @brief Determine contour shift angle.
         *
         * @param kappa Mean reversion speed.
         * @param theta Long-run variance.
         * @param sigma Vol-of-vol.
         * @param rho Correlation.
         * @param v0 Initial variance.
         * @param T Maturity.
         * @param w log(F/K).
         * @return Contour angle phi (radians).
         */
        static Real getPhi(Real kappa, 
            Real theta, 
            Real sigma, 
            Real rho, 
            Real v0, 
            Real T, 
            Real w) noexcept;

        /**
         * @brief Heston characteristic function (Albrecher et al., 2007 stable branch).
         *
         * @param kappa Mean reversion speed.
         * @param theta Long-run variance.
         * @param sigma Vol-of-vol.
         * @param rho Correlation.
         * @param v0 Initial variance.
         * @param T Maturity.
         * @param u Complex argument.
         * @return Characteristic function value psi(u).
         */
        static Complex<Real> charFunction(
            Real kappa, Real theta, Real sigma, Real rho, Real v0, Real T,
            const Complex<Real>& u) noexcept;

        /**
         * @brief Cached characteristic function path for efficient gradients.
         *
         * Same CF as charFunction(...), but returns intermediate values required
         * by callPriceWithGradient(...).
         */
        static CharFunCache charFunctionCal(
            Real kappa, Real theta, Real sigma, Real rho, Real v0, Real T,
            const Complex<Real>& u) noexcept;

    public:

        //--------------------------------------------------------------------------
        // Initialization
        //--------------------------------------------------------------------------
        Pricer() = delete;

        /**
         * @brief Construct a pricer with a provided quadrature engine and configuration.
         * @param quad Quadrature object used to integrate over [0, +inf).
         * @param config Pricer configuration (alphas, eps, tolerances, etc).
         *
         * @note The constructor validates alpha ranges (see implementation).
         */
        explicit Pricer(std::shared_ptr<const math::TanHSinH<N>> quad, const Config& config);

        //--------------------------------------------------------------------------
        // Pricing
        //--------------------------------------------------------------------------
        /**
         * @brief Price a European call using explicit Heston parameters.
         *
         * Andersen–Lake implementation with contour shift and residue correction.
         *
         * @param kappa Mean reversion speed.
         * @param theta Long-run variance.
         * @param sigma Vol-of-vol.
         * @param rho Correlation in [-1,1].
         * @param v0 Initial variance.
         * @param T Maturity.
         * @param F Forward.
         * @param r Risk-free rate used for discounting.
         * @param K Strike.
         * @return Discounted call price.
         */
        Real callPrice(Real kappa,
            Real theta, 
            Real sigma, 
            Real rho, 
            Real v0,
            Real T,
            Real F, 
            Real r, 
            Real K) const noexcept;

        /**
         * @brief Price a European call using parameters stored in the pricer.
         *
         * @throws (via UV_REQUIRE) if parameters are not set (see implementation).
         */
        Real callPrice(Real T,
            Real F, 
            Real r, 
            Real K) const;

        /**
         * @brief Price and gradient w.r.t. model parameters for calibration.
         *
         * Output layout:
         *  - out[0] price
         *  - out[1] dP/dkappa
         *  - out[2] dP/dtheta
         *  - out[3] dP/dsigma
         *  - out[4] dP/drho
         *  - out[5] dP/dv0
         */
        std::array<Real, 6> callPriceWithGradient(Real kappa,
            Real theta,
            Real sigma,
            Real rho,
            Real v0,
            Real T, 
            Real F, 
            Real r,
            Real K) const noexcept;

        //--------------------------------------------------------------------------
        // Setters
        //--------------------------------------------------------------------------
        /**
         * @brief Set internal Heston parameters.
         * @param params Parameter struct.
         */
        void setParams(const Params& params) noexcept;
    };

} // namespace uv::models::heston

#include "Pricer.inl"
