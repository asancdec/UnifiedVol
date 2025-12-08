
// SPDX-License-Identifier: Apache-2.0
/*
 * File:        double.hpp
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

#include "Math/Optimization/NLopt/Config.hpp"
#include "Utils/Aux/StopWatch.hpp"

#include <nlopt.hpp>
#include <array>
#include <cstddef>

namespace uv::math::opt::nlopt
{

    template <std::size_t N, ::nlopt::algorithm Algo>
    class Optimizer
    {
    private:

        //--------------------------------------------------------------------------
        // Type aliases
        //--------------------------------------------------------------------------
        using NloptFunction = double (*)(unsigned, const double*, double*, void*);


        //--------------------------------------------------------------------------
        // Member variables
        //--------------------------------------------------------------------------	
        // Config & engine
        Config<N>         config_;
        ::nlopt::opt      opt_;
        utils::StopWatch  timer_;

        // Bounds & guess
        std::array<double, N> lowerBounds_;
        std::array<double, N> upperBounds_;
        std::array<double, N> initGuess_;

        // Objective state
        NloptFunction userFn_;
        void* userData_;
        unsigned iterCount_;   

        //--------------------------------------------------------------------------
        // Static variables
        //--------------------------------------------------------------------------	
        static double ObjectiveThunk(unsigned n, const double* x, double* grad, void* p) noexcept;

    public:

        //--------------------------------------------------------------------------
        // Initialization
        //--------------------------------------------------------------------------
        // Constructors
        Optimizer() = delete;
        explicit Optimizer(const Config<N>& config);

        // Return new calibrator object with same settings
        Optimizer<N, Algo> fresh() const noexcept;

        //--------------------------------------------------------------------------
        // Optimization
        //--------------------------------------------------------------------------	
        // Set Initial Guess and Bounds
        void setGuessBounds(std::array<double, N> initGuess,
            std::array<double, N> lowerBounds,
            std::array<double, N> upperBounds) noexcept;

        // Add inequality constraints
        void addInequalityConstraint(NloptFunction c, void* data) noexcept;

        // Set objective function	
        void setMinObjective(NloptFunction f, void* data) noexcept;

        // Run optimization	
        Vector<double> optimize() noexcept;

        //--------------------------------------------------------------------------
        // Getters
        //--------------------------------------------------------------------------
        const double& eps() const noexcept;
        double tol() const noexcept;
    };
}

#include "Optimizer.inl"