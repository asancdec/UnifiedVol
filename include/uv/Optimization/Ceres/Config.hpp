// SPDX-License-Identifier: Apache-2.0
/*
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

#include <Base/Types.hpp>

#include <string_view>

namespace uv::opt::ceres
{

enum class TrustRegionStrategy
{
    LevenbergMarquardt,
    Dogleg
};

enum class LinearSolver
{
    DenseQR,
    DenseNormalCholesky,
    SparseNormalCholesky,
    SparseSchur
};

enum class Loss
{
    None,
    Huber,
    Cauchy,
    SoftL1
};

enum class GradientMode
{
    Analytic,
    NumericForward,
    NumericCentral
};

struct Config
{

    unsigned maxEval;
    double functionTol;
    double paramTol;
    double gradientTol;

    double lossScale;

    Vector<std::string_view> paramNames;

    bool verbose{false};
};

} // namespace uv::opt::ceres