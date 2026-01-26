// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Surface.hpp
 * Author:      Alvaro Sanchez de Carlos
 * Created:     2026-01-22
 *
 * Description:
 *   Local Volatility Surface data class
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

#include "Core/Matrix/Matrix.hpp"
#include "Core/Types.hpp"

#include <concepts>
#include <span>

namespace uv::models::localvol
{
template <std::floating_point T> class Surface
{
  private:
    Vector<T> tenors_;
    core::Matrix<T> logKF_;
    core::Matrix<T> localVar_;

    // Matrix<T> dydx_;

  public:
    Surface() = delete;

    explicit Surface(Vector<T> tenors, core::Matrix<T> logKF, core::Matrix<T> localVar);
};

} // namespace uv::models::localvol

#include "Surface.inl"