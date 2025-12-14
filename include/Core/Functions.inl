// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Functions.inl
 * Author:      Álvaro Sánchez de Carlos
 * Created:     2025-12-08
 *
 * Description:
 *   [Brief description of what this file declares or implements.]
 *
 * Copyright (c) 2025 Álvaro Sánchez de Carlos
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

#include "Utils/Aux/Errors.hpp"  

#include <algorithm>

namespace uv::core
{
    template<typename T>
    T minValue(const Vector<T>& x)
    {
        UV_REQUIRE
        (
            !x.empty(),
            ErrorCode::InvalidArgument,
            "minValue: input vector is empty"
        );

        return *std::min_element(x.begin(), x.end());
    }

    template<typename T>
    T maxValue(const Vector<T>& x)
    {
        UV_REQUIRE
        (
            !x.empty(),
            ErrorCode::InvalidArgument,
            "maxValue: input vector is empty"
        );

        return *std::max_element(x.begin(), x.end());
    }

} // namespace uv::core


