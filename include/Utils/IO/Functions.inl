// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Functions.inl
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

#include "Utils/IO/Log.hpp"

#include <iomanip>
#include <sstream>
#include <string_view>

namespace uv::utils
{
    template <typename HeaderVec, typename RowLabels, typename Matrix>
    void printMatrix(std::string_view title,
        const HeaderVec& header,
        const RowLabels& rowLabels,
        const Matrix& M,
        unsigned int headerPrec,
        unsigned int rowLabelPrec,
        unsigned int valuePrec) noexcept
    {
        std::ostringstream oss;
        oss << '\n' << title << '\t';

        // Header row
        oss << std::fixed << std::setprecision(headerPrec);
        for (const auto& h : header)
            oss << h << '\t';
        oss << '\n';

        // Body
        for (std::size_t i = 0; i < M.size(); ++i)
        {
            oss << std::fixed << std::setprecision(rowLabelPrec)
                << rowLabels[i] << '\t';

            oss << std::fixed << std::setprecision(valuePrec);
            for (const auto& v : M[i])
                oss << v << '\t';

            oss << '\n';
        }

        UV_INFO(oss.str());
    }
} // namespace uv::utils