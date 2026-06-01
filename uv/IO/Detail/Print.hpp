// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <concepts>
#include <string>

namespace uv::io
{

template <typename HeaderVec, typename RowLabels, typename Matrix> void printMatrix(
    std::string_view title,
    const HeaderVec& header,
    const RowLabels& rowLabels,
    const Matrix& M,
    unsigned int headerPrec = 2,
    unsigned int rowLabelPrec = 2,
    unsigned int valuePrec = 5
) noexcept;

template <typename Vector>
void printVector(const Vector& v, unsigned int valuePrec = 5) noexcept;

} // namespace uv::io

#include "IO/Detail/Print.inl"
