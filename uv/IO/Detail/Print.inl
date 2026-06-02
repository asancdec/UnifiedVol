// SPDX-License-Identifier: Apache-2.0

#include "Base/Macros/Inform.hpp"

#include <iomanip>
#include <sstream>
#include <string_view>

namespace uv::utils
{
template <typename HeaderVec, typename RowLabels, typename Matrix> void printMatrix(
    std::string_view title,
    const HeaderVec& header,
    const RowLabels& rowLabels,
    const Matrix& M,
    unsigned int headerPrec,
    unsigned int rowLabelPrec,
    unsigned int valuePrec
) noexcept
{
    std::ostringstream oss;
    oss << '\n' << title << '\t';

    oss << std::fixed << std::setprecision(headerPrec);
    for (const auto& h : header)
        oss << h << '\t';
    oss << '\n';

    for (std::size_t i = 0; i < M.rows(); ++i)
    {
        oss << std::fixed << std::setprecision(rowLabelPrec) << rowLabels[i] << '\t';

        oss << std::fixed << std::setprecision(valuePrec);
        for (const auto& v : M[i])
            oss << v << '\t';

        oss << '\n';
    }

    INFO(oss.str());
}

template <typename Vector>
void printVector(const Vector& v, unsigned int valuePrec) noexcept
{
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(valuePrec);

    for (const auto& x : v)
        oss << x << '\t';

    oss << '\n';

    INFO(oss.str());
}

} // namespace uv::utils