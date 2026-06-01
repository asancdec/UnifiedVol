// SPDX-License-Identifier: Apache-2.0

#include <cmath>

#include "Math/Functions/Detail/JackelDeclare.hpp"
#include "Math/Functions/Volatility.hpp"

namespace uv::math::vol::detail
{

double impliedVolJackelCall(double callPrice, double t, double dF, double F, double K)
{
    return implied_volatility_from_a_transformed_rational_guess(
        callPrice / dF,
        F,
        K,
        t,
        1
    );
}

} // namespace uv::math::vol::detail