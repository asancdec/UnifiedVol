// SPDX-License-Identifier: Apache-2.0

#include <cmath>

namespace uv::models::svi::detail
{
template <opt::nlopt::Algorithm Algo> void setMinObjective(
    opt::nlopt::Optimizer<4, Algo>& optimizer,
    const ObjectiveContexts& ctx
) noexcept
{
    optimizer.setMinObjective(&objectiveThunk, const_cast<ObjectiveContexts*>(&ctx));
}

} // namespace uv::models::svi::detail
