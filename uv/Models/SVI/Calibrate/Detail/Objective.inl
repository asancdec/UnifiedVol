// SPDX-License-Identifier: Apache-2.0

#include <cmath>

namespace uv::models::svi::detail
{
template <opt::nlopt::Algorithm Algo>
void setMinObjective(opt::nlopt::Optimizer<4, Algo>& optimizer, ObjectiveContexts& ctx)
{
    optimizer.setMinObjective(&objectiveThunk, &ctx);
}

} // namespace uv::models::svi::detail
