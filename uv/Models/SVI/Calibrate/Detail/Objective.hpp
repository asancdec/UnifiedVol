// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Models/SVI/Calibrate/Detail/Contexts.hpp"
#include "Optimization/NLopt/Optimizer.hpp"

namespace uv::models::svi::detail
{

template <opt::nlopt::Algorithm Algo>
void setMinObjective(opt::nlopt::Optimizer<4, Algo>& optimizer, ObjectiveContexts& ctx);

[[gnu::hot]] double
objectiveThunk(unsigned /*n*/, const double* x, double* grad, void* data) noexcept;

} // namespace uv::models::svi::detail

#include "Models/SVI/Calibrate/Detail/Objective.inl"
