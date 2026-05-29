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

#include "Base/Config.hpp"
#include "Base/Errors/Errors.hpp"
#include "Base/Types.hpp"
#include "Base/Utils/ScopedTimer.hpp"

#include "Core/Curve.hpp"
#include "Core/MarketData.hpp"
#include "Core/MarketState.hpp"
#include "Core/VolSurface.hpp"

#include "IO/Load.hpp"
#include "IO/Report.hpp"

#include "Optimization/Ceres/Config.hpp"
#include "Optimization/Ceres/Optimizer.hpp"
#include "Optimization/Cost.hpp"
#include "Optimization/NLopt/Config.hpp"
#include "Optimization/NLopt/Optimizer.hpp"

#include "Math/Functions/Black.hpp"
#include "Math/Integration/TanHSinH.hpp"

#include "Models/SVI/BuildSurface.hpp"
#include "Models/SVI/Calibrate/Calibrate.hpp"
#include "Models/SVI/Calibrate/Config.hpp"
#include "Models/SVI/Math.hpp"
#include "Models/SVI/Params.hpp"

#include "Models/Heston/BuildSurface.hpp"
#include "Models/Heston/Calibrate/Calibrate.hpp"
#include "Models/Heston/Params.hpp"
#include "Models/Heston/Price/Pricer.hpp"
