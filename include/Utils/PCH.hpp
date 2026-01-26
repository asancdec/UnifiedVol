// SPDX-License-Identifier: Apache-2.0
/*
 * File:        PCH.hpp
 * Author:      Alvaro Sanchez de Carlos
 * Created:     2025-12-08
 *
 * Description:
 *   Project precompiled header (PCH) that centralizes frequently used core,
 *   model, math, and optimization headers—including Ceres and NLopt—to reduce
 *   compile times and provide a single common include for implementation files.
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

// ---------- Core I/O ----------

#include "Utils/Aux/StopWatch.hpp"
#include "Utils/IO/CSV/Functions.hpp"
#include "Utils/IO/Log.hpp"

// ---------- Core Data ----------

#include "Core/Functions.hpp"
#include "Core/MarketData.hpp"
#include "Core/Matrix/Matrix.hpp"
#include "Core/Types.hpp"
#include "Core/VolSurface.hpp"

// ---------- Models ----------

#include "Models/Heston/Calibrator.hpp"
#include "Models/Heston/Config.hpp"
#include "Models/Heston/Pricer.hpp"
#include "Models/LocalVol/Calibrator.hpp"
#include "Models/LocalVol/Pricer.hpp"
#include "Models/SVI/Functions.hpp"
#include "Models/SVI/Params.hpp"

// ---------- Utils ----------

#include "Utils/Aux/Errors.hpp"

// ---------- Math ----------

#include "Math/Integration/TanHSinH.hpp"
#include "Math/Optimization/Ceres/Config.hpp"
#include "Math/Optimization/Ceres/Policy.hpp"

// ---------- Third-party ----------

#include <ceres/loss_function.h>
#include <ceres/types.h>
#include <nlopt.hpp>

// ---------- STL ----------

#include <cstdlib>
#include <exception>
#include <filesystem>
#include <format>
#include <iostream>
#include <memory>
#include <ratio>
