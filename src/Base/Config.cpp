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

#include <Base/Config.hpp>
#include <Base/Utils/Detail/Log.hpp>

namespace uv
{
void initialize(const Config& cfg)
{
    detail::applyLogConfig(cfg);
}
namespace detail
{
void applyLogConfig(const Config& cfg)
{
    utils::Log& log = utils::Log::instance();

    log.enableConsole(cfg.logToConsole);

    if (cfg.logToFile)
    {
        log.setFile(cfg.logFile);
    }
}
} // namespace detail
} // namespace uv