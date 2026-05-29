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

#include "Base/Utils/StopWatch.hpp"
#include <string_view>

namespace uv::utils
{

template <typename Period = std::milli> class ScopedTimer
{
  public:
    ScopedTimer() noexcept
    {
        watch_.StartStopWatch();
    }

    ScopedTimer(std::string_view label) noexcept
        : label_(label)
    {
        watch_.StartStopWatch();
    }

    ~ScopedTimer() noexcept
    {
        watch_.StopStopWatch();
        watch_.LogTime<Period>(label_);
    }

  private:
    std::string_view label_;
    uv::utils::StopWatch watch_;
};

} // namespace uv::utils