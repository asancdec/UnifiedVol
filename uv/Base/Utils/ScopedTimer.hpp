// SPDX-License-Identifier: Apache-2.0

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