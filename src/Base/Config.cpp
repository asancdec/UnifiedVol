// SPDX-License-Identifier: Apache-2.0

#include "Base/Config.hpp"
#include "Base/Utils/Detail/Log.hpp"

namespace uv
{
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

void initialize(const Config& cfg)
{
    detail::applyLogConfig(cfg);
}
} // namespace uv
