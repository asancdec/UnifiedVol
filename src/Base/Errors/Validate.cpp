// SPDX-License-Identifier: Apache-2.0

#include "Base/Errors/Validate.hpp"
#include "Base/Errors/Errors.hpp"

#include <format>

namespace uv::errors::validate
{

void sameSize(
    std::size_t a,
    std::size_t b,
    std::string_view what,
    std::source_location loc
)
{
    if (a != b) [[unlikely]]
    {
        raise(
            ErrorCode::InvalidArgument,
            std::format("{} size mismatch: {} != {}", what, a, b),
            loc
        );
    }
}

void state(bool ok, std::string_view message, std::source_location loc)
{
    if (!ok) [[unlikely]]
    {
        raise(ErrorCode::InvalidState, message, loc);
    }
}

void dirCreated(bool ok, const std::filesystem::path& dir, std::source_location loc)
{
    if (!ok) [[unlikely]]
    {
        raise(
            ErrorCode::FileIO,
            std::format("Failed to create directory: {}", dir.string()),
            loc
        );
    }
}

void fileOpened(bool ok, const std::filesystem::path& file, std::source_location loc)
{
    if (!ok) [[unlikely]]
    {
        raise(
            ErrorCode::FileIO,
            std::format("Unable to open file: {}", file.string()),
            loc
        );
    }
}
} // namespace uv::errors::validate
