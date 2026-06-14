// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

namespace uv::tests
{
inline std::filesystem::path
writeTempFile(const std::string& name, const std::string& contents)
{
    const auto path = std::filesystem::temp_directory_path() / name;
    std::ofstream out{path};
    if (!out)
        throw std::runtime_error("Could not open temp fixture file: " + path.string());
    out << contents;
    return path;
}
} // namespace uv::tests
