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

#include <source_location>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>

namespace uv::errors
{

enum class ErrorCode
{
    InvalidArgument,
    OutOfRange,
    FileIO,
    NotImplemented,
    DataFormat,
    InvalidState,
    CalibrationError,
    LinearAlgebra,
    Unreachable,
    Unknown,
};

std::string_view to_string(ErrorCode c) noexcept;

class UnifiedVolError final : public std::runtime_error
{
  public:
    UnifiedVolError(
        ErrorCode code,
        std::string_view message,
        std::source_location loc = std::source_location::current()
    );

    ErrorCode code() const noexcept;

    const std::string& file() const noexcept;

    const std::string& function() const noexcept;

    unsigned line() const noexcept;

  private:
    static std::string
    make_what(ErrorCode code, std::string_view msg, std::source_location loc);

    ErrorCode code_;
    std::string file_;
    std::string function_;
    unsigned line_;
};

[[noreturn]] void raise(
    ErrorCode code,
    std::string_view message,
    std::source_location loc = std::source_location::current()
);

template <typename E>
[[noreturn]] void unreachableEnum(
    E value,
    std::string_view enumName,
    std::source_location loc = std::source_location::current()
)
requires std::is_enum_v<E>;

} // namespace uv::errors

#include "Base/Errors/Detail/Errors.inl"
