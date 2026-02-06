// SPDX-License-Identifier: Apache-2.0
/*
 * Copyright (c) 2025 Álvaro Sánchez de Carlos
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

#include <concepts>
#include <cstddef>
#include <filesystem>
#include <optional>
#include <source_location>
#include <span>
#include <string_view>

namespace uv::errors
{
template <std::floating_point T>
void validateFinite(
    std::span<const T> xs,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <std::floating_point T>
void validateFinite(
    T x,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <std::floating_point T>
void validateNonNegative(
    std::span<const T> xs,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <std::floating_point T>
void validateNonNegative(
    T x,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <std::floating_point T>
void validatePositive(
    std::span<const T> xs,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <std::floating_point T>
void validatePositive(
    T x,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <std::floating_point T>
void validateEqualOrGreater(
    std::span<const T> threshold,
    std::span<const T> xs,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <std::floating_point T>
void validateEqualOrGreater(
    T threshold,
    std::span<const T> xs,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <std::floating_point T>
void validateEqualOrGreater(
    T threshold,
    T x,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <std::floating_point T>
void validateGreater(
    T threshold,
    std::span<const T> xs,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <std::floating_point T>
void validateGreater(
    T threshold,
    T x,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <std::floating_point T>
void validateStrictlyIncreasing(
    std::span<const T> xs,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <typename T>
void validateNonEmpty(
    std::span<const T> xs,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

template <typename T>
void validateSet(
    const std::optional<T>& x,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

void validateSameSize(
    std::size_t first,
    std::size_t second,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

void validateMinSize(
    std::size_t xSize,
    std::size_t minSize,
    std::string_view what,
    std::source_location loc = std::source_location::current()
);

void validateDirCreated(
    bool ok,
    const std::filesystem::path& dir,
    std::source_location loc = std::source_location::current()
);

void validateFileOpened(
    bool ok,
    const std::filesystem::path& file,
    std::source_location loc = std::source_location::current()
);

} // namespace uv::errors

#include <Base/Errors/Detail/Validate.inl>
