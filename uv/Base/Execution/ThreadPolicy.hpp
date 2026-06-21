// SPDX-License-Identifier: Apache-2.0

#pragma once

namespace uv::execution
{
namespace detail
{
int availableThreads() noexcept;
} // namespace detail

int requestThreads(int);
} // namespace uv::execution
