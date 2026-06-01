// SPDX-License-Identifier: Apache-2.0

#include <format>

namespace uv::errors
{
template <typename E> [[noreturn]] inline void
unreachableEnum(E value, std::string_view enumName, std::source_location loc)
requires std::is_enum_v<E>
{
    using U = std::underlying_type_t<E>;

    raise(
        ErrorCode::Unreachable,
        std::format("Unknown {} value: {}", enumName, static_cast<U>(value)),
        loc
    );
}
} // namespace uv::errors