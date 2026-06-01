// SPDX-License-Identifier: Apache-2.0

namespace uv
{
template <typename To, typename From> Vector<To> convertVector(const Vector<From>& x)
{
    Vector<To> out;
    out.reserve(x.size());

    for (const From& v : x)
        out.push_back(static_cast<To>(v));

    return out;
}

template <typename To, typename From> Vector<To> convertVector(std::span<const From> x)
{
    Vector<To> out;
    out.reserve(x.size());

    for (const From& v : x)
        out.push_back(static_cast<To>(v));

    return out;
}
} // namespace uv