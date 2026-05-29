// SPDX-License-Identifier: Apache-2.0
/*
 * Copyright (c) 2025 �lvaro S�nchez de Carlos
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