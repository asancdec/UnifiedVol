// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Base/Types.hpp"

#include <filesystem>
#include <map>
#include <string>

namespace uv::io::json
{
struct Value
{
    enum class Type
    {
        Null,
        Bool,
        Number,
        String,
        Object,
        Array
    };

    Type type{Type::Null};
    bool boolean{};
    double number{};
    std::string string;
    std::map<std::string, Value> object;
    Vector<Value> array;

    const Value& at(const std::string& key) const;
    bool asBool() const;
    double asNumber() const;
    const std::string& asString() const;
};

Value parse(std::string text);
Value read(const std::filesystem::path& path);
} // namespace uv::io::json
