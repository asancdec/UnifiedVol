// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Base/Types.hpp"

#include <cctype>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>

namespace uv::tests::json
{
struct Value
{
    enum class Type
    {
        Number,
        String,
        Object,
        Array
    };

    Type type{Type::Object};
    double number{};
    std::string string;
    std::map<std::string, Value> object;
    Vector<Value> array;

    const Value& at(const std::string& key) const
    {
        const auto it = object.find(key);
        if (it == object.end())
            throw std::runtime_error("Missing JSON key: " + key);
        return it->second;
    }

    double asNumber() const
    {
        if (type != Type::Number)
            throw std::runtime_error("Expected numeric JSON value");
        return number;
    }

    const std::string& asString() const
    {
        if (type != Type::String)
            throw std::runtime_error("Expected string JSON value");
        return string;
    }
};

class Parser
{
  public:
    explicit Parser(std::string text)
        : text_{std::move(text)}
    {
    }

    Value parse()
    {
        Value value{parseValue()};
        skipWhitespace();
        if (pos_ != text_.size())
            throw std::runtime_error("Unexpected trailing content in JSON");
        return value;
    }

  private:
    Value parseValue()
    {
        skipWhitespace();
        if (pos_ == text_.size())
            throw std::runtime_error("Unexpected end of JSON");

        if (text_[pos_] == '{')
            return parseObject();
        if (text_[pos_] == '[')
            return parseArray();
        if (text_[pos_] == '"')
            return parseStringValue();
        return parseNumber();
    }

    Value parseObject()
    {
        Value value;
        value.type = Value::Type::Object;
        expect('{');
        skipWhitespace();
        if (consume('}'))
            return value;

        while (true)
        {
            const std::string key{parseString()};
            expect(':');
            value.object.emplace(key, parseValue());
            skipWhitespace();
            if (consume('}'))
                return value;
            expect(',');
        }
    }

    Value parseArray()
    {
        Value value;
        value.type = Value::Type::Array;
        expect('[');
        skipWhitespace();
        if (consume(']'))
            return value;

        while (true)
        {
            value.array.emplace_back(parseValue());
            skipWhitespace();
            if (consume(']'))
                return value;
            expect(',');
        }
    }

    Value parseNumber()
    {
        skipWhitespace();
        const char* begin{text_.c_str() + pos_};
        char* end{};
        const double number{std::strtod(begin, &end)};
        if (end == begin)
            throw std::runtime_error("Expected numeric JSON value");
        pos_ += static_cast<std::size_t>(end - begin);

        Value value;
        value.type = Value::Type::Number;
        value.number = number;
        return value;
    }

    Value parseStringValue()
    {
        Value value;
        value.type = Value::Type::String;
        value.string = parseString();
        return value;
    }

    std::string parseString()
    {
        skipWhitespace();
        expect('"');
        std::string out;
        while (pos_ < text_.size())
        {
            const char c{text_[pos_++]};
            if (c == '"')
                return out;
            if (c == '\\')
            {
                if (pos_ == text_.size())
                    throw std::runtime_error("Unfinished escape in JSON string");
                out.push_back(text_[pos_++]);
                continue;
            }
            out.push_back(c);
        }
        throw std::runtime_error("Unterminated JSON string");
    }

    bool consume(const char c)
    {
        skipWhitespace();
        if (pos_ < text_.size() && text_[pos_] == c)
        {
            ++pos_;
            return true;
        }
        return false;
    }

    void expect(const char c)
    {
        if (!consume(c))
            throw std::runtime_error("Unexpected token in JSON");
    }

    void skipWhitespace()
    {
        while (pos_ < text_.size() &&
               std::isspace(static_cast<unsigned char>(text_[pos_])) != 0)
            ++pos_;
    }

    std::string text_;
    std::size_t pos_{};
};

inline Value read(const std::filesystem::path& path)
{
    std::ifstream in{path};
    if (!in)
        throw std::runtime_error("Could not open JSON file: " + path.string());

    std::string text{
        std::istreambuf_iterator<char>{in},
        std::istreambuf_iterator<char>{}
    };
    return Parser{std::move(text)}.parse();
}
} // namespace uv::tests::json
