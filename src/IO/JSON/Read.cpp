// SPDX-License-Identifier: Apache-2.0

#include "IO/JSON/Read.hpp"

#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <utility>

namespace uv::io::json
{
namespace
{
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
        if (startsWith("true"))
            return parseBool(true);
        if (startsWith("false"))
            return parseBool(false);
        if (startsWith("null"))
            return parseNull();
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
            const bool inserted{value.object.emplace(key, parseValue()).second};
            if (!inserted)
                throw std::runtime_error("Duplicate JSON key: " + key);
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

    Value parseBool(const bool boolean)
    {
        pos_ += boolean ? 4U : 5U;
        Value value;
        value.type = Value::Type::Bool;
        value.boolean = boolean;
        return value;
    }

    Value parseNull()
    {
        pos_ += 4U;
        Value value;
        value.type = Value::Type::Null;
        return value;
    }

    Value parseNumber()
    {
        skipWhitespace();
        const std::size_t beginPos{pos_};
        consumeJsonNumber();

        const char* begin{text_.c_str() + beginPos};
        char* end{};
        const double number{std::strtod(begin, &end)};
        if (end != text_.c_str() + pos_ || !std::isfinite(number))
            throw std::runtime_error("Expected finite numeric JSON value");

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
            if (static_cast<unsigned char>(c) < 0x20U)
                throw std::runtime_error("Control character in JSON string");
            if (c == '\\')
            {
                out.push_back(parseEscape());
                continue;
            }
            out.push_back(c);
        }
        throw std::runtime_error("Unterminated JSON string");
    }

    char parseEscape()
    {
        if (pos_ == text_.size())
            throw std::runtime_error("Unfinished escape in JSON string");

        const char c{text_[pos_++]};
        switch (c)
        {
        case '"':
        case '\\':
        case '/':
            return c;
        case 'b':
            return '\b';
        case 'f':
            return '\f';
        case 'n':
            return '\n';
        case 'r':
            return '\r';
        case 't':
            return '\t';
        default:
            throw std::runtime_error("Unsupported escape in JSON string");
        }
    }

    void consumeJsonNumber()
    {
        if (pos_ < text_.size() && text_[pos_] == '-')
            ++pos_;

        if (pos_ == text_.size())
            throw std::runtime_error("Expected numeric JSON value");

        if (text_[pos_] == '0')
        {
            ++pos_;
            if (pos_ < text_.size() && isDigit(text_[pos_]))
                throw std::runtime_error("Leading zero in JSON number");
        }
        else if (isDigitOneToNine(text_[pos_]))
        {
            while (pos_ < text_.size() && isDigit(text_[pos_]))
                ++pos_;
        }
        else
        {
            throw std::runtime_error("Expected numeric JSON value");
        }

        if (pos_ < text_.size() && text_[pos_] == '.')
        {
            ++pos_;
            if (pos_ == text_.size() || !isDigit(text_[pos_]))
                throw std::runtime_error("Expected digit after JSON number decimal point"
                );
            while (pos_ < text_.size() && isDigit(text_[pos_]))
                ++pos_;
        }

        if (pos_ < text_.size() && (text_[pos_] == 'e' || text_[pos_] == 'E'))
        {
            ++pos_;
            if (pos_ < text_.size() && (text_[pos_] == '+' || text_[pos_] == '-'))
                ++pos_;
            if (pos_ == text_.size() || !isDigit(text_[pos_]))
                throw std::runtime_error("Expected digit in JSON number exponent");
            while (pos_ < text_.size() && isDigit(text_[pos_]))
                ++pos_;
        }
    }

    static bool isDigit(const char c) noexcept
    {
        return c >= '0' && c <= '9';
    }

    static bool isDigitOneToNine(const char c) noexcept
    {
        return c >= '1' && c <= '9';
    }

    bool startsWith(const std::string& literal) const
    {
        return text_.compare(pos_, literal.size(), literal) == 0;
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
} // namespace

const Value& Value::at(const std::string& key) const
{
    if (type != Type::Object)
        throw std::runtime_error("Expected JSON object");

    const auto it = object.find(key);
    if (it == object.end())
        throw std::runtime_error("Missing JSON key: " + key);
    return it->second;
}

bool Value::asBool() const
{
    if (type != Type::Bool)
        throw std::runtime_error("Expected boolean JSON value");
    return boolean;
}

double Value::asNumber() const
{
    if (type != Type::Number)
        throw std::runtime_error("Expected numeric JSON value");
    return number;
}

const std::string& Value::asString() const
{
    if (type != Type::String)
        throw std::runtime_error("Expected string JSON value");
    return string;
}

Value parse(std::string text)
{
    return Parser{std::move(text)}.parse();
}

Value read(const std::filesystem::path& path)
{
    std::ifstream in{path};
    if (!in)
        throw std::runtime_error("Could not open JSON file: " + path.string());

    std::string text{
        std::istreambuf_iterator<char>{in},
        std::istreambuf_iterator<char>{}
    };
    return parse(std::move(text));
}
} // namespace uv::io::json
