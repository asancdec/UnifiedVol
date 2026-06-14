// SPDX-License-Identifier: Apache-2.0

#include "IO/JSON/Read.hpp"
#include "Support/TempFile.hpp"

#include <filesystem>
#include <gtest/gtest.h>
#include <stdexcept>

TEST(UnitIOJSONRead, ParsesObjectsArraysAndScalarValues)
{
    const auto root = uv::io::json::parse(R"({
      "name": "surface",
      "enabled": true,
      "missing": null,
      "values": [1.25, -2.0e-3, {"nested": "yes"}]
    })");

    ASSERT_EQ(root.type, uv::io::json::Value::Type::Object);
    EXPECT_EQ(root.at("name").asString(), "surface");
    EXPECT_TRUE(root.at("enabled").asBool());
    EXPECT_EQ(root.at("missing").type, uv::io::json::Value::Type::Null);

    const auto& values = root.at("values");
    ASSERT_EQ(values.type, uv::io::json::Value::Type::Array);
    ASSERT_EQ(values.array.size(), 3U);
    EXPECT_DOUBLE_EQ(values.array[0].asNumber(), 1.25);
    EXPECT_DOUBLE_EQ(values.array[1].asNumber(), -2.0e-3);
    EXPECT_EQ(values.array[2].at("nested").asString(), "yes");
}

TEST(UnitIOJSONRead, ParsesWhitespaceAndFalseBoolean)
{
    const auto root = uv::io::json::parse(" \n\t { \"enabled\" : false } \r\n ");

    EXPECT_FALSE(root.at("enabled").asBool());
}

TEST(UnitIOJSONRead, ParsesEscapedStringCharacters)
{
    const auto root =
        uv::io::json::parse(R"({"text": "line\nquote: \" slash: \/ tab:\t cr:\r"})");

    EXPECT_EQ(root.at("text").asString(), "line\nquote: \" slash: / tab:\t cr:\r");
}

TEST(UnitIOJSONRead, ReadsJsonFromFile)
{
    const auto path = uv::tests::writeTempFile(
        "unifiedvol_json_read_fixture.json",
        R"({"threshold": 0.05})"
    );

    EXPECT_DOUBLE_EQ(uv::io::json::read(path).at("threshold").asNumber(), 0.05);
}

TEST(UnitIOJSONRead, RejectsMalformedJsonAndTypeMismatches)
{
    EXPECT_THROW(
        static_cast<void>(uv::io::json::parse(R"({"x": [1, 2,]})")),
        std::runtime_error
    );

    const auto root = uv::io::json::parse(R"({"x": "not a number"})");
    EXPECT_THROW(static_cast<void>(root.at("x").asNumber()), std::runtime_error);
    EXPECT_THROW(static_cast<void>(root.at("missing")), std::runtime_error);

    EXPECT_THROW(static_cast<void>(root.at("x").asBool()), std::runtime_error);
    EXPECT_THROW(
        static_cast<void>(uv::io::json::parse("1").asString()),
        std::runtime_error
    );
    EXPECT_THROW(
        static_cast<void>(uv::io::json::parse("true").asNumber()),
        std::runtime_error
    );
    EXPECT_THROW(
        static_cast<void>(uv::io::json::parse("[1]").at("x")),
        std::runtime_error
    );
}

TEST(UnitIOJSONRead, RejectsMissingEmptyAndMalformedFiles)
{
    EXPECT_THROW(
        static_cast<void>(uv::io::json::read(
            std::filesystem::temp_directory_path() /
            "unifiedvol_json_missing_fixture.json"
        )),
        std::runtime_error
    );

    const auto emptyPath =
        uv::tests::writeTempFile("unifiedvol_json_empty_fixture.json", "");
    EXPECT_THROW(static_cast<void>(uv::io::json::read(emptyPath)), std::runtime_error);
}

TEST(UnitIOJSONRead, RejectsDuplicateObjectKeys)
{
    EXPECT_THROW(
        static_cast<void>(uv::io::json::parse(R"({"x": 1, "x": 2})")),
        std::runtime_error
    );
}

TEST(UnitIOJSONRead, RejectsTrailingContent)
{
    EXPECT_THROW(
        static_cast<void>(uv::io::json::parse(R"({"x": 1} trailing)")),
        std::runtime_error
    );
}

TEST(UnitIOJSONRead, RejectsUnsupportedUnicodeEscapes)
{
    EXPECT_THROW(
        static_cast<void>(uv::io::json::parse(R"({"text": "\u1234"})")),
        std::runtime_error
    );
}

TEST(UnitIOJSONRead, RejectsInvalidNumberForms)
{
    EXPECT_THROW(static_cast<void>(uv::io::json::parse("01")), std::runtime_error);
    EXPECT_THROW(static_cast<void>(uv::io::json::parse("+1")), std::runtime_error);
    EXPECT_THROW(static_cast<void>(uv::io::json::parse("1.")), std::runtime_error);
    EXPECT_THROW(static_cast<void>(uv::io::json::parse("1e")), std::runtime_error);
}
