// SPDX-License-Identifier: Apache-2.0

#include "IO/JSON/Read.hpp"
#include "Base/Errors/Errors.hpp"
#include "Support/TempFile.hpp"

#include <filesystem>
#include <gtest/gtest.h>

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
        uv::errors::UnifiedVolError
    );

    const auto root = uv::io::json::parse(R"({"x": "not a number"})");
    EXPECT_THROW(static_cast<void>(root.at("x").asNumber()), uv::errors::UnifiedVolError);
    EXPECT_THROW(static_cast<void>(root.at("missing")), uv::errors::UnifiedVolError);

    EXPECT_THROW(static_cast<void>(root.at("x").asBool()), uv::errors::UnifiedVolError);
    EXPECT_THROW(
        static_cast<void>(uv::io::json::parse("1").asString()),
        uv::errors::UnifiedVolError
    );
    EXPECT_THROW(
        static_cast<void>(uv::io::json::parse("true").asNumber()),
        uv::errors::UnifiedVolError
    );
    EXPECT_THROW(
        static_cast<void>(uv::io::json::parse("[1]").at("x")),
        uv::errors::UnifiedVolError
    );
}

TEST(UnitIOJSONRead, RejectsMissingEmptyAndMalformedFiles)
{
    EXPECT_THROW(
        static_cast<void>(uv::io::json::read(
            std::filesystem::temp_directory_path() /
            "unifiedvol_json_missing_fixture.json"
        )),
        uv::errors::UnifiedVolError
    );

    const auto emptyPath =
        uv::tests::writeTempFile("unifiedvol_json_empty_fixture.json", "");
    EXPECT_THROW(
        static_cast<void>(uv::io::json::read(emptyPath)),
        uv::errors::UnifiedVolError
    );
}

TEST(UnitIOJSONRead, RejectsDuplicateObjectKeys)
{
    EXPECT_THROW(
        static_cast<void>(uv::io::json::parse(R"({"x": 1, "x": 2})")),
        uv::errors::UnifiedVolError
    );
}

TEST(UnitIOJSONRead, RejectsTrailingContent)
{
    EXPECT_THROW(
        static_cast<void>(uv::io::json::parse(R"({"x": 1} trailing)")),
        uv::errors::UnifiedVolError
    );
}

TEST(UnitIOJSONRead, RejectsUnsupportedUnicodeEscapes)
{
    EXPECT_THROW(
        static_cast<void>(uv::io::json::parse(R"({"text": "\u1234"})")),
        uv::errors::UnifiedVolError
    );
}

TEST(UnitIOJSONRead, RejectsInvalidNumberForms)
{
    EXPECT_THROW(
        static_cast<void>(uv::io::json::parse("01")),
        uv::errors::UnifiedVolError
    );
    EXPECT_THROW(
        static_cast<void>(uv::io::json::parse("+1")),
        uv::errors::UnifiedVolError
    );
    EXPECT_THROW(
        static_cast<void>(uv::io::json::parse("1.")),
        uv::errors::UnifiedVolError
    );
    EXPECT_THROW(
        static_cast<void>(uv::io::json::parse("1e")),
        uv::errors::UnifiedVolError
    );
}

TEST(UnitIOJSONRead, UsesSpecificUnifiedVolErrorCodes)
{
    try
    {
        static_cast<void>(uv::io::json::parse(R"({"x": [1, 2,]})"));
        FAIL() << "Expected malformed JSON to throw";
    }
    catch (const uv::errors::UnifiedVolError& e)
    {
        EXPECT_EQ(e.code(), uv::errors::ErrorCode::DataFormat);
    }

    try
    {
        static_cast<void>(uv::io::json::read(
            std::filesystem::temp_directory_path() /
            "unifiedvol_json_missing_code_fixture.json"
        ));
        FAIL() << "Expected missing JSON file to throw";
    }
    catch (const uv::errors::UnifiedVolError& e)
    {
        EXPECT_EQ(e.code(), uv::errors::ErrorCode::FileIO);
    }
}
