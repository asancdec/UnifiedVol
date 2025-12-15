// SPDX-License-Identifier: Apache-2.0
/*
 * File:        CSVRead.cpp
 * Author:      Alvaro Sanchez de Carlos
 * Created:     2025-12-08
 *
 * Description:
 *   [Brief description of what this file declares or implements.]
 *
 * Copyright (c) 2025 Alvaro Sanchez de Carlos
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


#include "Utils/Aux/Errors.hpp"       
#include "Utils/IO/CSVRead.hpp"
#include "Core/MarketData.hpp"
#include "Core/VolSurface.hpp"

#include <fstream>
#include <sstream>
#include <iostream>
#include <cctype>
#include <string>

namespace uv::utils
{
    core::VolSurface readVolSurface(const std::string& filename, const core::MarketData& mktData)
    {
        // Open the CSV file for reading
        std::ifstream file(filename);
        if (!file.is_open())
        {
            raise(ErrorCode::FileIO, "Could not open file " + filename);
        }

        // Small helpers (kept local to preserve file style)
        auto trim = [](std::string& s)
            {
                auto issp = [](unsigned char c) { return std::isspace(c); };
                while (!s.empty() && issp(static_cast<unsigned char>(s.front()))) s.erase(s.begin());
                while (!s.empty() && issp(static_cast<unsigned char>(s.back())))  s.pop_back();
            };

        auto parseCell = [&](const std::string& raw, Real& out, bool strict = true) -> bool
            {
                std::string s = raw;
                trim(s);
                if (s.empty()) return false;

                bool percent = false;
                if (!s.empty() && s.back() == '%')
                {
                    percent = true; s.pop_back(); trim(s);
                    if (s.empty()) { if (strict) raise(ErrorCode::DataFormat, "Lonely % cell"); else return false; }
                }

                size_t idx = 0;
                try { out = std::stod(s, &idx); }
                catch (...) { if (strict) raise(ErrorCode::DataFormat, "Non-numeric cell: \"" + raw + "\""); else return false; }

                if (idx != s.size())
                {
                    if (strict) raise(ErrorCode::DataFormat, "Trailing junk in \"" + raw + "\"");
                    return false;
                }

                if (percent) out *= 0.01;
                return true;
            };

        auto splitComma = [](const std::string& line)
            {
                std::vector<std::string> cells;
                std::stringstream ss(line);
                std::string cell;
                while (std::getline(ss, cell, ',')) cells.push_back(cell);
                if (!line.empty() && line.back() == ',') cells.emplace_back(""); // trailing comma
                return cells;
            };

        // --------------------------------------------------------------------------
        // Read header row
        // Expect: [label or empty], mny1, mny2, ...
        // --------------------------------------------------------------------------
        std::string line{};
        if (!std::getline(file, line))
        {
            raise(ErrorCode::DataFormat, "CSV file is empty: " + filename);
        }

        const auto header = splitComma(line);
        if (header.size() < 2)
        {
            raise(ErrorCode::DataFormat, "Header must have at least 2 columns (blank/label + one moneyness)");
        }

        Vector<Real> mny;
        mny.reserve(header.size() - 1);

        for (std::size_t j = 1; j < header.size(); ++j) // skip first header cell
        {
            Real v{};
            if (!parseCell(header[j], v, /*strict=*/true))
                raise(ErrorCode::DataFormat, "Non-numeric moneyness at header col " + std::to_string(j + 1));
            mny.push_back(v);
        }

        if (mny.empty())
        {
            raise(ErrorCode::DataFormat, "No moneyness columns found in header");
        }

        // --------------------------------------------------------------------------
        // Read data rows: maturity, vol_1, vol_2, ..., vol_N
        // --------------------------------------------------------------------------
        Vector<Real> tenors{};
        Matrix<Real> vols{};

        std::size_t lineNo = 1; // already consumed header
        while (std::getline(file, line))
        {
            ++lineNo;

            std::string probe = line;
            trim(probe);
            if (probe.empty()) continue; // skip blank lines

            const auto cells = splitComma(line);
            if (cells.size() < 2)
            {
                raise(ErrorCode::DataFormat, "Row " + std::to_string(lineNo) + " has fewer than 2 columns");
            }

            // First column is maturity
            Real T{};
            if (!parseCell(cells[0], T, /*strict=*/true))
            {
                raise(ErrorCode::DataFormat, "Missing/invalid maturity at row " + std::to_string(lineNo));
            }

            // Remaining must match mny.size()
            if (cells.size() - 1 < mny.size())
            {
                raise(ErrorCode::DataFormat,
                    "Row " + std::to_string(lineNo) + " has only " +
                    std::to_string(cells.size() - 1) + " vol columns; expected " +
                    std::to_string(mny.size()));
            }

            Vector<Real> row;
            row.reserve(mny.size());
            for (std::size_t j = 0; j < mny.size(); ++j)
            {
                Real sigma{};
                if (!parseCell(cells[1 + j], sigma, /*strict=*/true))
                {
                    raise(ErrorCode::DataFormat,
                        "Non-numeric vol at row " + std::to_string(lineNo) +
                        ", col " + std::to_string(1 + j + 1));
                }
                row.push_back(sigma);
            }

            tenors.push_back(T);
            vols.push_back(std::move(row));
        }

        if (tenors.empty())
        {
            raise(ErrorCode::DataFormat, "CSV file has no data rows: " + filename);
        }

        // Final sanity: shapes must match
        if (vols.size() != tenors.size())
        {
            raise(ErrorCode::DataFormat, "Internal error: vols rows != maturities count");
        }
        for (std::size_t i = 0; i < vols.size(); ++i)
        {
            if (vols[i].size() != mny.size())
            {
                raise(ErrorCode::DataFormat,
                    "Ragged row at data row " + std::to_string(i + 2) + // +2 for header + 1-based
                    ": got " + std::to_string(vols[i].size()) +
                    " vols, expected " + std::to_string(mny.size()));
            }
        }

        // Construct VolSurface using the MarketData-based ctor
        return core::VolSurface(tenors, mny, vols, mktData);
    }
} // namespace uv::utils