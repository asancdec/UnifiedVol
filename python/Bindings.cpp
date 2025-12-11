// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Bindings.cpp
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


#include "Math/Interpolation.hpp"
#include "Utils/Types.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;
using ld = long double;


static py::array_t<ld, py::array::c_style | py::array::forcecast>
require_1d_ndarray(py::handle obj, const char* name) {
    if (!py::isinstance<py::array>(obj))
        throw py::type_error(std::string(name) + " must be a numpy.ndarray (1-D).");
    auto arr = py::array_t<ld, py::array::c_style | py::array::forcecast>::ensure(obj);
    if (!arr || arr.ndim() != 1)
        throw py::value_error(std::string(name) + " must be a 1-D numpy.ndarray convertible to long double.");
    return arr;
}

static py::array_t<ld, py::array::c_style | py::array::forcecast>
require_2d_ndarray(py::handle obj, const char* name) {
    if (!py::isinstance<py::array>(obj))
        throw py::type_error(std::string(name) + " must be a numpy.ndarray (2-D).");
    auto arr = py::array_t<ld, py::array::c_style | py::array::forcecast>::ensure(obj);
    if (!arr || arr.ndim() != 2)
        throw py::value_error(std::string(name) + " must be a 2-D numpy.ndarray convertible to long double.");
    return arr;
}

static std::vector<ld> vec_from_1d(py::array_t<ld, py::array::c_style | py::array::forcecast> a) {
    const ld* ptr = static_cast<const ld*>(a.data());
    return std::vector<ld>(ptr, ptr + a.size());
}

static std::vector<std::vector<ld>> mat_from_2d(py::array_t<ld, py::array::c_style | py::array::forcecast> a) {
    ssize_t r = a.shape(0), c = a.shape(1);
    const ld* ptr = static_cast<const ld*>(a.data());
    std::vector<std::vector<ld>> M(static_cast<size_t>(r), std::vector<ld>(static_cast<size_t>(c)));
    for (ssize_t i = 0; i < r; ++i)
        std::copy(ptr + i * c, ptr + i * c + c, M[static_cast<size_t>(i)].begin());
    return M;
}

static py::array_t<ld> make_1d(py::ssize_t n, const std::vector<ld>& v) {
    py::array_t<ld> out(n);
    std::copy(v.data(), v.data() + v.size(), out.mutable_data());
    return out;
}

static py::array_t<ld> make_2d(py::ssize_t r, py::ssize_t c, const std::vector<std::vector<ld>>& M) {
    py::array_t<ld> out({ r, c });
    ld* dst = out.mutable_data();
    for (py::ssize_t i = 0; i < r; ++i)
        std::copy(M[static_cast<size_t>(i)].begin(), M[static_cast<size_t>(i)].end(), dst + i * c);
    return out;
}

PYBIND11_MODULE(interp, m) {
    m.doc() = "Minimal numpy-only bindings for uv::math::interp";

    m.def("pchip_interp",
        [](py::handle x_obj, py::handle xs_obj, py::handle ys_obj) {
            auto x = require_1d_ndarray(x_obj, "x");
            auto xs = require_1d_ndarray(xs_obj, "xs");
            auto ys = require_1d_ndarray(ys_obj, "ys");
            auto out = uv::math::interp::pchipInterp<ld>(vec_from_1d(x), vec_from_1d(xs), vec_from_1d(ys));
            return make_1d(x.size(), out);
        },
        py::arg("x"), py::arg("xs"), py::arg("ys"));

    m.def("pchip_interp_2d",
        [](py::handle x_obj, py::handle y_obj, py::handle xs_obj, py::handle ys_obj, py::handle zs_obj) {
            auto x = require_1d_ndarray(x_obj, "x");
            auto y = require_1d_ndarray(y_obj, "y");
            auto xs = require_1d_ndarray(xs_obj, "xs");
            auto ys = require_1d_ndarray(ys_obj, "ys");
            auto zs = require_2d_ndarray(zs_obj, "zs");

            // shape check: zs must be (len(ys), len(xs))
            if (zs.shape(0) != ys.size() || zs.shape(1) != xs.size())
                throw py::value_error("zs shape must be (len(ys), len(xs)).");

            auto out = uv::math::interp::pchipInterp2D<ld>(
                vec_from_1d(x), vec_from_1d(y), vec_from_1d(xs), vec_from_1d(ys), mat_from_2d(zs)
            );

            // result expected as matrix (rows=len(y), cols=len(x))
            return make_2d(static_cast<py::ssize_t>(out.size()), static_cast<py::ssize_t>(out.empty() ? 0 : out[0].size()), out);
        },
        py::arg("x"), py::arg("y"), py::arg("xs"), py::arg("ys"), py::arg("zs"));

    m.def("hermite_interp",
        [](py::handle x_obj, py::handle xs_obj, py::handle ys_obj, py::handle dydx_obj) {
            auto x = require_1d_ndarray(x_obj, "x");
            auto xs = require_1d_ndarray(xs_obj, "xs");
            auto ys = require_1d_ndarray(ys_obj, "ys");
            auto ddx = require_1d_ndarray(dydx_obj, "dydx");
            auto out = uv::math::interp::hermiteSplineInterp<ld>(vec_from_1d(x), vec_from_1d(xs), vec_from_1d(ys), vec_from_1d(ddx));
            return make_1d(x.size(), out);
        },
        py::arg("x"), py::arg("xs"), py::arg("ys"), py::arg("dydx"));

    m.def("pchip_derivatives",
        [](py::handle xs_obj, py::handle ys_obj) {
            auto xs = require_1d_ndarray(xs_obj, "xs");
            auto ys = require_1d_ndarray(ys_obj, "ys");
            auto out = uv::math::interp::pchipDerivatives<ld>(vec_from_1d(xs), vec_from_1d(ys));
            return make_1d(xs.size(), out);
        },
        py::arg("xs"), py::arg("ys"));
}