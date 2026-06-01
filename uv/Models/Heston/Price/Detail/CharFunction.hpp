// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "Base/Types.hpp"

#include <concepts>

namespace uv::models::heston::price::detail
{

template <std::floating_point T> struct CharFunCache
{
    Complex<T> logPsi;
    Complex<T> A;
    Complex<T> B;
    Complex<T> beta;
    Complex<T> D;
    Complex<T> betaMinusD;

    Complex<T> uu;
    Complex<T> eDT;
    Complex<T> oneMinusEDT;
    Complex<T> g;
    Complex<T> Q;
    Complex<T> R;
    Complex<T> S;
    Complex<T> denomG;
};

template <std::floating_point T> [[gnu::hot]] Complex<T> charFunction(
    T kappa,
    T kappaThetaDivSigma2,
    T sigma2,
    T v0,
    T t,
    Complex<T> tDivTwo,
    Complex<T> sigmaRho,
    Complex<T> u
) noexcept;

template <std::floating_point T> [[gnu::hot]] CharFunCache<T> charFunctionCached(
    T kappa,
    T kappaThetaDivSigma2,
    T sigma2,
    T v0,
    T t,
    Complex<T> tDivTwo,
    Complex<T> sigmaRho,
    Complex<T> u
) noexcept;

} // namespace uv::models::heston::price::detail

#include "Models/Heston/Price/Detail/CharFunction.inl"