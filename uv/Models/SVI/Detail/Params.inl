// SPDX-License-Identifier: Apache-2.0

#include "Models/SVI/Math.hpp"

namespace uv::models::svi
{

template <std::floating_point T>
Params<T>::Params(T t_, T a_, T b_, T rho_, T m_, T sigma_) noexcept
    : t(t_),
      a(a_),
      b(b_),
      rho(rho_),
      m(m_),
      sigma(sigma_)
{
}

template <std::floating_point T>
Params<T>::Params(T t_, std::span<const double> params, double atmTotalVariance) noexcept
    : t(t_),
      b(static_cast<T>(params[0])),
      rho(static_cast<T>(params[1])),
      m(static_cast<T>(params[2])),
      sigma(static_cast<T>(params[3]))
{
    a = detail::aParam<T>(static_cast<T>(atmTotalVariance), b, rho, m, sigma);
}

template <std::floating_point T> template <std::floating_point U>
Params<U> Params<T>::as() const noexcept
{
    return Params<U>{
        static_cast<U>(t),
        static_cast<U>(a),
        static_cast<U>(b),
        static_cast<U>(rho),
        static_cast<U>(m),
        static_cast<U>(sigma)
    };
}
} // namespace uv::models::svi