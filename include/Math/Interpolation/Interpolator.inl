// SPDX-License-Identifier: Apache-2.0
/*
 * File:        Interpolator.inl
 * Author:      Alvaro Sanchez de Carlos
 * Created:     2026-01-22
 *
 * Description:
 *   Interpolator class implementation.
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

namespace uv::math::interp
{
	template
	<
		std::floating_point T,
		class DerivPolicy,
		class EvalPolicy
	>
	void Interpolator<T, DerivPolicy, EvalPolicy>::operator()
	(
		std::span<const T> x,
		std::span<const T> xs,
		std::span<const T> ys,
		std::span<const T> dydx,
		std::span<T> y,
		bool doValidate
	) const
		requires HasDerivatives<DerivPolicy, T> && HasEvaluate<EvalPolicy, T>
	{
		eval
		(
			x,
			xs,
			ys,
			dydx,
			y,
			doValidate
		);
	}

	template
	<
		std::floating_point T,
		class DerivPolicy,
		class EvalPolicy
	>
	Vector<T> Interpolator<T, DerivPolicy, EvalPolicy>::operator()
	(
		std::span<const T> x,
		std::span<const T> xs,
		std::span<const T> ys,
		bool doValidate
	) const
		requires HasDerivatives<DerivPolicy, T> && HasEvaluate<EvalPolicy, T>
	{
		Vector<T> dydx(xs.size());
		deriv(xs, ys, dydx, doValidate);

		Vector<T> y(x.size());

		(*this)
		(
			x,
			xs,
			ys,
			dydx,
			y,
			doValidate
		);

		return y;
	}

	template
	<
		std::floating_point T,
		class DerivPolicy,
		class EvalPolicy
	>
	Vector<T> Interpolator<T, DerivPolicy, EvalPolicy>::operator()
	(
		std::span<const T> x,
		std::span<const T> xs,
		std::span<const T> ys,
		std::span<const T> dydx,
		bool doValidate
	) const
		requires HasDerivatives<DerivPolicy, T> && HasEvaluate<EvalPolicy, T>
	{
		Vector<T> y(x.size());

		(*this)
		(
			x,
			xs,
			ys,
			dydx,
			y,
			doValidate
		);

		return y;
	}

	template
	<
		std::floating_point T,
		class DerivPolicy,
		class EvalPolicy
	>
	T Interpolator<T, DerivPolicy, EvalPolicy>::operator()
	(
		T x,
		std::span<const T> xs,
		std::span<const T> ys,
		bool doValidate
	) const
		requires HasDerivatives<DerivPolicy, T> && HasEvaluate<EvalPolicy, T>
	{
		Vector<T> dydx(xs.size());
		deriv(xs, ys, dydx, doValidate);

		const std::array<T, 1> xIn{ x };
		std::array<T, 1> y{};

		(*this)
		(
			xIn,
			xs,
			ys,
			dydx,
			y,
			doValidate
		);

		return y.front();
	}

	template
	<
		std::floating_point T,
		class DerivPolicy,
		class EvalPolicy
	>
	T Interpolator<T, DerivPolicy, EvalPolicy>::operator()
	(
		T x,
		std::span<const T> xs,
		std::span<const T> ys,
		std::span<const T> dydx,
		bool doValidate
	) const
		requires HasDerivatives<DerivPolicy, T> && HasEvaluate<EvalPolicy, T>
	{
		const std::array<T, 1> xIn{ x };
		std::array<T, 1> y{};

		(*this)
		(
			xIn,
			xs,
			ys,
			dydx,
			y,
			doValidate
		);

		return y.front();
	}

} // namespace uv::math::interp


