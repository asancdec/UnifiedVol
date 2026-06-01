// SPDX-License-Identifier: Apache-2.0

#pragma once

extern "C"
{
    double implied_volatility_from_a_transformed_rational_guess(
        double price,
        double forward,
        double strike,
        double time_to_expiry,
        int is_call
    );
}