#pragma once

#include "bachelier/option.hpp"

/// The Bachelier (normal, or arithmetic Brownian motion) model.
///
/// Louis Bachelier's 1900 model assumes the underlying moves by absolute
/// amounts rather than proportional ones. Black-Scholes replaced it because
/// prices cannot go negative - but that assumption is exactly what breaks in
/// negative-rate markets and in spread and some commodity options, where the
/// underlying genuinely can be at or below zero. The normal model is the
/// standard choice there, and unlike Black-Scholes it prices a non-positive
/// forward or strike without complaint.
///
/// Volatility here is an *absolute* (normal) volatility quoted in price units
/// per sqrt(year), not the proportional volatility Black-Scholes takes. Near
/// the money the two are related by sigma_normal ~ sigma_lognormal * S.
///
/// Note: `normal.hpp` holds the normal *distribution* helpers; this header is
/// the normal *pricing model*.
namespace bachelier::normal_model {

/// Undiscounted-forward price of a European option under the normal model.
double price(const OptionSpec& spec);

/// Analytic sensitivities. Delta and gamma are with respect to spot, vega with
/// respect to the absolute volatility.
Greeks greeks(const OptionSpec& spec);

/// C - P under the normal model, which must equal e^(-rT)·(F - K).
double parity_residual(const OptionSpec& spec);

}  // namespace bachelier::normal_model
