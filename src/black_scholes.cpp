#include "bachelier/black_scholes.hpp"

#include <algorithm>
#include <cmath>

#include "bachelier/normal.hpp"

namespace bachelier {

void OptionSpec::validate_normal() const {
    if (!std::isfinite(spot)) throw std::invalid_argument("spot must be finite");
    if (!std::isfinite(strike)) throw std::invalid_argument("strike must be finite");
    if (time_to_expiry < 0.0) throw std::invalid_argument("time_to_expiry must be >= 0");
    if (volatility < 0.0) throw std::invalid_argument("volatility must be >= 0");
    if (!std::isfinite(rate)) throw std::invalid_argument("rate must be finite");
    if (!std::isfinite(dividend_yield)) throw std::invalid_argument("dividend_yield must be finite");
}

void OptionSpec::validate() const {
    validate_normal();
    // Lognormal models additionally require strict positivity: the model is
    // built on log(S/K), which is undefined otherwise.
    if (!(spot > 0.0)) throw std::invalid_argument("spot must be > 0");
    if (!(strike > 0.0)) throw std::invalid_argument("strike must be > 0");
}

namespace black_scholes {
namespace {

struct Ds {
    double d1;
    double d2;
    double sqrt_t;
};

Ds compute_ds(const OptionSpec& s) {
    const double sqrt_t = std::sqrt(s.time_to_expiry);
    const double vol_sqrt_t = s.volatility * sqrt_t;
    const double d1 = (std::log(s.spot / s.strike) +
                       (s.rate - s.dividend_yield + 0.5 * s.volatility * s.volatility) * s.time_to_expiry) /
                      vol_sqrt_t;
    return {d1, d1 - vol_sqrt_t, sqrt_t};
}

/// Payoff at expiry, used when T == 0 or σ == 0 makes the d-terms degenerate.
double intrinsic(const OptionSpec& s) {
    const double df_r = std::exp(-s.rate * s.time_to_expiry);
    const double df_q = std::exp(-s.dividend_yield * s.time_to_expiry);
    const double fwd = s.spot * df_q;
    const double disc_k = s.strike * df_r;
    return s.type == OptionType::Call ? std::max(fwd - disc_k, 0.0) : std::max(disc_k - fwd, 0.0);
}

bool degenerate(const OptionSpec& s) {
    return s.time_to_expiry <= 0.0 || s.volatility <= 0.0;
}

}  // namespace

double forward(const OptionSpec& spec) {
    spec.validate();
    return spec.spot * std::exp((spec.rate - spec.dividend_yield) * spec.time_to_expiry);
}

double price(const OptionSpec& spec) {
    spec.validate();
    if (degenerate(spec)) return intrinsic(spec);

    const auto [d1, d2, sqrt_t] = compute_ds(spec);
    const double df_r = std::exp(-spec.rate * spec.time_to_expiry);
    const double df_q = std::exp(-spec.dividend_yield * spec.time_to_expiry);

    if (spec.type == OptionType::Call) {
        return spec.spot * df_q * norm_cdf(d1) - spec.strike * df_r * norm_cdf(d2);
    }
    return spec.strike * df_r * norm_cdf(-d2) - spec.spot * df_q * norm_cdf(-d1);
}

Greeks greeks(const OptionSpec& spec) {
    spec.validate();
    Greeks g{};
    if (degenerate(spec)) {
        // At expiry the option is a step function: delta is 0 or 1 (or -1 for a
        // put) and the higher-order sensitivities are undefined. Report the
        // limiting delta and leave the rest at zero rather than emit NaN.
        const bool itm = spec.type == OptionType::Call ? spec.spot > spec.strike : spec.spot < spec.strike;
        g.delta = itm ? (spec.type == OptionType::Call ? 1.0 : -1.0) : 0.0;
        return g;
    }

    const auto [d1, d2, sqrt_t] = compute_ds(spec);
    const double df_r = std::exp(-spec.rate * spec.time_to_expiry);
    const double df_q = std::exp(-spec.dividend_yield * spec.time_to_expiry);
    const double pdf_d1 = norm_pdf(d1);

    // Gamma and vega are identical for calls and puts.
    g.gamma = df_q * pdf_d1 / (spec.spot * spec.volatility * sqrt_t);
    g.vega = spec.spot * df_q * pdf_d1 * sqrt_t;

    const double theta_common = -(spec.spot * df_q * pdf_d1 * spec.volatility) / (2.0 * sqrt_t);

    if (spec.type == OptionType::Call) {
        g.delta = df_q * norm_cdf(d1);
        g.theta = theta_common - spec.rate * spec.strike * df_r * norm_cdf(d2) +
                  spec.dividend_yield * spec.spot * df_q * norm_cdf(d1);
        g.rho = spec.strike * spec.time_to_expiry * df_r * norm_cdf(d2);
    } else {
        g.delta = df_q * (norm_cdf(d1) - 1.0);
        g.theta = theta_common + spec.rate * spec.strike * df_r * norm_cdf(-d2) -
                  spec.dividend_yield * spec.spot * df_q * norm_cdf(-d1);
        g.rho = -spec.strike * spec.time_to_expiry * df_r * norm_cdf(-d2);
    }
    return g;
}

double parity_residual(const OptionSpec& spec) {
    OptionSpec call = spec;
    call.type = OptionType::Call;
    OptionSpec put = spec;
    put.type = OptionType::Put;
    return price(call) - price(put);
}

}  // namespace black_scholes
}  // namespace bachelier
