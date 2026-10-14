#include "bachelier/normal_model.hpp"

#include <algorithm>
#include <cmath>

#include "bachelier/normal.hpp"

namespace bachelier::normal_model {
namespace {

struct Terms {
    double df;       // e^(-rT), discount factor
    double fwd;      // F = S·e^((r-q)T)
    double v;        // sigma·sqrt(T), the absolute move scale over the life
    double d;        // (F - K) / v
    double n_d;      // phi(d)
    double N_d;      // Phi(d)
};

Terms terms(const OptionSpec& s) {
    Terms t{};
    t.df = std::exp(-s.rate * s.time_to_expiry);
    t.fwd = s.spot * std::exp((s.rate - s.dividend_yield) * s.time_to_expiry);
    t.v = s.volatility * std::sqrt(s.time_to_expiry);
    t.d = (t.fwd - s.strike) / t.v;
    t.n_d = norm_pdf(t.d);
    t.N_d = norm_cdf(t.d);
    return t;
}

bool degenerate(const OptionSpec& s) {
    return s.time_to_expiry <= 0.0 || s.volatility <= 0.0;
}

/// Discounted intrinsic on the forward, used when the diffusion term vanishes.
double intrinsic(const OptionSpec& s) {
    const double df = std::exp(-s.rate * s.time_to_expiry);
    const double fwd = s.spot * std::exp((s.rate - s.dividend_yield) * s.time_to_expiry);
    const double sign = s.type == OptionType::Call ? 1.0 : -1.0;
    return df * std::max(sign * (fwd - s.strike), 0.0);
}

}  // namespace

double price(const OptionSpec& spec) {
    spec.validate_normal();
    if (degenerate(spec)) return intrinsic(spec);

    const Terms t = terms(spec);
    // Undiscounted value on the forward. The call and put forms differ only in
    // which tail of the normal is integrated; the phi(d) term is shared.
    const double undiscounted = spec.type == OptionType::Call
                                    ? (t.fwd - spec.strike) * t.N_d + t.v * t.n_d
                                    : (spec.strike - t.fwd) * norm_cdf(-t.d) + t.v * t.n_d;
    return t.df * undiscounted;
}

Greeks greeks(const OptionSpec& spec) {
    spec.validate_normal();
    Greeks g{};
    if (degenerate(spec)) {
        const double fwd = spec.spot * std::exp((spec.rate - spec.dividend_yield) * spec.time_to_expiry);
        const bool itm = spec.type == OptionType::Call ? fwd > spec.strike : fwd < spec.strike;
        g.delta = itm ? (spec.type == OptionType::Call ? 1.0 : -1.0) : 0.0;
        return g;
    }

    const Terms t = terms(spec);
    const double T = spec.time_to_expiry;
    const double carry = std::exp((spec.rate - spec.dividend_yield) * T);  // dF/dS

    // Gamma and vega do not depend on the option type: both derive from the
    // density at d, which is symmetric between the two payoffs.
    g.gamma = t.df * (t.n_d / t.v) * carry * carry;
    g.vega = t.df * t.n_d * std::sqrt(T);

    // Undiscounted value and its derivative in F. The phi(d) terms cancel in
    // dU/dF, leaving just the probability of finishing in the money.
    const double u = spec.type == OptionType::Call
                         ? (t.fwd - spec.strike) * t.N_d + t.v * t.n_d
                         : (spec.strike - t.fwd) * norm_cdf(-t.d) + t.v * t.n_d;
    const double du_df = spec.type == OptionType::Call ? t.N_d : -norm_cdf(-t.d);

    g.delta = t.df * du_df * carry;

    // theta is -dV/dT, matching the Black-Scholes convention in this library.
    // T enters through the discount factor, the forward, and the move scale.
    const double du_dT = du_df * (spec.rate - spec.dividend_yield) * t.fwd +
                         t.n_d * spec.volatility / (2.0 * std::sqrt(T));
    g.theta = spec.rate * t.df * u - t.df * du_dT;

    // rho: r moves both the discount factor and the forward.
    g.rho = T * t.df * (du_df * t.fwd - u);
    return g;
}

double parity_residual(const OptionSpec& spec) {
    OptionSpec call = spec;
    call.type = OptionType::Call;
    OptionSpec put = spec;
    put.type = OptionType::Put;
    return price(call) - price(put);
}

}  // namespace bachelier::normal_model
