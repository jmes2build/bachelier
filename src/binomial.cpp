#include "bachelier/binomial.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace bachelier::binomial {
namespace {

double payoff(const OptionSpec& spec, double underlying) {
    return spec.type == OptionType::Call ? std::max(underlying - spec.strike, 0.0)
                                         : std::max(spec.strike - underlying, 0.0);
}

struct Lattice {
    double value{};   // price at the root
    double s_up{};    // underlying after one up move
    double s_down{};  // underlying after one down move
    double v_up{};    // option value at those nodes, for delta/gamma
    double v_down{};
    double s_uu{};
    double s_ud{};
    double s_dd{};
    double v_uu{};
    double v_ud{};
    double v_dd{};
};

Lattice build(const OptionSpec& spec, Exercise exercise, int steps) {
    spec.validate();
    if (steps < 1) throw std::invalid_argument("steps must be >= 1");

    if (spec.time_to_expiry <= 0.0 || spec.volatility <= 0.0) {
        Lattice flat{};
        flat.value = payoff(spec, spec.spot);
        return flat;
    }

    const double dt = spec.time_to_expiry / steps;
    const double u = std::exp(spec.volatility * std::sqrt(dt));
    const double d = 1.0 / u;
    const double growth = std::exp((spec.rate - spec.dividend_yield) * dt);
    const double p = (growth - d) / (u - d);
    const double disc = std::exp(-spec.rate * dt);

    // A risk-neutral probability outside [0,1] means the lattice is arbitrageable,
    // which happens when steps is too small for the drift relative to vol.
    if (!(p >= 0.0 && p <= 1.0)) {
        throw std::invalid_argument(
            "risk-neutral probability outside [0,1]; increase steps or reduce the rate/vol ratio");
    }

    // Terminal underlying prices, then fold backwards.
    std::vector<double> values(static_cast<std::size_t>(steps) + 1);
    for (int i = 0; i <= steps; ++i) {
        const double s_t = spec.spot * std::pow(u, 2 * i - steps);
        values[static_cast<std::size_t>(i)] = payoff(spec, s_t);
    }

    Lattice out{};
    for (int step = steps - 1; step >= 0; --step) {
        for (int i = 0; i <= step; ++i) {
            const auto idx = static_cast<std::size_t>(i);
            double cont = disc * (p * values[idx + 1] + (1.0 - p) * values[idx]);
            if (exercise == Exercise::American) {
                const double s_t = spec.spot * std::pow(u, 2 * i - step);
                cont = std::max(cont, payoff(spec, s_t));
            }
            values[idx] = cont;
        }
        // Capture the two levels nearest the root on the way past, so delta and
        // gamma come from the lattice itself rather than a re-priced bump.
        if (step == 2) {
            out.s_uu = spec.spot * u * u;
            out.s_ud = spec.spot;
            out.s_dd = spec.spot * d * d;
            out.v_uu = values[2];
            out.v_ud = values[1];
            out.v_dd = values[0];
        }
        if (step == 1) {
            out.s_up = spec.spot * u;
            out.s_down = spec.spot * d;
            out.v_up = values[1];
            out.v_down = values[0];
        }
    }
    out.value = values[0];
    return out;
}

}  // namespace

double price(const OptionSpec& spec, Exercise exercise, int steps) {
    return build(spec, exercise, steps).value;
}

Greeks greeks(const OptionSpec& spec, Exercise exercise, int steps) {
    const Lattice lat = build(spec, exercise, steps);
    Greeks g{};

    if (lat.s_up > 0.0 && lat.s_up != lat.s_down) {
        g.delta = (lat.v_up - lat.v_down) / (lat.s_up - lat.s_down);
    }
    if (lat.s_uu > 0.0 && lat.s_uu != lat.s_dd) {
        const double d_up = (lat.v_uu - lat.v_ud) / (lat.s_uu - lat.s_ud);
        const double d_down = (lat.v_ud - lat.v_dd) / (lat.s_ud - lat.s_dd);
        g.gamma = (d_up - d_down) / (0.5 * (lat.s_uu - lat.s_dd));
    }

    // Central differences for the rest. The bumps are absolute and chosen to be
    // large enough to dominate lattice discretisation noise.
    const auto bump = [&](auto mutate, double h) {
        OptionSpec up = spec, down = spec;
        mutate(up, h);
        mutate(down, -h);
        return (price(up, exercise, steps) - price(down, exercise, steps)) / (2.0 * h);
    };

    g.vega = bump([](OptionSpec& s, double h) { s.volatility += h; }, 0.01);
    g.rho = bump([](OptionSpec& s, double h) { s.rate += h; }, 1e-4);

    // theta is -dV/dT. Central, like vega and rho above - the previous
    // one-sided difference was only first-order accurate while everything
    // else here was second-order.
    constexpr double theta_h = 1e-3;
    if (spec.time_to_expiry > 2.0 * theta_h) {
        g.theta = -bump([](OptionSpec& s, double h) { s.time_to_expiry += h; }, theta_h);
    }
    return g;
}

double early_exercise_premium(const OptionSpec& spec, int steps) {
    return price(spec, Exercise::American, steps) - price(spec, Exercise::European, steps);
}

}  // namespace bachelier::binomial
