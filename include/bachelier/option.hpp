#pragma once

#include <stdexcept>
#include <string>

namespace bachelier {

enum class OptionType { Call, Put };

enum class Exercise { European, American };

/// Contract and market inputs for a vanilla equity option.
///
/// Rates and volatility are continuously compounded annualised decimals,
/// e.g. 5% is 0.05. `time_to_expiry` is in years.
struct OptionSpec {
    double spot{};            // S  - underlying price
    double strike{};          // K  - strike price
    double time_to_expiry{};  // T  - years to expiry
    double rate{};            // r  - risk-free rate
    double volatility{};      // σ  - annualised volatility
    double dividend_yield{};  // q  - continuous dividend yield
    OptionType type{OptionType::Call};

    /// Throws std::invalid_argument if the inputs cannot be priced under a
    /// lognormal model, which requires a strictly positive spot and strike.
    void validate() const;

    /// Validation for the normal (Bachelier) model, which permits a
    /// non-positive spot and strike - that is the whole reason to reach for
    /// it. Everything else is checked as in validate().
    void validate_normal() const;
};

/// First- and second-order sensitivities.
///
/// `vega`, `theta` and `rho` are reported in raw (per unit) terms:
/// divide vega and rho by 100 for a 1% move, and theta by 365 for one day.
struct Greeks {
    double delta{};
    double gamma{};
    double vega{};
    double theta{};
    double rho{};
};

inline std::string to_string(OptionType t) {
    return t == OptionType::Call ? "call" : "put";
}

inline std::string to_string(Exercise e) {
    return e == Exercise::European ? "european" : "american";
}

}  // namespace bachelier
