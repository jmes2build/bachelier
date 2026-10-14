// bachelier - options pricing and risk from the command line.

#include <charconv>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <map>
#include <optional>
#include <string>
#include <string_view>

#include "bachelier/binomial.hpp"
#include "bachelier/black_scholes.hpp"
#include "bachelier/implied_vol.hpp"
#include "bachelier/normal_model.hpp"

namespace {

using namespace bachelier;

void usage() {
    std::puts(R"(bachelier - options pricing and risk

USAGE
  bachelier price   [options]            Black-Scholes price and Greeks
  bachelier normal  [options]            Bachelier (normal) price and Greeks
  bachelier tree    [options] [--steps N] [--american]
                                      Binomial lattice price and Greeks
  bachelier iv      [options] --price P  Solve for implied volatility

OPTIONS
  -s, --spot        Underlying price           (required)
  -k, --strike      Strike price               (required)
  -t, --time        Years to expiry            (required)
  -r, --rate        Risk-free rate, decimal    (default 0.0)
  -v, --vol         Volatility, decimal        (required except for `iv`)
                    For `normal`, an absolute volatility in price units
  -q, --dividend    Dividend yield, decimal    (default 0.0)
      --put         Price a put                (default is a call)
      --steps N     Lattice steps              (default 512)
      --american    American exercise          (tree only)
      --price P     Target price               (iv only)

EXAMPLES
  bachelier price -s 100 -k 100 -t 1 -r 0.05 -v 0.2
  bachelier tree  -s 100 -k 95 -t 0.5 -r 0.03 -v 0.35 --put --american
  bachelier iv    -s 100 -k 100 -t 1 -r 0.05 --price 10.4506
  bachelier normal -s -15 -k -10 -t 1 -r 0.05 -v 8      (negative underlying))");
}

std::optional<double> to_double(std::string_view sv) {
    try {
        std::size_t used = 0;
        const double v = std::stod(std::string(sv), &used);
        if (used != sv.size()) return std::nullopt;
        return v;
    } catch (...) {
        return std::nullopt;
    }
}

struct Args {
    OptionSpec spec{};
    Exercise exercise{Exercise::European};
    int steps{512};
    std::optional<double> target_price;
    bool have_vol{false};
};

/// Returns nullopt and prints a message if the arguments are malformed.
std::optional<Args> parse(int argc, char** argv) {
    Args a{};
    a.spec.rate = 0.0;
    a.spec.dividend_yield = 0.0;

    const auto need_value = [&](int& i, std::string_view flag) -> std::optional<double> {
        if (i + 1 >= argc) {
            std::fprintf(stderr, "error: %.*s requires a value\n", static_cast<int>(flag.size()), flag.data());
            return std::nullopt;
        }
        auto v = to_double(argv[++i]);
        if (!v) {
            std::fprintf(stderr, "error: %.*s expects a number, got '%s'\n", static_cast<int>(flag.size()), flag.data(),
                         argv[i]);
        }
        return v;
    };

    for (int i = 2; i < argc; ++i) {
        std::string_view arg = argv[i];
        if (arg == "-s" || arg == "--spot") {
            auto v = need_value(i, arg); if (!v) return std::nullopt; a.spec.spot = *v;
        } else if (arg == "-k" || arg == "--strike") {
            auto v = need_value(i, arg); if (!v) return std::nullopt; a.spec.strike = *v;
        } else if (arg == "-t" || arg == "--time") {
            auto v = need_value(i, arg); if (!v) return std::nullopt; a.spec.time_to_expiry = *v;
        } else if (arg == "-r" || arg == "--rate") {
            auto v = need_value(i, arg); if (!v) return std::nullopt; a.spec.rate = *v;
        } else if (arg == "-v" || arg == "--vol") {
            auto v = need_value(i, arg); if (!v) return std::nullopt; a.spec.volatility = *v; a.have_vol = true;
        } else if (arg == "-q" || arg == "--dividend") {
            auto v = need_value(i, arg); if (!v) return std::nullopt; a.spec.dividend_yield = *v;
        } else if (arg == "--price") {
            auto v = need_value(i, arg); if (!v) return std::nullopt; a.target_price = *v;
        } else if (arg == "--steps") {
            auto v = need_value(i, arg); if (!v) return std::nullopt; a.steps = static_cast<int>(*v);
        } else if (arg == "--put") {
            a.spec.type = OptionType::Put;
        } else if (arg == "--american") {
            a.exercise = Exercise::American;
        } else {
            std::fprintf(stderr, "error: unknown option '%s'\n", argv[i]);
            return std::nullopt;
        }
    }
    return a;
}

void print_greeks(const Greeks& g) {
    std::printf("  delta      %12.6f\n", g.delta);
    std::printf("  gamma      %12.6f\n", g.gamma);
    std::printf("  vega       %12.6f   (%.6f per 1%% vol)\n", g.vega, g.vega / 100.0);
    std::printf("  theta      %12.6f   (%.6f per day)\n", g.theta, g.theta / 365.0);
    std::printf("  rho        %12.6f   (%.6f per 1%% rate)\n", g.rho, g.rho / 100.0);
}

void print_spec(const OptionSpec& s) {
    std::printf("  %s  S=%.4f  K=%.4f  T=%.4f  r=%.4f  q=%.4f  vol=%.4f\n", to_string(s.type).c_str(), s.spot,
                s.strike, s.time_to_expiry, s.rate, s.dividend_yield, s.volatility);
}

int cmd_price(const Args& a) {
    if (!a.have_vol) {
        std::fputs("error: --vol is required for `price`\n", stderr);
        return 2;
    }
    print_spec(a.spec);
    std::printf("\n  price      %12.6f\n", black_scholes::price(a.spec));
    std::printf("  forward    %12.6f\n\n", black_scholes::forward(a.spec));
    print_greeks(black_scholes::greeks(a.spec));
    return 0;
}

int cmd_normal(const Args& a) {
    if (!a.have_vol) {
        std::fputs("error: --vol is required for `normal`\n", stderr);
        return 2;
    }
    print_spec(a.spec);
    std::printf("  (normal model - volatility is absolute, in price units)\n");
    std::printf("\n  price      %12.6f\n", normal_model::price(a.spec));
    std::printf("  forward    %12.6f\n\n",
                a.spec.spot * std::exp((a.spec.rate - a.spec.dividend_yield) * a.spec.time_to_expiry));
    print_greeks(normal_model::greeks(a.spec));
    return 0;
}

int cmd_tree(const Args& a) {
    if (!a.have_vol) {
        std::fputs("error: --vol is required for `tree`\n", stderr);
        return 2;
    }
    print_spec(a.spec);
    std::printf("  %s exercise, %d steps\n\n", to_string(a.exercise).c_str(), a.steps);

    const double lattice = binomial::price(a.spec, a.exercise, a.steps);
    std::printf("  price      %12.6f\n", lattice);
    if (a.exercise == Exercise::European) {
        const double closed = black_scholes::price(a.spec);
        std::printf("  black-scholes %9.6f   (lattice error %.2e)\n", closed, std::abs(lattice - closed));
    } else {
        const double euro = binomial::price(a.spec, Exercise::European, a.steps);
        std::printf("  european   %12.6f\n", euro);
        std::printf("  early-exercise premium %.6f\n", lattice - euro);
    }
    std::printf("\n");
    print_greeks(binomial::greeks(a.spec, a.exercise, a.steps));
    return 0;
}

int cmd_iv(const Args& a) {
    if (!a.target_price) {
        std::fputs("error: --price is required for `iv`\n", stderr);
        return 2;
    }
    const auto bounds = implied_vol::price_bounds(a.spec);
    const auto result = implied_vol::solve(a.spec, *a.target_price);
    if (!result) {
        std::fprintf(stderr,
                     "error: price %.6f admits no implied volatility.\n"
                     "       no-arbitrage bounds for this contract are [%.6f, %.6f].\n",
                     *a.target_price, bounds.lower, bounds.upper);
        return 1;
    }
    print_spec(a.spec);
    std::printf("\n  target price  %12.6f\n", *a.target_price);
    std::printf("  implied vol   %12.6f   (%.4f%%)\n", result->volatility, result->volatility * 100.0);
    std::printf("  iterations    %12d\n", result->iterations);
    std::printf("  residual      %12.3e\n\n", result->residual);

    OptionSpec solved = a.spec;
    solved.volatility = result->volatility;
    print_greeks(black_scholes::greeks(solved));
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        usage();
        return 2;
    }
    const std::string_view cmd = argv[1];
    if (cmd == "-h" || cmd == "--help" || cmd == "help") {
        usage();
        return 0;
    }

    auto args = parse(argc, argv);
    if (!args) return 2;

    try {
        if (cmd == "price") return cmd_price(*args);
        if (cmd == "normal") return cmd_normal(*args);
        if (cmd == "tree") return cmd_tree(*args);
        if (cmd == "iv") return cmd_iv(*args);
        std::fprintf(stderr, "error: unknown command '%s'\n\n", argv[1]);
        usage();
        return 2;
    } catch (const std::exception& e) {
        // stdout is block-buffered and stderr is not, so without this flush a
        // partially-printed result appears after the error that aborted it.
        std::fflush(stdout);
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
}
