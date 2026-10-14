# bachelier

[![CI](https://github.com/jmes2build/bachelier/actions/workflows/ci.yml/badge.svg)](https://github.com/jmes2build/bachelier/actions/workflows/ci.yml)

A small options pricing and risk library in C++20, with a command-line front end.
No external dependencies — just a compiler and `make`.

Named for [Louis Bachelier](https://en.wikipedia.org/wiki/Louis_Bachelier), whose
1900 thesis *Théorie de la spéculation* gave the first mathematical treatment of
option pricing — five years before Einstein's paper on Brownian motion, and
seventy before Black, Scholes and Merton. The library implements both his
arithmetic model and the lognormal Black-Scholes-Merton one that replaced it.

```
$ bachelier price -s 100 -k 100 -t 1 -r 0.05 -v 0.2
  call  S=100.0000  K=100.0000  T=1.0000  r=0.0500  q=0.0000  vol=0.2000

  price         10.450584
  forward      105.127110

  delta          0.636831
  gamma          0.018762
  vega          37.524035   (0.375240 per 1% vol)
  theta         -6.414028   (-0.017573 per day)
  rho           53.232482   (0.532325 per 1% rate)
```

## What it does

- **Black-Scholes-Merton** pricing for European calls and puts, with a continuous dividend yield.
- **Analytic Greeks** — delta, gamma, vega, theta, rho.
- **Implied volatility** by Newton-Raphson with a bracketed bisection fallback.
- **Cox-Ross-Rubinstein binomial lattice**, including American exercise and the early-exercise premium.
- **The Bachelier (normal) model**, which prices a negative underlying or strike — the case Black-Scholes cannot express at all.

## Build

```sh
make          # build the CLI into build/bin/bachelier
make test     # build and run the test suite
make run      # build and price a sample option
```

CMake is also supported:

```sh
cmake -B build -S . && cmake --build build && ctest --test-dir build
```

## Usage

```
bachelier price   [options]                       Black-Scholes price and Greeks
bachelier tree    [options] [--steps N] [--american]   Lattice price and Greeks
bachelier iv      [options] --price P             Solve for implied volatility
```

| Flag | Meaning | Default |
| --- | --- | --- |
| `-s, --spot` | Underlying price | required |
| `-k, --strike` | Strike price | required |
| `-t, --time` | Years to expiry | required |
| `-r, --rate` | Risk-free rate (decimal) | `0.0` |
| `-v, --vol` | Volatility (decimal) | required except for `iv` |
| `-q, --dividend` | Dividend yield (decimal) | `0.0` |
| `--put` | Price a put | call |
| `--steps N` | Lattice steps | `512` |
| `--american` | American exercise (`tree` only) | European |
| `--price P` | Target price (`iv` only) | — |

An American put on a dividend-free underlying is worth more than its European
twin, because exercising early recovers the strike sooner:

```
$ bachelier tree -s 90 -k 100 -t 1 -r 0.05 -v 0.2 --put --american --steps 1000
  price         11.493351
  european      10.215080
  early-exercise premium 1.278271
```

## Using it as a library

```cmake
find_package(bachelier 0.2 REQUIRED)
target_link_libraries(your_target PRIVATE bachelier::bachelier)
```

`cmake --install` lays down the headers, the static library and a package config,
so `find_package` works from any prefix. Pulling the project in with
`add_subdirectory` or `FetchContent` builds only the library — the CLI, the tests
and the install rules switch themselves off when it is not the top-level project.

## Why two models

Black-Scholes assumes the underlying moves proportionally, so it is built on
`log(S/K)` and cannot represent a price at or below zero. That assumption broke
when rates went negative, and it never held for spread options, where the
underlying is a difference between two prices and is negative roughly half the
time.

The Bachelier model assumes absolute moves instead. Its volatility is quoted in
price units per √year rather than as a percentage, and near the money the two
are related by `σ_normal ≈ σ_lognormal · S`.

```
$ bachelier normal -s -15 -k -10 -t 1 -r 0.05 -v 8
  price          1.048904
  delta          0.235414

$ bachelier price  -s -15 -k -10 -t 1 -r 0.05 -v 8
error: spot must be > 0
```

The type system enforces the distinction: `OptionSpec::validate()` requires a
positive spot and strike, while `validate_normal()` only requires them to be
finite.

## Model

With `S` spot, `K` strike, `T` years to expiry, `r` the risk-free rate,
`q` the dividend yield and `σ` volatility:

```
d₁ = [ ln(S/K) + (r − q + σ²/2)·T ] / (σ√T)
d₂ = d₁ − σ√T

call = S·e^(−qT)·N(d₁) − K·e^(−rT)·N(d₂)
put  = K·e^(−rT)·N(−d₂) − S·e^(−qT)·N(−d₁)
```

`N` is the standard normal CDF, computed with `erfc` rather than a polynomial
approximation — it stays accurate in the tails, which matters when solving for
implied volatility on far out-of-the-money options.

The lattice uses `u = e^(σ√Δt)`, `d = 1/u`, and risk-neutral probability
`p = (e^((r−q)Δt) − d) / (u − d)`, folded backwards with an early-exercise
check at each node for American options.

## Testing

`make test` runs 539 checks. Most of them verify properties that must hold by
construction rather than comparing against hard-coded numbers:

- **Put-call parity** — `C − P = S·e^(−qT) − K·e^(−rT)` across 240 combinations
  of spot, volatility, maturity and dividend yield.
- **Greeks against finite differences** of the price function, including gamma
  as a second derivative.
- **Lattice convergence** — the binomial price approaches the closed form
  monotonically as steps increase, to within 5e-3 at 4000 steps.
- **American bounds** — the American price is never below the European one, and
  the early-exercise premium on a non-dividend-paying American call is zero
  (Merton's result).
- **Implied volatility round-trips** across strikes and volatility levels, and
  returns no solution for prices that violate the no-arbitrage bounds.
- **Normal-model parity** — `C − P = e^(−rT)·(F − K)` across 180 combinations
  including negative spots and strikes, and its analytic Greeks are checked
  against finite differences of its own price function.

### A note on the implied-vol solver

Newton-Raphson alone is not enough here. Far from the money vega underflows, so
a whole band of volatilities reprice to within any reasonable price tolerance —
converging on price would return whichever value the iteration happened to land
on. The solver therefore keeps a bracket at all times, bisects whenever vega is
too small to give a meaningful step, and terminates on the width of the bracket
in volatility rather than on the price residual alone.

## Licence

MIT — see [LICENSE](LICENSE).
