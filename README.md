# Options Pricing Engine (C++)

![tests](https://github.com/cathy1ee/options-pricing-engine/actions/workflows/test.yml/badge.svg)

Prices European options three ways and checks each one against the others: closed-form Black-Scholes, analytic Greeks, and Monte Carlo simulation with antithetic variance reduction. It has no dependencies beyond the C++17 standard library.

## Features

| Component | What it does |
|---|---|
| `bs_price` | Closed-form Black-Scholes price for European calls and puts. |
| `bs_greeks` | Analytic delta, gamma, vega, theta and rho. |
| `implied_vol` | Recovers volatility from a market price, using Newton's method with the analytic vega and falling back to bisection. |
| `mc_price_parallel` | `mc_price` across several threads, each with its own random stream. |
| `mc_price_sobol` | Quasi-Monte Carlo pricing with scrambled Sobol points, which converges faster than random draws. |
| `mc_price_control_variate` | Monte Carlo pricing with the discounted underlying as a control variate. |
| `binomial_price` | European and American prices on a Cox-Ross-Rubinstein binomial tree. |
| `mc_price` | Monte Carlo pricing under risk-neutral geometric Brownian motion. It returns the price **and its standard error**, with optional antithetic variates. |

## Model

- **Underlying:** geometric Brownian motion with continuous dividend yield `q`, `S_T = S·exp((r − q − σ²/2)T + σ√T·Z)` with `Z ~ N(0,1)`. Every pricer takes `q` from `Option`, where it defaults to 0.
- **Black-Scholes:** `C = S·e^(−qT)·N(d1) − K·e^(−rT)·N(d2)`, where `d1 = [ln(S/K) + (r − q + σ²/2)T] / (σ√T)` and `d2 = d1 − σ√T`.
- **Binomial tree:** the spot moves up by `u = e^(σ√Δt)` or down by `1/u` each step, with risk-neutral up probability `p = (e^((r−q)Δt) − 1/u) / (u − 1/u)`. American prices take the larger of holding and exercising at each node.
- **Antithetic variates:** each draw `Z` is paired with `−Z`, and the two payoffs are averaged into one sample. The standard error is computed over these pairs, since they are the independent samples (N paths give N/2 samples).
- **Control variate:** each sample becomes `payoff − β(S_T − S·e^((r−q)T))`, which keeps the mean because `E[S_T] = S·e^((r−q)T)`. The first 10% of the paths are a pilot run that sets `β` and are not used in the price, which keeps it unbiased.
- **Sobol points:** each path needs one normal draw, so only the first Sobol dimension is used, and draws come from an inverse normal CDF. The paths are split into 32 independently scrambled replications. Sobol points are not independent, so the standard error comes from the spread of the replications' estimates.
- **Parallel Monte Carlo:** samples are split evenly across threads, each seeded from the main seed with SplitMix64, and the threads' results are merged in a fixed order. A seed and thread count always give the same result, but changing the thread count changes it.
- **Units:** vega and rho are per 1.00 change in volatility or rate. Theta is per year.
- **Edge cases:** `T = 0` and `σ = 0` return the discounted intrinsic value on the forward, with gamma 0 by convention. Invalid inputs (`S ≤ 0`, `K ≤ 0`, `T < 0`, `σ < 0`, NaN) and Monte Carlo runs with fewer than 2 samples throw `std::invalid_argument`.

## Build and run

```bash
make test     # correctness suite
make bench    # convergence and speed benchmark
make          # both
```

Requires g++ or clang++ with C++17. The build uses `-Wall -Wextra` and produces no warnings.

## Validation

`make test` runs 7,258 checks, and all of them pass. Most run across a grid of 480 option configurations, which covers strikes from 80 to 120 (in, at and out of the money), volatilities of 10%, 25% and 50%, expiries from 0.1 to 2 years, rates of 0% and 5%, dividend yields of 0% and 3%, and both calls and puts.

| Check | What it verifies | Tolerance |
|---|---|---|
| Textbook values | Hull's example: S=42, K=40, r=10%, σ=20%, T=0.5 gives a call of 4.76 and a put of 0.81 | ±0.005 |
| Put-call parity | `C − P = S·e^(−qT) − K·e^(−rT)` on every configuration | 1e-10 |
| Dividends | A yield matches a lower spot with no dividends, American calls can beat European calls, and the zero-volatility exercise time matches a brute-force search | 1e-10, 1e-6 |
| Greeks | Analytic Greeks against central finite differences of the price | 1e-4 relative plus 5e-6 absolute (max absolute error observed: 1.35e-5) |
| Monte Carlo | Plain and antithetic MC prices (200k paths) against Black-Scholes | within 4 standard errors |
| Control variate | MC price (200k paths) against Black-Scholes, and a standard error that matches the spread over 100 seeds and beats plain MC | within 4 standard errors, 25% for spread |
| Sobol | MC price (65,536 paths) against Black-Scholes, and a standard error that matches the spread over 100 seeds, beats plain MC, and shrinks faster than 1/√N | within 5 standard errors, 25% for spread |
| Parallel | MC price (200k paths, 4 threads) against Black-Scholes, repeatable for a seed, and matching a single-threaded price and standard error | within 4 standard errors, 5% on standard error |
| Standard error | Matches the spread of prices over 100 seeds, shrinks like 1/√N, and is smaller with antithetic variates | 25% for spread, 10% for scaling |
| Implied volatility | Pricing and inverting recovers σ on every configuration, or the price where vega is too small to pin down σ | 1e-8 on σ, or 1e-12 on price |
| Binomial tree | European prices (500 steps) against Black-Scholes, American calls equal European calls without dividends, American prices at or above the European price and intrinsic value, and Hull's American put at 5 to 500 steps | ±0.02, 1e-9, 1e-12, ±0.002 |
| Limits | Prices and Greeks at `T = 0` and `σ = 0` against their closed-form limits | 1e-12 |
| Invalid inputs | Every pricer rejects bad options, and Monte Carlo rejects too few paths | exact |
| Low volatility | MC standard error at `σ = 0` and `σ = 1e-8` | exact at 0, 5% at 1e-8 |

## Results

From `make bench`, for an at-the-money 1-year call (S=K=100, r=5%, σ=20%, Black-Scholes price 10.4506). The error is the RMSE over 20 seeds.

| Paths | Method | RMSE vs BS | Error % of price | ms / run |
|---:|---|---:|---:|---:|
| 1,000 | plain | 0.3701 | 3.54% | 0.04 |
| 1,000 | antithetic | 0.2637 | 2.52% | 0.03 |
| 10,000 | plain | 0.1512 | 1.45% | 0.40 |
| 10,000 | antithetic | 0.1109 | 1.06% | 0.26 |
| 100,000 | plain | 0.0483 | 0.46% | 4.1 |
| 100,000 | antithetic | 0.0350 | 0.34% | 2.6 |
| 1,000,000 | plain | 0.0117 | 0.11% | 41 |
| 1,000,000 | antithetic | 0.0100 | 0.10% | 26 |
| 10,000,000 | plain | 0.0063 | 0.06% | 403 |
| 10,000,000 | antithetic | 0.0032 | 0.03% | 260 |

- **Convergence rate:** the fitted slope of log(RMSE) against log(N) is −0.465, close to the theoretical −0.5 (error shrinks like 1/√N). The gap comes from using only 20 seeds per point.
- **Variance reduction:** at 1M paths, antithetic variates cut the estimator's variance by 50% (standard error 0.0147 → 0.0104). Reaching the same accuracy without them would take about twice as many paths.
- **Control variate:** at 1M paths it cuts variance by 84% (standard error 0.0147 → 0.0059), even after spending 10% of the paths on the pilot. Deep out of the money the underlying says little about the payoff, and the standard error can be up to about 5% worse than plain.
- **Sobol:** at 2^20 (about 1M) paths the standard error is 0.00022, about 67 times smaller than plain Monte Carlo at 1M paths, and it shrinks faster than 1/√N.
- **Speed:** antithetic runs are also about 35% faster at the same path count, since each normal draw is used twice.
- **Parallel:** 10M antithetic paths take 260 ms on 1 thread, 130 ms on 2 (2.0x) and 78 ms on 4 (3.3x). 8 threads were no faster than 4 on this machine.
- **Closed form:** 17.7M prices per second (about 57 ns each), single-threaded.

Timings come from one run on a GitHub Actions `ubuntu-latest` runner with g++ and `-O2`; your numbers will differ.

## Layout

```
include/pricing.hpp     public API: Option, Greeks, MCResult
src/binomial.cpp        Cox-Ross-Rubinstein binomial tree
src/black_scholes.cpp   closed-form price and Greeks
src/implied_vol.cpp     implied volatility solver
src/monte_carlo.cpp     Monte Carlo pricer with standard error and antithetic variates
tests/test_pricing.cpp  self-contained test suite (no framework)
bench/convergence.cpp   convergence, variance-reduction and speed benchmark
Makefile
```

## Limitations and next steps

- Constant volatility, and dividends only as a continuous yield. American exercise is only available on the binomial tree.

Reference: Hull, J. C. *Options, Futures, and Other Derivatives.*
