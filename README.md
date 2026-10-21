# Options Pricing Engine (C++)

Prices European options three ways and checks each one against the others: closed-form Black-Scholes, analytic Greeks, and Monte Carlo simulation with antithetic variance reduction. It has no dependencies beyond the C++17 standard library.

## Features

| Component | What it does |
|---|---|
| `bs_price` | Closed-form Black-Scholes price for European calls and puts. |
| `bs_greeks` | Analytic delta, gamma, vega, theta and rho. |
| `mc_price` | Monte Carlo pricing under risk-neutral geometric Brownian motion. It returns the price **and its standard error**, with optional antithetic variates. |

## Model

- **Underlying:** geometric Brownian motion, `S_T = S·exp((r − σ²/2)T + σ√T·Z)` with `Z ~ N(0,1)`.
- **Black-Scholes:** `C = S·N(d1) − K·e^(−rT)·N(d2)`, where `d1 = [ln(S/K) + (r + σ²/2)T] / (σ√T)` and `d2 = d1 − σ√T`.
- **Antithetic variates:** each draw `Z` is paired with `−Z`, and the two payoffs are averaged into one sample. The standard error is computed over these pairs, since they are the independent samples (N paths give N/2 samples).
- **Units:** vega and rho are per 1.00 change in volatility or rate. Theta is per year.

## Build and run

```bash
make test     # correctness suite
make bench    # convergence and speed benchmark
make          # both
```

Requires g++ or clang++ with C++17. The build uses `-Wall -Wextra` and produces no warnings.

## Validation

`make test` runs 1,562 checks across 240 option configurations, and all of them pass. The grid covers strikes from 80 to 120 (in, at and out of the money), volatilities of 10%, 25% and 50%, expiries from 0.1 to 2 years, rates of 0% and 5%, and both calls and puts.

| Check | What it verifies | Tolerance |
|---|---|---|
| Textbook values | Hull's example: S=42, K=40, r=10%, σ=20%, T=0.5 gives a call of 4.76 and a put of 0.81 | ±0.005 |
| Put-call parity | `C − P = S − K·e^(−rT)` on every configuration | 1e-10 |
| Greeks | Analytic Greeks against central finite differences of the price | 1e-4 relative (max error observed: 1.35e-5) |
| Monte Carlo | MC price (200k paths) against Black-Scholes | within 4 standard errors |

## Results

From `make bench`, for an at-the-money 1-year call (S=K=100, r=5%, σ=20%, Black-Scholes price 10.4506). The error is the RMSE over 20 seeds.

| Paths | Method | RMSE vs BS | Error % of price | ms / run |
|---:|---|---:|---:|---:|
| 1,000 | plain | 0.3701 | 3.54% | 0.04 |
| 1,000 | antithetic | 0.2637 | 2.52% | 0.03 |
| 10,000 | plain | 0.1512 | 1.45% | 0.43 |
| 10,000 | antithetic | 0.1109 | 1.06% | 0.26 |
| 100,000 | plain | 0.0483 | 0.46% | 4.6 |
| 100,000 | antithetic | 0.0350 | 0.34% | 3.3 |
| 1,000,000 | plain | 0.0117 | 0.11% | 44 |
| 1,000,000 | antithetic | 0.0100 | 0.10% | 25 |
| 10,000,000 | plain | 0.0063 | 0.06% | 430 |
| 10,000,000 | antithetic | 0.0032 | 0.03% | 258 |

- **Convergence rate:** the fitted slope of log(RMSE) against log(N) is −0.465, close to the theoretical −0.5 (error shrinks like 1/√N). The gap comes from using only 20 seeds per point.
- **Variance reduction:** at 1M paths, antithetic variates cut the estimator's variance by 50% (standard error 0.0147 → 0.0104). Reaching the same accuracy without them would take about twice as many paths.
- **Speed:** antithetic runs are also about 40% faster at the same path count, since each normal draw is used twice.
- **Closed form:** 18.6M prices per second (about 54 ns each), single-threaded.

Timings come from a single Linux machine with `-O2`; your numbers will differ.

## Layout

```
include/pricing.hpp     public API: Option, Greeks, MCResult
src/black_scholes.cpp   closed-form price and Greeks
src/monte_carlo.cpp     Monte Carlo pricer with standard error and antithetic variates
tests/test_pricing.cpp  self-contained test suite (no framework)
bench/convergence.cpp   convergence, variance-reduction and speed benchmark
Makefile
```

## Limitations and next steps

- Only European exercise, constant volatility and no dividends.
- **Implied volatility:** invert Black-Scholes with Newton's method, using the analytic vega, and fall back to bisection when vega is small.
- **American options:** a binomial tree (Cox-Ross-Rubinstein), checked against Black-Scholes for European options as a consistency test.
- **More variance reduction:** use the discounted underlying as a control variate, or quasi-random (Sobol) draws.
- **Parallel Monte Carlo:** independent RNG streams per thread.

Reference: Hull, J. C. *Options, Futures, and Other Derivatives.*
