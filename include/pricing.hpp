#pragma once
#include <cstdint>

// Inputs for a European option.
struct Option {
    double S;      // spot price
    double K;      // strike
    double r;      // risk-free rate (continuously compounded, annual)
    double sigma;  // volatility (annual)
    double T;      // time to expiry in years
    bool is_call;
    double q = 0.0;  // continuous dividend yield (annual)
};

struct Greeks {
    double delta;  // dV/dS
    double gamma;  // d2V/dS2
    double vega;   // dV/dsigma      (per 1.00 change in vol)
    double theta;  // dV/dt          (per year, i.e. -dV/dT)
    double rho;    // dV/dr          (per 1.00 change in rate)
};

// Throws std::invalid_argument unless S > 0, K > 0, sigma >= 0, T >= 0 and q is finite. NaN is rejected.
void validate_option(const Option& o);

// ---- Closed-form Black-Scholes ----
// Both functions validate their input. T = 0 and sigma = 0 are valid: the price is the discounted
// intrinsic value on the forward and delta is a step (+-0.5 exactly at the present value of the
// strike). Gamma is 0 by convention, and theta leaves out the sigma-driven decay term.
double bs_price(const Option& o);
Greeks bs_greeks(const Option& o);

// ---- Implied volatility ----
// Returns the volatility at which bs_price matches market_price, ignoring o.sigma. A price at the
// lower no-arbitrage bound returns 0. Throws std::invalid_argument if T = 0 or the price is outside
// the no-arbitrage bounds. Deep out of the money the price may barely depend on volatility, so the
// result then matches the price rather than any particular volatility.
double implied_vol(const Option& o, double market_price);

// ---- Binomial tree ----
enum class Exercise { European, American };

// Prices an option on a Cox-Ross-Rubinstein tree with the given number of steps. Throws
// std::invalid_argument for invalid options, for steps = 0, or when there are too few steps for
// the up probability to lie between 0 and 1. T = 0 and sigma = 0 are valid, as for bs_price.
double binomial_price(const Option& o, std::uint32_t steps, Exercise exercise);

// ---- Monte Carlo ----
struct MCResult {
    double price;
    double std_error;  // standard error of the estimate
};

// Simulates n_paths terminal prices under geometric Brownian motion.
// antithetic = true pairs each normal draw Z with -Z to reduce variance. An odd n_paths is then
// rounded down to an even number. Throws std::invalid_argument for invalid options, or for fewer
// than 2 samples (n_paths < 2, or n_paths < 4 with antithetic), where a standard error is undefined.
MCResult mc_price(const Option& o, std::uint64_t n_paths, std::uint64_t seed, bool antithetic);

// mc_price across n_threads threads, each with its own random stream. The result is reproducible
// for a fixed seed and thread count, but differs between thread counts, and from mc_price, because
// the streams differ. Throws std::invalid_argument for invalid options, n_threads = 0, or fewer
// than 2 samples.
MCResult mc_price_parallel(const Option& o, std::uint64_t n_paths, std::uint64_t seed, bool antithetic,
                           unsigned n_threads);

// Quasi-Monte Carlo with Sobol points, which fill the distribution more evenly than random draws
// and so converge faster than 1/sqrt(N). The paths are split into 32 independently scrambled
// replications, and the standard error comes from the spread of their estimates. A power of 2
// paths gives the most even coverage, and any remainder after dividing by 32 is dropped. Throws
// std::invalid_argument for invalid options or n_paths < 32.
MCResult mc_price_sobol(const Option& o, std::uint64_t n_paths, std::uint64_t seed);

// Monte Carlo with the discounted underlying as a control variate. The first 10% of the paths are a
// pilot run that sets the coefficient and do not count towards the price, which keeps it unbiased.
// Deep out of the money the underlying says little about the payoff, so the pilot can make the
// standard error slightly larger than plain Monte Carlo. Throws std::invalid_argument for invalid
// options or n_paths < 4.
MCResult mc_price_control_variate(const Option& o, std::uint64_t n_paths, std::uint64_t seed);
