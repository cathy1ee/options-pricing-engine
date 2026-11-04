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
};

struct Greeks {
    double delta;  // dV/dS
    double gamma;  // d2V/dS2
    double vega;   // dV/dsigma      (per 1.00 change in vol)
    double theta;  // dV/dt          (per year, i.e. -dV/dT)
    double rho;    // dV/dr          (per 1.00 change in rate)
};

// Throws std::invalid_argument unless S > 0, K > 0, sigma >= 0 and T >= 0. NaN is rejected.
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
