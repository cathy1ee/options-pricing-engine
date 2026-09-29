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

// ---- Closed-form Black-Scholes ----
double bs_price(const Option& o);
Greeks bs_greeks(const Option& o);

// ---- Monte Carlo ----
struct MCResult {
    double price;
    double std_error;  // standard error of the estimate
};

// Simulates n_paths terminal prices under geometric Brownian motion.
// antithetic = true pairs each normal draw Z with -Z to reduce variance.
MCResult mc_price(const Option& o, std::uint64_t n_paths, std::uint64_t seed, bool antithetic);
