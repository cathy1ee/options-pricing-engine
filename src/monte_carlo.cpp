#include "pricing.hpp"
#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>

MCResult mc_price(const Option& o, std::uint64_t n_paths, std::uint64_t seed, bool antithetic) {
    validate_option(o);

    // With antithetic variates one sample is the average of the payoffs at Z and -Z,
    // so n_paths paths produce n_paths / 2 independent samples.
    const std::uint64_t n_samples = antithetic ? n_paths / 2 : n_paths;
    if (n_samples < 2) throw std::invalid_argument("need at least 2 samples for a standard error");

    std::mt19937_64 rng(seed);
    std::normal_distribution<double> normal(0.0, 1.0);

    // Under the risk-neutral measure, S_T = S * exp((r - sigma^2/2) T + sigma sqrt(T) Z).
    const double drift = (o.r - 0.5 * o.sigma * o.sigma) * o.T;
    const double vol = o.sigma * std::sqrt(o.T);
    const double disc = std::exp(-o.r * o.T);

    auto payoff = [&](double z) {
        const double st = o.S * std::exp(drift + vol * z);
        return o.is_call ? std::max(st - o.K, 0.0) : std::max(o.K - st, 0.0);
    };

    // Track the mean and squared deviations of each *sample* with Welford's algorithm so we can
    // report a standard error. Unlike a running sum of squares it cannot lose all its precision,
    // or go negative, when the payoffs barely vary.
    double mean = 0.0, sq_dev = 0.0;
    for (std::uint64_t i = 0; i < n_samples; ++i) {
        const double z = normal(rng);
        const double x = antithetic ? 0.5 * (payoff(z) + payoff(-z)) : payoff(z);
        const double delta = x - mean;
        mean += delta / (i + 1);
        sq_dev += delta * (x - mean);
    }

    const double var = sq_dev / (n_samples - 1);
    return {disc * mean, disc * std::sqrt(var / n_samples)};
}
