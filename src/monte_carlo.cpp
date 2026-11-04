#include "pricing.hpp"
#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>

namespace {
// Terminal prices and payoffs under the risk-neutral measure, where
// S_T = S * exp((r - sigma^2/2) T + sigma sqrt(T) Z).
struct Model {
    const Option& o;
    double drift, vol, disc;

    explicit Model(const Option& option)
        : o(option),
          drift((o.r - 0.5 * o.sigma * o.sigma) * o.T),
          vol(o.sigma * std::sqrt(o.T)),
          disc(std::exp(-o.r * o.T)) {}

    double terminal(double z) const { return o.S * std::exp(drift + vol * z); }
    double payoff(double st) const { return o.is_call ? std::max(st - o.K, 0.0) : std::max(o.K - st, 0.0); }
};

// Tracks the mean and squared deviations of each *sample* with Welford's algorithm so we can
// report a standard error. Unlike a running sum of squares it cannot lose all its precision,
// or go negative, when the payoffs barely vary.
struct RunningStats {
    std::uint64_t n = 0;
    double mean = 0.0, sq_dev = 0.0;

    void add(double x) {
        ++n;
        const double delta = x - mean;
        mean += delta / n;
        sq_dev += delta * (x - mean);
    }
    double std_error() const { return std::sqrt(sq_dev / (n - 1) / n); }
};
}  // namespace

MCResult mc_price(const Option& o, std::uint64_t n_paths, std::uint64_t seed, bool antithetic) {
    validate_option(o);

    // With antithetic variates one sample is the average of the payoffs at Z and -Z,
    // so n_paths paths produce n_paths / 2 independent samples.
    const std::uint64_t n_samples = antithetic ? n_paths / 2 : n_paths;
    if (n_samples < 2) throw std::invalid_argument("need at least 2 samples for a standard error");

    std::mt19937_64 rng(seed);
    std::normal_distribution<double> normal(0.0, 1.0);
    const Model model(o);

    auto payoff_at = [&](double z) { return model.payoff(model.terminal(z)); };

    RunningStats stats;
    for (std::uint64_t i = 0; i < n_samples; ++i) {
        const double z = normal(rng);
        stats.add(antithetic ? 0.5 * (payoff_at(z) + payoff_at(-z)) : payoff_at(z));
    }
    return {model.disc * stats.mean, model.disc * stats.std_error()};
}

MCResult mc_price_control_variate(const Option& o, std::uint64_t n_paths, std::uint64_t seed) {
    validate_option(o);
    if (n_paths < 4) throw std::invalid_argument("need at least 4 paths for a control variate");

    std::mt19937_64 rng(seed);
    std::normal_distribution<double> normal(0.0, 1.0);
    const Model model(o);

    // A pilot run estimates the coefficient. Reusing these paths in the price would bias it,
    // because the coefficient would then depend on the samples it adjusts.
    const std::uint64_t n_pilot = std::max<std::uint64_t>(n_paths / 10, 2);
    double mean_x = 0.0, mean_st = 0.0, cov = 0.0, var_st = 0.0;
    for (std::uint64_t i = 0; i < n_pilot; ++i) {
        const double st = model.terminal(normal(rng));
        const double x = model.payoff(st);
        const double delta_x = x - mean_x, delta_st = st - mean_st;
        mean_x += delta_x / (i + 1);
        mean_st += delta_st / (i + 1);
        cov += delta_x * (st - mean_st);
        var_st += delta_st * (st - mean_st);
    }
    const double beta = var_st > 0.0 ? cov / var_st : 0.0;

    // E[S_T] is known exactly, so subtracting beta * (S_T - E[S_T]) leaves the mean unchanged.
    const double forward = o.S * std::exp(o.r * o.T);
    RunningStats stats;
    for (std::uint64_t i = n_pilot; i < n_paths; ++i) {
        const double st = model.terminal(normal(rng));
        stats.add(model.payoff(st) - beta * (st - forward));
    }
    return {model.disc * stats.mean, model.disc * stats.std_error()};
}
