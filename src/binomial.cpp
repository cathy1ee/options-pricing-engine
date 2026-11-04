#include "pricing.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

double binomial_price(const Option& o, std::uint32_t steps, Exercise exercise) {
    validate_option(o);
    if (steps == 0) throw std::invalid_argument("binomial tree needs at least 1 step");

    const bool american = exercise == Exercise::American;
    auto payoff = [&](double spot) { return std::max(o.is_call ? spot - o.K : o.K - spot, 0.0); };

    // With no randomness left exercising at time t is worth +-(S e^{-qt} - K e^{-rt}) today. Its only
    // turning point is at t_best, so the best time to exercise is now, at expiry, or then. With r or q
    // zero or equal t_best comes out infinite or NaN, and the range check drops it.
    if (o.sigma * std::sqrt(o.T) == 0.0) {
        if (!american) return bs_price(o);
        auto value_at = [&](double t) {
            const double forward_gain = o.S * std::exp(-o.q * t) - o.K * std::exp(-o.r * t);
            return std::max(o.is_call ? forward_gain : -forward_gain, 0.0);
        };
        const double t_best = std::log(o.q * o.S / (o.r * o.K)) / (o.q - o.r);
        const double ends = std::max(value_at(0.0), value_at(o.T));
        return t_best > 0.0 && t_best < o.T ? std::max(ends, value_at(t_best)) : ends;
    }

    // Cox-Ross-Rubinstein: the spot moves up by e^{sigma sqrt(dt)} or down by its inverse each step.
    const std::size_t n = steps;
    const double dt = o.T / n;
    const double jump = o.sigma * std::sqrt(dt);
    const double up = std::exp(jump), down = 1.0 / up;
    const double p = (std::exp((o.r - o.q) * dt) - down) / (up - down);  // risk-neutral up probability
    if (!(p > 0.0 && p < 1.0)) throw std::invalid_argument("too few steps for this rate and volatility");
    const double disc = std::exp(-o.r * dt);

    // spots[n + m] is the spot after m net up moves. A node's spot depends only on its net moves,
    // so one array serves the whole tree.
    std::vector<double> spots(2 * n + 1);
    for (std::size_t k = 0; k < spots.size(); ++k) spots[k] = o.S * std::exp(((double)k - (double)n) * jump);

    // Roll back from expiry. values[j] is the node with j up moves at the current step.
    std::vector<double> values(n + 1);
    for (std::size_t j = 0; j <= n; ++j) values[j] = payoff(spots[2 * j]);
    for (std::size_t i = n; i-- > 0;) {
        for (std::size_t j = 0; j <= i; ++j) {
            const double held = disc * (p * values[j + 1] + (1.0 - p) * values[j]);
            values[j] = american ? std::max(held, payoff(spots[n + 2 * j - i])) : held;
        }
    }
    return values[0];
}
