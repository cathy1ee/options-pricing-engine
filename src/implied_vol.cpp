#include "pricing.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

double implied_vol(const Option& o, double market_price) {
    Option x = o;  // the solver's working copy: o.sigma is ignored
    x.sigma = 0.0;
    validate_option(x);
    if (!(o.T > 0.0)) throw std::invalid_argument("implied volatility needs a positive time to expiry");

    // No-arbitrage bounds. At sigma = 0 the price is the lower bound, and it rises towards the
    // upper bound as sigma grows. A model price can round to just below the lower bound.
    const double lower = bs_price(x);
    const double upper = o.is_call ? o.S * std::exp(-o.q * o.T) : o.K * std::exp(-o.r * o.T);
    const double slack = 1e-12 * (o.S + o.K);
    if (!(market_price >= lower - slack && market_price < upper))
        throw std::invalid_argument("market price is outside the no-arbitrage bounds");
    if (market_price <= lower) return 0.0;

    auto error_at = [&](double sigma) {
        x.sigma = sigma;
        return bs_price(x) - market_price;
    };

    // Bracket the root. The price increases with sigma, so the error changes sign once.
    double lo = 0.0, hi = 1.0;
    while (error_at(hi) < 0.0) {
        lo = hi;
        hi *= 2.0;
        if (hi > 1e6) throw std::invalid_argument("market price needs an implied volatility above 1e6");
    }

    // Start at the Manaster-Koehler guess, where vega is largest for this moneyness.
    const double moneyness = std::fabs(std::log(o.S / o.K) + (o.r - o.q) * o.T);
    double sigma = moneyness > 0.0 ? std::min(std::max(std::sqrt(2.0 * moneyness / o.T), 0.1), hi) : 0.2;
    if (!(lo < sigma && sigma < hi)) sigma = 0.5 * (lo + hi);

    // Newton's method with the analytic vega, bisecting whenever a step leaves the bracket. Deep out of
    // the money vega can be too small for a useful step, or 0, which also lands in the fallback.
    for (int i = 0; i < 100; ++i) {
        const double err = error_at(sigma);
        if (err == 0.0) break;
        (err > 0.0 ? hi : lo) = sigma;

        double next = sigma - err / bs_greeks(x).vega;
        if (!(lo < next && next < hi)) next = 0.5 * (lo + hi);

        const bool converged = std::fabs(next - sigma) < 1e-14 * std::max(sigma, 1.0);
        sigma = next;
        if (converged) break;
    }
    return sigma;
}
