#include "pricing.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace {
constexpr double INV_SQRT_2PI = 0.3989422804014327;

double norm_pdf(double x) { return INV_SQRT_2PI * std::exp(-0.5 * x * x); }
double norm_cdf(double x) { return 0.5 * std::erfc(-x / std::sqrt(2.0)); }

// d1 and d2 from the Black-Scholes formula. Requires sigma * sqrt(T) > 0.
void d1_d2(const Option& o, double& d1, double& d2) {
    const double vol_sqrt_t = o.sigma * std::sqrt(o.T);
    d1 = (std::log(o.S / o.K) + (o.r - o.q + 0.5 * o.sigma * o.sigma) * o.T) / vol_sqrt_t;
    d2 = d1 - vol_sqrt_t;
}

// The limit of N(d1) and N(d2) as sigma * sqrt(T) -> 0: a step at the present value of the strike.
double step(double spot, double disc_k) { return spot > disc_k ? 1.0 : spot < disc_k ? 0.0 : 0.5; }
}  // namespace

void validate_option(const Option& o) {
    // Written as !(x > 0) rather than x <= 0 so that NaN is rejected too.
    if (!(o.S > 0.0)) throw std::invalid_argument("spot must be positive");
    if (!(o.K > 0.0)) throw std::invalid_argument("strike must be positive");
    if (!(o.sigma >= 0.0)) throw std::invalid_argument("volatility must not be negative");
    if (!(o.T >= 0.0)) throw std::invalid_argument("time to expiry must not be negative");
    if (!std::isfinite(o.q)) throw std::invalid_argument("dividend yield must be finite");
}

double bs_price(const Option& o) {
    validate_option(o);
    const double disc_k = o.K * std::exp(-o.r * o.T);  // present value of strike
    const double disc_s = o.S * std::exp(-o.q * o.T);  // spot net of dividends paid before expiry

    // With no randomness left the option is worth its discounted intrinsic value on the forward.
    if (o.sigma * std::sqrt(o.T) == 0.0) return std::max(o.is_call ? disc_s - disc_k : disc_k - disc_s, 0.0);

    double d1, d2;
    d1_d2(o, d1, d2);
    if (o.is_call) return disc_s * norm_cdf(d1) - disc_k * norm_cdf(d2);
    return disc_k * norm_cdf(-d2) - disc_s * norm_cdf(-d1);
}

Greeks bs_greeks(const Option& o) {
    validate_option(o);
    const double sqrt_t = std::sqrt(o.T);
    const double disc = std::exp(-o.r * o.T);
    const double div_disc = std::exp(-o.q * o.T);

    Greeks g{};
    double n_d1, n_d2, n_neg_d2, decay;
    if (o.sigma * sqrt_t == 0.0) {
        // Gamma stays 0 by convention, even at the kink where it diverges. Theta drops the
        // sigma-driven decay term for the same reason.
        n_d1 = n_d2 = step(o.S * div_disc, o.K * disc);
        n_neg_d2 = 1.0 - n_d2;
        decay = 0.0;
        g.vega = n_d1 == 0.5 ? o.S * div_disc * INV_SQRT_2PI * sqrt_t : 0.0;
    } else {
        double d1, d2;
        d1_d2(o, d1, d2);
        const double pdf_d1 = norm_pdf(d1);
        n_d1 = norm_cdf(d1);
        n_d2 = norm_cdf(d2);
        n_neg_d2 = norm_cdf(-d2);
        decay = -o.S * div_disc * pdf_d1 * o.sigma / (2.0 * sqrt_t);
        // Gamma and vega are identical for calls and puts.
        g.gamma = div_disc * pdf_d1 / (o.S * o.sigma * sqrt_t);
        g.vega = o.S * div_disc * pdf_d1 * sqrt_t;
    }

    if (o.is_call) {
        g.delta = div_disc * n_d1;
        g.theta = decay - o.r * o.K * disc * n_d2 + o.q * o.S * div_disc * n_d1;
        g.rho = o.K * o.T * disc * n_d2;
    } else {
        g.delta = div_disc * (n_d1 - 1.0);
        g.theta = decay + o.r * o.K * disc * n_neg_d2 - o.q * o.S * div_disc * (1.0 - n_d1);
        g.rho = -o.K * o.T * disc * n_neg_d2;
    }
    return g;
}
