#include "pricing.hpp"
#include <cmath>

namespace {
constexpr double INV_SQRT_2PI = 0.3989422804014327;

double norm_pdf(double x) { return INV_SQRT_2PI * std::exp(-0.5 * x * x); }
double norm_cdf(double x) { return 0.5 * std::erfc(-x / std::sqrt(2.0)); }

// d1 and d2 from the Black-Scholes formula.
void d1_d2(const Option& o, double& d1, double& d2) {
    const double vol_sqrt_t = o.sigma * std::sqrt(o.T);
    d1 = (std::log(o.S / o.K) + (o.r + 0.5 * o.sigma * o.sigma) * o.T) / vol_sqrt_t;
    d2 = d1 - vol_sqrt_t;
}
}  // namespace

double bs_price(const Option& o) {
    double d1, d2;
    d1_d2(o, d1, d2);
    const double disc_k = o.K * std::exp(-o.r * o.T);  // present value of strike
    if (o.is_call) return o.S * norm_cdf(d1) - disc_k * norm_cdf(d2);
    return disc_k * norm_cdf(-d2) - o.S * norm_cdf(-d1);
}

Greeks bs_greeks(const Option& o) {
    double d1, d2;
    d1_d2(o, d1, d2);
    const double sqrt_t = std::sqrt(o.T);
    const double disc = std::exp(-o.r * o.T);
    const double pdf_d1 = norm_pdf(d1);

    Greeks g{};
    // Gamma and vega are identical for calls and puts.
    g.gamma = pdf_d1 / (o.S * o.sigma * sqrt_t);
    g.vega = o.S * pdf_d1 * sqrt_t;

    const double decay = -o.S * pdf_d1 * o.sigma / (2.0 * sqrt_t);
    if (o.is_call) {
        g.delta = norm_cdf(d1);
        g.theta = decay - o.r * o.K * disc * norm_cdf(d2);
        g.rho = o.K * o.T * disc * norm_cdf(d2);
    } else {
        g.delta = norm_cdf(d1) - 1.0;
        g.theta = decay + o.r * o.K * disc * norm_cdf(-d2);
        g.rho = -o.K * o.T * disc * norm_cdf(-d2);
    }
    return g;
}
