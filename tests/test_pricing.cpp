// Self-contained tests: no framework needed. Build with `make test`.
#include "pricing.hpp"
#include <cmath>
#include <cstdio>
#include <vector>

static int n_checks = 0, n_failed = 0;

static void check(bool ok, const char* what, const Option& o, double got, double want) {
    ++n_checks;
    if (!ok) {
        ++n_failed;
        std::printf("FAIL %-18s S=%.0f K=%.0f r=%.2f vol=%.2f T=%.2f %s: got %.8f want %.8f\n", what, o.S,
                    o.K, o.r, o.sigma, o.T, o.is_call ? "call" : "put", got, want);
    }
}

// Grid of market inputs covering in/at/out of the money, low/high vol, short/long expiry.
static std::vector<Option> test_grid() {
    std::vector<Option> grid;
    for (double K : {80.0, 90.0, 100.0, 110.0, 120.0})
        for (double sigma : {0.10, 0.25, 0.50})
            for (double T : {0.1, 0.5, 1.0, 2.0})
                for (double r : {0.0, 0.05})
                    for (bool call : {true, false}) grid.push_back({100.0, K, r, sigma, T, call});
    return grid;
}

int main() {
    const auto grid = test_grid();
    double max_greek_err = 0.0;

    // 1) Known textbook value (Hull): S=42, K=40, r=10%, vol=20%, T=0.5 -> call 4.76, put 0.81
    {
        Option c{42, 40, 0.10, 0.20, 0.5, true}, p = c;
        p.is_call = false;
        check(std::fabs(bs_price(c) - 4.76) < 0.005, "hull call", c, bs_price(c), 4.76);
        check(std::fabs(bs_price(p) - 0.81) < 0.005, "hull put", p, bs_price(p), 0.81);
    }

    for (const Option& o : grid) {
        // 2) Put-call parity: C - P = S - K e^{-rT}
        if (o.is_call) {
            Option p = o;
            p.is_call = false;
            const double lhs = bs_price(o) - bs_price(p);
            const double rhs = o.S - o.K * std::exp(-o.r * o.T);
            check(std::fabs(lhs - rhs) < 1e-10, "put-call parity", o, lhs, rhs);
        }

        // 3) Analytic Greeks vs central finite differences of the price.
        const Greeks g = bs_greeks(o);
        auto bump = [&](double Option::*field, double h) {
            Option up = o, dn = o;
            up.*field += h;
            dn.*field -= h;
            return std::make_pair(bs_price(up), bs_price(dn));
        };
        const double hS = 1e-3 * o.S, hv = 1e-4, hT = 1e-5, hr = 1e-5;
        auto [su, sd] = bump(&Option::S, hS);
        auto [vu, vd] = bump(&Option::sigma, hv);
        auto [tu, td] = bump(&Option::T, hT);
        auto [ru, rd] = bump(&Option::r, hr);
        const double fd[5] = {(su - sd) / (2 * hS), (su - 2 * bs_price(o) + sd) / (hS * hS),
                              (vu - vd) / (2 * hv), -(tu - td) / (2 * hT), (ru - rd) / (2 * hr)};
        const double an[5] = {g.delta, g.gamma, g.vega, g.theta, g.rho};
        const char* names[5] = {"delta", "gamma", "vega", "theta", "rho"};
        for (int i = 0; i < 5; ++i) {
            const double err = std::fabs(fd[i] - an[i]);
            max_greek_err = std::max(max_greek_err, err);
            check(err < 1e-4 * std::max(1.0, std::fabs(an[i])), names[i], o, an[i], fd[i]);
        }

        // 4) Monte Carlo agrees with Black-Scholes within 4 standard errors (plus 1e-6 absolute, since
        //    deep out-of-the-money options worth ~1e-8 may see zero paths finish in the money).
        const MCResult mc = mc_price(o, 200000, 42, true);
        const double bs = bs_price(o);
        check(std::fabs(mc.price - bs) < 4 * mc.std_error + 1e-6, "mc vs bs", o, mc.price, bs);
    }

    std::printf("%d/%d checks passed across %zu option configurations (max |analytic - FD| Greek error: %.2e)\n",
                n_checks - n_failed, n_checks, grid.size(), max_greek_err);
    return n_failed == 0 ? 0 : 1;
}
