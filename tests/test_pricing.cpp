// Self-contained tests: no framework needed. Build with `make test`.
#include "pricing.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

static int n_checks = 0, n_failed = 0;

template <class F>
static bool throws_invalid_argument(F f) {
    try {
        f();
    } catch (const std::invalid_argument&) {
        return true;
    }
    return false;
}

static void check(bool ok, const char* what, const Option& o, double got, double want) {
    ++n_checks;
    if (!ok) {
        ++n_failed;
        std::printf("FAIL %-18s S=%.0f K=%.0f r=%.2f vol=%.2f T=%.2f q=%.2f %s: got %.8f want %.8f\n", what, o.S,
                    o.K, o.r, o.sigma, o.T, o.q, o.is_call ? "call" : "put", got, want);
    }
}

// Grid of market inputs covering in/at/out of the money, low/high vol, short/long expiry, with and
// without dividends.
static std::vector<Option> test_grid() {
    std::vector<Option> grid;
    for (double K : {80.0, 90.0, 100.0, 110.0, 120.0})
        for (double sigma : {0.10, 0.25, 0.50})
            for (double T : {0.1, 0.5, 1.0, 2.0})
                for (double r : {0.0, 0.05})
                    for (bool call : {true, false})
                        for (double q : {0.0, 0.03}) grid.push_back({100.0, K, r, sigma, T, call, q});
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
        // 2) Put-call parity: C - P = S e^{-qT} - K e^{-rT}
        if (o.is_call) {
            Option p = o;
            p.is_call = false;
            const double lhs = bs_price(o) - bs_price(p);
            const double rhs = o.S * std::exp(-o.q * o.T) - o.K * std::exp(-o.r * o.T);
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
            check(err < 5e-6 + 1e-4 * std::fabs(an[i]), names[i], o, an[i], fd[i]);
        }

        // 4) Monte Carlo agrees with Black-Scholes within 4 standard errors (plus 1e-6 absolute, since
        //    deep out-of-the-money options worth ~1e-8 may see zero paths finish in the money).
        const MCResult mc = mc_price(o, 200000, 42, true);
        const double bs = bs_price(o);
        check(std::fabs(mc.price - bs) < 4 * mc.std_error + 1e-6, "mc vs bs", o, mc.price, bs);
        const MCResult plain = mc_price(o, 200000, 42, false);
        check(std::fabs(plain.price - bs) < 4 * plain.std_error + 1e-6, "plain mc vs bs", o, plain.price, bs);
    }

    // 5) Limits with no randomness left: T = 0 or sigma = 0.
    {
        struct Limit {
            Option o;
            double price, delta;
        };
        const double pv_k = 100.0 * std::exp(-0.05);
        const Limit limits[] = {
            {{100, 100, 0.05, 0.20, 0, true}, 0.0, 0.5},
            {{100, 100, 0.05, 0.20, 0, false}, 0.0, -0.5},
            {{110, 100, 0.05, 0.20, 0, true}, 10.0, 1.0},
            {{110, 100, 0.05, 0.20, 0, false}, 0.0, 0.0},
            {{90, 100, 0.05, 0.20, 0, true}, 0.0, 0.0},
            {{90, 100, 0.05, 0.20, 0, false}, 10.0, -1.0},
            {{110, 100, 0.05, 0, 1, true}, 110.0 - pv_k, 1.0},
            {{110, 100, 0.05, 0, 1, false}, 0.0, 0.0},
            {{90, 100, 0.05, 0, 1, true}, 0.0, 0.0},
            {{90, 100, 0.05, 0, 1, false}, pv_k - 90.0, -1.0},
            {{100, 100, 0.0, 0, 1, true}, 0.0, 0.5},
            {{100, 100, 0.0, 0, 1, false}, 0.0, -0.5},
        };
        for (const Limit& l : limits) {
            const double price = bs_price(l.o);
            const Greeks g = bs_greeks(l.o);
            check(std::fabs(price - l.price) < 1e-12, "limit price", l.o, price, l.price);
            check(std::fabs(g.delta - l.delta) < 1e-12, "limit delta", l.o, g.delta, l.delta);
            check(g.gamma == 0.0, "limit gamma", l.o, g.gamma, 0.0);
            check(std::isfinite(g.vega) && std::isfinite(g.theta) && std::isfinite(g.rho), "limit finite", l.o,
                  g.theta, 0.0);
        }
    }

    // 6) Invalid inputs are rejected by every pricer.
    {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        const Option bad[] = {{0, 100, 0.05, 0.20, 1, true},    {-1, 100, 0.05, 0.20, 1, true},
                              {100, 0, 0.05, 0.20, 1, true},    {100, 100, 0.05, 0.20, -1, true},
                              {100, 100, 0.05, -0.20, 1, true}, {100, 100, 0.05, nan, 1, true},
                              {100, 100, 0.05, 0.20, 1, true, nan}};
        for (const Option& o : bad) {
            check(throws_invalid_argument([&] { bs_price(o); }), "bad bs_price", o, 0, 0);
            check(throws_invalid_argument([&] { bs_greeks(o); }), "bad bs_greeks", o, 0, 0);
            check(throws_invalid_argument([&] { mc_price(o, 1000, 42, false); }), "bad mc_price", o, 0, 0);
        }
    }

    // 7) Monte Carlo path counts: too few samples throw, the minimum works, odd antithetic rounds down.
    {
        const Option o{100, 100, 0.05, 0.20, 1, true};
        const struct {
            std::uint64_t n_paths;
            bool antithetic;
        } too_few[] = {{0, false}, {1, false}, {0, true}, {1, true}, {2, true}, {3, true}};
        for (const auto& c : too_few)
            check(throws_invalid_argument([&] { mc_price(o, c.n_paths, 42, c.antithetic); }), "too few paths", o,
                  (double)c.n_paths, 0);

        const MCResult plain = mc_price(o, 2, 42, false), anti = mc_price(o, 4, 42, true);
        check(std::isfinite(plain.price) && std::isfinite(plain.std_error), "min paths plain", o, plain.std_error, 0);
        check(std::isfinite(anti.price) && std::isfinite(anti.std_error), "min paths anti", o, anti.std_error, 0);

        const double odd = mc_price(o, 1001, 42, true).price, even = mc_price(o, 1000, 42, true).price;
        check(odd == even, "odd antithetic", o, odd, even);
    }

    // 8) Standard error stays accurate when the payoffs barely vary (a running sum of squares does not).
    {
        for (const Option& o : {Option{110, 100, 0.05, 0, 1, true}, Option{100, 50, 0.05, 0, 2, true}})
            for (bool anti : {false, true}) {
                const double se = mc_price(o, 1000000, 42, anti).std_error;
                check(se == 0.0, "zero vol std error", o, se, 0.0);
            }

        // Deep in the money the discounted payoff is S (1 + sigma Z) - K e^{-rT}, so se = S sigma / sqrt(N).
        const Option o{1000, 50, 0.05, 1e-8, 1, true};
        const double se = mc_price(o, 1000000, 42, false).std_error, want = 1000 * 1e-8 / 1000;
        check(std::fabs(se - want) < 0.05 * want, "low vol std error", o, se, want);
    }

    // 9) The reported standard error means what it says.
    for (const Option& o : {Option{100, 100, 0.05, 0.20, 1, true}, Option{100, 90, 0.05, 0.25, 0.5, false},
                            Option{100, 80, 0.0, 0.50, 2, true}}) {
        // It matches the observed spread of prices across seeds.
        const int n_seeds = 100;
        double sum = 0.0, sum_sq_dev = 0.0, sum_se = 0.0;
        std::vector<double> prices;
        for (int s = 0; s < n_seeds; ++s) {
            const MCResult m = mc_price(o, 10000, 1000 + s, false);
            prices.push_back(m.price);
            sum += m.price;
            sum_se += m.std_error;
        }
        for (double p : prices) sum_sq_dev += (p - sum / n_seeds) * (p - sum / n_seeds);
        const double spread = std::sqrt(sum_sq_dev / (n_seeds - 1)), se = sum_se / n_seeds;
        check(std::fabs(spread / se - 1.0) < 0.25, "se vs spread", o, se, spread);

        // It shrinks like 1/sqrt(N).
        const MCResult small = mc_price(o, 10000, 7, false), large = mc_price(o, 100000, 7, false);
        const double shrink = small.std_error / large.std_error;
        check(std::fabs(shrink - std::sqrt(10.0)) < 0.1 * std::sqrt(10.0), "se scaling", o, shrink, std::sqrt(10.0));

        // Antithetic variates reduce it for these monotonic payoffs.
        const MCResult anti = mc_price(o, 100000, 7, true);
        check(anti.std_error < large.std_error, "antithetic se", o, anti.std_error, large.std_error);

        // The same seed gives the same result.
        check(mc_price(o, 1000, 3, false).price == mc_price(o, 1000, 3, false).price, "plain seed", o, 0, 0);
        check(mc_price(o, 1000, 3, true).price == mc_price(o, 1000, 3, true).price, "antithetic seed", o, 0, 0);
    }

    // 10) Implied volatility round trips. Where vega is tiny the price no longer pins down the
    //     volatility, so matching the price is the most that can be asked for.
    for (const Option& o : grid) {
        const double price = bs_price(o);
        Option back = o;
        back.sigma = implied_vol(o, price);
        check(std::fabs(back.sigma - o.sigma) < 1e-8 || std::fabs(bs_price(back) - price) < 1e-12, "implied vol", o,
              back.sigma, o.sigma);
    }

    // 11) Implied volatility away from the grid, and at and beyond its bounds.
    {
        for (const Option& o : {Option{100, 100, 0.05, 3.0, 1, true}, Option{100, 100, 0.05, 0.01, 0.01, true},
                                Option{100, 150, 0.05, 0.20, 0.5, true}}) {
            const double iv = implied_vol(o, bs_price(o));
            check(std::fabs(iv - o.sigma) < 1e-8, "implied vol extreme", o, iv, o.sigma);
        }

        const Option c{100, 100, 0.05, 0.20, 1, true}, p{100, 100, 0.05, 0.20, 1, false};
        Option no_vol = c;
        no_vol.sigma = std::numeric_limits<double>::quiet_NaN();
        const double iv = implied_vol(no_vol, bs_price(c));
        check(std::fabs(iv - c.sigma) < 1e-8, "implied vol ignores sigma", c, iv, c.sigma);

        Option flat = c;
        flat.sigma = 0.0;
        const double at_lower = implied_vol(c, bs_price(flat));
        check(at_lower == 0.0, "implied vol at lower bound", c, at_lower, 0.0);

        const Option expired{100, 100, 0.05, 0.20, 0, true};
        check(throws_invalid_argument([&] { implied_vol(c, 100.0); }), "iv at upper bound", c, 100.0, 0);
        check(throws_invalid_argument([&] { implied_vol(c, 1.0); }), "iv below lower bound", c, 1.0, 0);
        check(throws_invalid_argument([&] { implied_vol(p, 96.0); }), "iv above put bound", p, 96.0, 0);
        check(throws_invalid_argument([&] { implied_vol(p, -1.0); }), "iv negative price", p, -1.0, 0);
        check(throws_invalid_argument([&] { implied_vol(c, no_vol.sigma); }), "iv nan price", c, 0, 0);
        check(throws_invalid_argument([&] { implied_vol(expired, 1.0); }), "iv at expiry", expired, 1.0, 0);
    }

    // 12) Binomial tree: European prices approach Black-Scholes, and American prices respect their bounds.
    for (const Option& o : grid) {
        const double european = binomial_price(o, 500, Exercise::European), bs = bs_price(o);
        check(std::fabs(european - bs) < 0.02, "tree vs bs", o, european, bs);

        const double eu = binomial_price(o, 200, Exercise::European), am = binomial_price(o, 200, Exercise::American);
        if (o.is_call && o.q == 0.0) {
            // Without dividends early exercise never pays for a call.
            check(std::fabs(am - eu) < 1e-9, "american call", o, am, eu);
        } else {
            const double intrinsic = o.is_call ? o.S - o.K : o.K - o.S;
            check(am >= eu - 1e-12, "american >= european", o, am, eu);
            check(am >= intrinsic - 1e-12, "american >= intrinsic", o, am, intrinsic);
        }
    }

    // 13) Binomial tree: convergence, Hull's American put, limits and invalid inputs.
    {
        const Option atm{100, 100, 0.05, 0.20, 1, true};
        const double bs = bs_price(atm);
        const double coarse = std::fabs(binomial_price(atm, 100, Exercise::European) - bs);
        const double fine = std::fabs(binomial_price(atm, 1600, Exercise::European) - bs);
        check(fine < coarse / 10, "tree convergence", atm, fine, coarse);

        // Hull: S=K=50, r=10%, vol=40%, T=5 months.
        const Option hull{50, 50, 0.10, 0.40, 0.4167, false};
        const struct {
            std::uint32_t steps;
            double price;
        } hull_values[] = {{5, 4.49}, {30, 4.263}, {50, 4.272}, {100, 4.278}, {500, 4.283}};
        for (const auto& h : hull_values) {
            const double am = binomial_price(hull, h.steps, Exercise::American);
            check(std::fabs(am - h.price) < 0.002, "hull american put", hull, am, h.price);
        }

        const double pv_k = 100.0 * std::exp(-0.05);
        const Option flat_put{90, 100, 0.05, 0, 1, false}, flat_call{110, 100, 0.05, 0, 1, true};
        const Option expired_put{90, 100, 0.05, 0.20, 0, false};
        const double flat_am = binomial_price(flat_put, 10, Exercise::American);
        const double flat_eu = binomial_price(flat_put, 10, Exercise::European);
        const double flat_call_am = binomial_price(flat_call, 10, Exercise::American);
        const double expired_am = binomial_price(expired_put, 10, Exercise::American);
        check(std::fabs(flat_am - 10.0) < 1e-12, "american put at zero vol", flat_put, flat_am, 10.0);
        check(std::fabs(flat_eu - (pv_k - 90.0)) < 1e-12, "european put at zero vol", flat_put, flat_eu, pv_k - 90.0);
        check(std::fabs(flat_call_am - (110.0 - pv_k)) < 1e-12, "american call at zero vol", flat_call, flat_call_am,
              110.0 - pv_k);
        check(std::fabs(expired_am - 10.0) < 1e-12, "american put at expiry", expired_put, expired_am, 10.0);

        const Option bad{0, 100, 0.05, 0.20, 1, true}, low_vol{100, 100, 0.05, 0.01, 10, true};
        check(throws_invalid_argument([&] { binomial_price(atm, 0, Exercise::European); }), "tree no steps", atm, 0, 0);
        check(throws_invalid_argument([&] { binomial_price(bad, 10, Exercise::American); }), "bad tree", bad, 0, 0);
        check(throws_invalid_argument([&] { binomial_price(low_vol, 1, Exercise::European); }), "tree too coarse",
              low_vol, 0, 0);
    }

    // 14) Control variate agrees with Black-Scholes within 4 standard errors, as in check 4.
    for (const Option& o : grid) {
        const MCResult cv = mc_price_control_variate(o, 200000, 42);
        const double bs = bs_price(o);
        check(std::fabs(cv.price - bs) < 4 * cv.std_error + 1e-6, "control variate vs bs", o, cv.price, bs);
    }

    // 15) Control variate standard error: honest, smaller than plain, and reproducible.
    for (const Option& o : {Option{100, 100, 0.05, 0.20, 1, true}, Option{100, 90, 0.05, 0.25, 0.5, false},
                            Option{100, 80, 0.0, 0.50, 2, true}}) {
        const int n_seeds = 100;
        double sum = 0.0, sum_sq_dev = 0.0, sum_se = 0.0;
        std::vector<double> prices;
        for (int s = 0; s < n_seeds; ++s) {
            const MCResult m = mc_price_control_variate(o, 10000, 1000 + s);
            prices.push_back(m.price);
            sum += m.price;
            sum_se += m.std_error;
        }
        for (double p : prices) sum_sq_dev += (p - sum / n_seeds) * (p - sum / n_seeds);
        const double spread = std::sqrt(sum_sq_dev / (n_seeds - 1)), se = sum_se / n_seeds;
        check(std::fabs(spread / se - 1.0) < 0.25, "cv se vs spread", o, se, spread);

        const MCResult cv = mc_price_control_variate(o, 100000, 7), plain = mc_price(o, 100000, 7, false);
        check(cv.std_error < plain.std_error, "cv se vs plain", o, cv.std_error, plain.std_error);
        check(mc_price_control_variate(o, 1000, 3).price == mc_price_control_variate(o, 1000, 3).price, "cv seed", o, 0,
              0);
    }

    // 16) Control variate at zero volatility, at its minimum path count, and with invalid inputs.
    {
        const Option flat{110, 100, 0.05, 0, 1, true}, o{100, 100, 0.05, 0.20, 1, true}, bad{0, 100, 0.05, 0.20, 1, true};
        const double flat_se = mc_price_control_variate(flat, 10000, 42).std_error;
        check(flat_se == 0.0, "cv zero vol", flat, flat_se, 0.0);

        for (std::uint64_t n : {0ULL, 1ULL, 2ULL, 3ULL})
            check(throws_invalid_argument([&] { mc_price_control_variate(o, n, 42); }), "cv too few paths", o, (double)n, 0);

        const MCResult fewest = mc_price_control_variate(o, 4, 42);
        check(std::isfinite(fewest.price) && std::isfinite(fewest.std_error), "cv min paths", o, fewest.std_error, 0);
        check(throws_invalid_argument([&] { mc_price_control_variate(bad, 1000, 42); }), "bad cv", bad, 0, 0);
    }

    // 17) Sobol agrees with Black-Scholes. The bound is 5 standard errors rather than 4 because the
    //     standard error comes from only 32 replications, which gives it heavier tails.
    for (const Option& o : grid) {
        const MCResult sobol = mc_price_sobol(o, 1 << 16, 42);
        const double bs = bs_price(o);
        check(std::fabs(sobol.price - bs) < 5 * sobol.std_error + 1e-6, "sobol vs bs", o, sobol.price, bs);
    }

    // 18) Sobol standard error: honest, smaller than plain, faster than 1/sqrt(N), and reproducible.
    for (const Option& o : {Option{100, 100, 0.05, 0.20, 1, true}, Option{100, 90, 0.05, 0.25, 0.5, false},
                            Option{100, 80, 0.0, 0.50, 2, true}}) {
        const int n_seeds = 100;
        double sum = 0.0, sum_sq_dev = 0.0, sum_se = 0.0;
        std::vector<double> prices;
        for (int s = 0; s < n_seeds; ++s) {
            const MCResult m = mc_price_sobol(o, 1 << 12, 1000 + s);
            prices.push_back(m.price);
            sum += m.price;
            sum_se += m.std_error;
        }
        for (double p : prices) sum_sq_dev += (p - sum / n_seeds) * (p - sum / n_seeds);
        const double spread = std::sqrt(sum_sq_dev / (n_seeds - 1)), se = sum_se / n_seeds;
        check(std::fabs(spread / se - 1.0) < 0.25, "sobol se vs spread", o, se, spread);

        const MCResult sobol = mc_price_sobol(o, 1 << 16, 7), plain = mc_price(o, 1 << 16, 7, false);
        check(sobol.std_error < plain.std_error, "sobol se vs plain", o, sobol.std_error, plain.std_error);

        // 16x the paths would shrink a random estimate's error 4x.
        const double coarse = mc_price_sobol(o, 1 << 12, 7).std_error;
        check(sobol.std_error < coarse / 5, "sobol convergence", o, sobol.std_error, coarse);

        check(mc_price_sobol(o, 1024, 3).price == mc_price_sobol(o, 1024, 3).price, "sobol seed", o, 0, 0);
    }

    // 19) Sobol at zero volatility, at its minimum path count, and with invalid inputs.
    {
        const Option flat{110, 100, 0.05, 0, 1, true}, o{100, 100, 0.05, 0.20, 1, true}, bad{0, 100, 0.05, 0.20, 1, true};
        const double flat_se = mc_price_sobol(flat, 1024, 42).std_error;
        check(flat_se == 0.0, "sobol zero vol", flat, flat_se, 0.0);

        check(throws_invalid_argument([&] { mc_price_sobol(o, 0, 42); }), "sobol no paths", o, 0, 0);
        check(throws_invalid_argument([&] { mc_price_sobol(o, 31, 42); }), "sobol too few paths", o, 31, 0);

        const MCResult fewest = mc_price_sobol(o, 32, 42);
        check(std::isfinite(fewest.price) && std::isfinite(fewest.std_error), "sobol min paths", o, fewest.std_error, 0);
        check(throws_invalid_argument([&] { mc_price_sobol(bad, 1024, 42); }), "bad sobol", bad, 0, 0);
    }

    // 20) Parallel Monte Carlo agrees with Black-Scholes within 4 standard errors, as in check 4.
    for (const Option& o : grid) {
        const MCResult par = mc_price_parallel(o, 200000, 42, true, 4);
        const double bs = bs_price(o);
        check(std::fabs(par.price - bs) < 4 * par.std_error + 1e-6, "parallel vs bs", o, par.price, bs);
    }

    // 21) Parallel results are reproducible, agree with one thread, and merge to the right standard error.
    for (const Option& o : {Option{100, 100, 0.05, 0.20, 1, true}, Option{100, 90, 0.05, 0.25, 0.5, false},
                            Option{100, 80, 0.0, 0.50, 2, true}}) {
        for (bool anti : {false, true}) {
            const double first = mc_price_parallel(o, 100000, 3, anti, 4).price;
            const double second = mc_price_parallel(o, 100000, 3, anti, 4).price;
            check(first == second, "parallel seed", o, first, second);
        }

        const MCResult single = mc_price(o, 200000, 7, false), par = mc_price_parallel(o, 200000, 7, false, 3);
        const double se = std::sqrt(single.std_error * single.std_error + par.std_error * par.std_error);
        check(std::fabs(single.price - par.price) < 4 * se, "parallel vs single", o, par.price, single.price);
        check(std::fabs(par.std_error / single.std_error - 1.0) < 0.05, "parallel se", o, par.std_error,
              single.std_error);
    }

    // 22) Parallel with more threads than samples, at zero volatility, and with invalid inputs.
    {
        const Option flat{110, 100, 0.05, 0, 1, true}, o{100, 100, 0.05, 0.20, 1, true}, bad{0, 100, 0.05, 0.20, 1, true};
        const MCResult crowded = mc_price_parallel(o, 10, 42, false, 16);
        check(std::isfinite(crowded.price) && std::isfinite(crowded.std_error), "parallel idle threads", o,
              crowded.std_error, 0);

        const double flat_se = mc_price_parallel(flat, 10000, 42, false, 4).std_error;
        check(flat_se == 0.0, "parallel zero vol", flat, flat_se, 0.0);

        check(throws_invalid_argument([&] { mc_price_parallel(o, 1000, 42, false, 0); }), "parallel no threads", o, 0, 0);
        check(throws_invalid_argument([&] { mc_price_parallel(o, 1, 42, false, 4); }), "parallel too few", o, 1, 0);
        check(throws_invalid_argument([&] { mc_price_parallel(bad, 1000, 42, false, 4); }), "bad parallel", bad, 0, 0);
    }

    // 23) A dividend yield is equivalent to a lower spot with no dividends.
    for (const Option& o : grid) {
        if (o.q == 0.0) continue;
        Option no_div = o;
        no_div.S = o.S * std::exp(-o.q * o.T);
        no_div.q = 0.0;
        check(std::fabs(bs_price(o) - bs_price(no_div)) < 1e-10, "dividend as spot", o, bs_price(o), bs_price(no_div));
    }

    // 24) Dividends in the limits, and the early exercise they make worthwhile.
    {
        const Option flat{110, 100, 0.05, 0, 1, true, 0.02};
        const double want = 110.0 * std::exp(-0.02) - 100.0 * std::exp(-0.05);
        check(std::fabs(bs_price(flat) - want) < 1e-12, "dividend limit price", flat, bs_price(flat), want);
        const double delta = bs_greeks(flat).delta;
        check(std::fabs(delta - std::exp(-0.02)) < 1e-12, "dividend limit delta", flat, delta, std::exp(-0.02));

        const Option deep{100, 80, 0.05, 0.20, 1, true, 0.10};
        const double am = binomial_price(deep, 200, Exercise::American);
        const double eu = binomial_price(deep, 200, Exercise::European);
        check(am > eu + 0.01, "american call with dividends", deep, am, eu);

        // At zero volatility the best time to exercise can fall strictly between now and expiry.
        const Option timed{100, 100, 0.05, 0, 40, true, 0.02};
        double best = 0.0;
        for (int i = 0; i <= 100000; ++i) {
            const double t = 40.0 * i / 100000;
            best = std::max(best, 100.0 * std::exp(-0.02 * t) - 100.0 * std::exp(-0.05 * t));
        }
        const double timed_am = binomial_price(timed, 10, Exercise::American);
        check(std::fabs(timed_am - best) < 1e-6, "american exercise time", timed, timed_am, best);
    }

    std::printf("%d/%d checks passed across %zu option configurations (max |analytic - FD| Greek error: %.2e)\n",
                n_checks - n_failed, n_checks, grid.size(), max_greek_err);
    return n_failed == 0 ? 0 : 1;
}
