// Measures Monte Carlo convergence and speed against the Black-Scholes benchmark.
// Build and run with `make bench`.
#include "pricing.hpp"
#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>

int main() {
    const Option o{100.0, 100.0, 0.05, 0.20, 1.0, true};  // at-the-money 1y call
    const double bs = bs_price(o);
    const int n_seeds = 20;  // repeat each run with different seeds to measure typical error
    std::printf("Black-Scholes price: %.6f\n\n", bs);
    std::printf("%10s | %-10s | %12s | %12s | %10s\n", "paths", "method", "RMSE vs BS", "RMSE % of px",
                "ms / run");
    std::printf("-----------+------------+--------------+--------------+-----------\n");

    std::vector<double> log_n, log_err;
    for (std::uint64_t n : {1000ULL, 10000ULL, 100000ULL, 1000000ULL, 10000000ULL}) {
        for (bool anti : {false, true}) {
            double sq_err = 0.0;
            auto t0 = std::chrono::steady_clock::now();
            for (int s = 0; s < n_seeds; ++s) {
                const double e = mc_price(o, n, 1000 + s, anti).price - bs;
                sq_err += e * e;
            }
            auto t1 = std::chrono::steady_clock::now();
            const double rmse = std::sqrt(sq_err / n_seeds);
            const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count() / n_seeds;
            std::printf("%10llu | %-10s | %12.6f | %11.4f%% | %10.2f\n", (unsigned long long)n,
                        anti ? "antithetic" : "plain", rmse, 100 * rmse / bs, ms);
            if (!anti) {
                log_n.push_back(std::log10((double)n));
                log_err.push_back(std::log10(rmse));
            }
        }
    }

    // Least-squares slope of log(error) vs log(paths); theory says -0.5 (error ~ 1/sqrt(N)).
    double mx = 0, my = 0;
    for (size_t i = 0; i < log_n.size(); ++i) mx += log_n[i], my += log_err[i];
    mx /= log_n.size(), my /= log_n.size();
    double num = 0, den = 0;
    for (size_t i = 0; i < log_n.size(); ++i) num += (log_n[i] - mx) * (log_err[i] - my), den += (log_n[i] - mx) * (log_n[i] - mx);
    std::printf("\nConvergence slope (log RMSE vs log N): %.3f  (theory: -0.500)\n", num / den);

    // Variance reduction: same standard error needs (se_plain/se_anti)^2 fewer samples.
    const MCResult p = mc_price(o, 1000000, 7, false), a = mc_price(o, 1000000, 7, true);
    std::printf("Std error at 1M paths: plain %.5f, antithetic %.5f -> variance reduced %.1f%%\n", p.std_error,
                a.std_error, 100 * (1 - (a.std_error * a.std_error) / (p.std_error * p.std_error)));

    // Closed-form throughput.
    volatile double sink = 0;
    const int n_bs = 5000000;
    auto t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < n_bs; ++i) {
        Option x = o;
        x.K = 80.0 + (i % 400) * 0.1;
        sink = sink + bs_price(x);
    }
    auto t1 = std::chrono::steady_clock::now();
    const double sec = std::chrono::duration<double>(t1 - t0).count();
    std::printf("Closed-form: %.1fM prices/sec (%.0f ns each)\n", n_bs / sec / 1e6, 1e9 * sec / n_bs);
    return 0;
}
