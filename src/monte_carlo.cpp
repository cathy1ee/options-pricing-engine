#include "pricing.hpp"
#include <algorithm>
#include <cmath>
#include <future>
#include <random>
#include <stdexcept>
#include <vector>

namespace {
// Terminal prices and payoffs under the risk-neutral measure, where
// S_T = S * exp((r - q - sigma^2/2) T + sigma sqrt(T) Z).
struct Model {
    const Option& o;
    double drift, vol, disc;

    explicit Model(const Option& option)
        : o(option),
          drift((o.r - o.q - 0.5 * o.sigma * o.sigma) * o.T),
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
    // Combines two sets of samples, as if every sample had gone through add (Chan et al.).
    void merge(const RunningStats& other) {
        if (other.n == 0) return;
        if (n == 0) {
            *this = other;
            return;
        }
        const double total = (double)n + (double)other.n;
        const double delta = other.mean - mean;
        mean += delta * ((double)other.n / total);
        sq_dev += other.sq_dev + delta * delta * ((double)n * (double)other.n / total);
        n += other.n;
    }

    double std_error() const { return std::sqrt(sq_dev / (n - 1) / n); }
};

// Runs n_samples samples on a random stream of their own.
RunningStats simulate(const Model& model, std::uint64_t n_samples, std::uint64_t seed, bool antithetic) {
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> normal(0.0, 1.0);
    auto payoff_at = [&](double z) { return model.payoff(model.terminal(z)); };

    RunningStats stats;
    for (std::uint64_t i = 0; i < n_samples; ++i) {
        const double z = normal(rng);
        stats.add(antithetic ? 0.5 * (payoff_at(z) + payoff_at(-z)) : payoff_at(z));
    }
    return stats;
}

// Derives a seed for one stream with SplitMix64. Consecutive seeds can start mt19937_64 in related
// states, so each stream's seed is scrambled. Overlap between streams is then extremely unlikely,
// though not impossible, since mt19937_64 has no way to guarantee disjoint streams.
std::uint64_t stream_seed(std::uint64_t seed, std::uint64_t stream) {
    std::uint64_t z = seed + (stream + 1) * 0x9E3779B97F4A7C15ULL;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

// Inverse of the standard normal CDF for p strictly between 0 and 1: Acklam's rational
// approximation, polished to near machine precision with one Halley step.
double norm_inv(double p) {
    constexpr double a[] = {-3.969683028665376e+01, 2.209460984245205e+02,  -2.759285104469687e+02,
                            1.383577518672690e+02,  -3.066479806614716e+01, 2.506628277459239e+00};
    constexpr double b[] = {-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
                            6.680131188771972e+01,  -1.328068155288572e+01};
    constexpr double c[] = {-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
                            -2.549732539343734e+00, 4.374664141464968e+00,  2.938163982698783e+00};
    constexpr double d[] = {7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
                            3.754408661907416e+00};
    constexpr double tail = 0.02425;

    auto tail_value = [&](double q) {
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
               ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
    };

    double x;
    if (p < tail) {
        x = tail_value(std::sqrt(-2.0 * std::log(p)));
    } else if (p > 1.0 - tail) {
        x = -tail_value(std::sqrt(-2.0 * std::log(1.0 - p)));
    } else {
        const double q = p - 0.5, r = q * q;
        x = (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
            (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0);
    }

    const double err = 0.5 * std::erfc(-x / std::sqrt(2.0)) - p;
    const double u = err * std::sqrt(2.0 * 3.14159265358979323846) * std::exp(0.5 * x * x);
    return x - u / (1.0 + 0.5 * x * u);
}

std::uint64_t reverse_bits(std::uint64_t x) {
    x = ((x >> 1) & 0x5555555555555555ULL) | ((x & 0x5555555555555555ULL) << 1);
    x = ((x >> 2) & 0x3333333333333333ULL) | ((x & 0x3333333333333333ULL) << 2);
    x = ((x >> 4) & 0x0F0F0F0F0F0F0F0FULL) | ((x & 0x0F0F0F0F0F0F0F0FULL) << 4);
    x = ((x >> 8) & 0x00FF00FF00FF00FFULL) | ((x & 0x00FF00FF00FF00FFULL) << 8);
    x = ((x >> 16) & 0x0000FFFF0000FFFFULL) | ((x & 0x0000FFFF0000FFFFULL) << 16);
    return (x >> 32) | (x << 32);
}
}  // namespace

double inverse_normal_cdf(double p) {
    if (!(p > 0.0 && p < 1.0)) throw std::invalid_argument("probability must be strictly between 0 and 1");
    return norm_inv(p);
}

MCResult mc_price(const Option& o, std::uint64_t n_paths, std::uint64_t seed, bool antithetic) {
    validate_option(o);

    // With antithetic variates one sample is the average of the payoffs at Z and -Z,
    // so n_paths paths produce n_paths / 2 independent samples.
    const std::uint64_t n_samples = antithetic ? n_paths / 2 : n_paths;
    if (n_samples < 2) throw std::invalid_argument("need at least 2 samples for a standard error");

    const Model model(o);
    const RunningStats stats = simulate(model, n_samples, seed, antithetic);
    return {model.disc * stats.mean, model.disc * stats.std_error()};
}

MCResult mc_price_parallel(const Option& o, std::uint64_t n_paths, std::uint64_t seed, bool antithetic,
                           unsigned n_threads) {
    validate_option(o);
    if (n_threads == 0) throw std::invalid_argument("need at least 1 thread");

    const std::uint64_t n_samples = antithetic ? n_paths / 2 : n_paths;
    if (n_samples < 2) throw std::invalid_argument("need at least 2 samples for a standard error");

    // Samples are split as evenly as possible, and each thread gets its own random stream.
    const Model model(o);
    std::vector<std::future<RunningStats>> workers;
    for (unsigned t = 0; t < n_threads; ++t) {
        const std::uint64_t share = n_samples / n_threads + (t < n_samples % n_threads ? 1 : 0);
        const std::uint64_t worker_seed = stream_seed(seed, t);
        workers.push_back(std::async(std::launch::async, [&model, share, worker_seed, antithetic] {
            return simulate(model, share, worker_seed, antithetic);
        }));
    }

    // Merging in thread order keeps the result the same on every run.
    RunningStats stats;
    for (std::future<RunningStats>& worker : workers) stats.merge(worker.get());
    return {model.disc * stats.mean, model.disc * stats.std_error()};
}

MCResult mc_price_sobol(const Option& o, std::uint64_t n_paths, std::uint64_t seed) {
    validate_option(o);

    // Sobol points are not independent, so their own spread says nothing about the error. Instead
    // each replication scrambles the whole sequence with a random digital shift, which makes the
    // replications' estimates independent and unbiased, and the standard error comes from those.
    constexpr std::uint64_t n_replications = 32;
    const std::uint64_t n_points = n_paths / n_replications;
    if (n_points < 1) throw std::invalid_argument("need at least 32 paths for Sobol replications");

    std::mt19937_64 rng(seed);
    const Model model(o);

    RunningStats stats;
    for (std::uint64_t rep = 0; rep < n_replications; ++rep) {
        const std::uint64_t shift = rng();
        double sum = 0.0;
        for (std::uint64_t i = 0; i < n_points; ++i) {
            // Each path has one normal draw, so only Sobol's first dimension is needed, and that is
            // the base-2 van der Corput sequence: the bits of i reversed. Keeping 52 bits leaves
            // room for the + 0.5 to stay exact, so u can never round to 0 or 1.
            const double u = (((reverse_bits(i) ^ shift) >> 12) + 0.5) * 0x1p-52;
            sum += model.payoff(model.terminal(norm_inv(u)));
        }
        stats.add(sum / n_points);
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
    const double forward = o.S * std::exp((o.r - o.q) * o.T);
    RunningStats stats;
    for (std::uint64_t i = n_pilot; i < n_paths; ++i) {
        const double st = model.terminal(normal(rng));
        stats.add(model.payoff(st) - beta * (st - forward));
    }
    return {model.disc * stats.mean, model.disc * stats.std_error()};
}
