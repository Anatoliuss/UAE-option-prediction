#pragma once
#include <random>
#include <cmath>
#include <cstdint>
#include <utility>
#include <limits>
#include <algorithm>

inline std::pair<double,double> mc_euro_option_price(
    double S, double K, double r, double q, double sigma, double T,
    bool is_call,
    std::uint64_t n_paths = 100'000,
    std::uint64_t seed    = 42,
    bool antithetic       = true)
{
    if (S <= 0.0 || K <= 0.0)
        return {std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN()};

    if (T <= 0.0) {
        double payoff = is_call ? std::max(0.0, S - K)
                                : std::max(0.0, K - S);
        return {payoff, 0.0};
    }

    if (sigma <= 0.0) {
        const double ST = S * std::exp((r - q) * T);
        const double disc = std::exp(-r * T);
        double payoff = is_call ? std::max(0.0, ST - K)
                                : std::max(0.0, K - ST);
        return {disc * payoff, 0.0};
    }

    if (n_paths == 0)
        return {std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN()};

    std::mt19937_64 gen(seed);
    std::normal_distribution<double> norm01(0.0, 1.0);

    const double drift = (r - q - 0.5 * sigma * sigma) * T;
    const double volT  = sigma * std::sqrt(T);
    const double disc  = std::exp(-r * T);

    auto payoff_disc = [&](double Z){
        const double ST = S * std::exp(drift + volT * Z);
        const double payoff = is_call ? std::max(0.0, ST - K)
                                      : std::max(0.0, K - ST);
        return disc * payoff;
    };

    double mean = 0.0;
    double M2   = 0.0;    
    std::uint64_t n = 0;

    const bool use_pairs = antithetic && (n_paths >= 2);
    const std::uint64_t loops = use_pairs ? n_paths / 2 : n_paths;

    for (std::uint64_t i = 0; i < loops; ++i) {
        const double z  = norm01(gen);
        if (use_pairs) {
            const double y1 = payoff_disc(z);
            const double y2 = payoff_disc(-z);
            const double y  = 0.5 * (y1 + y2);

            ++n;
            const double delta  = y - mean;
            mean += delta / static_cast<double>(n);
            const double delta2 = y - mean;
            M2   += delta * delta2;
        } else {
            const double y = payoff_disc(z);
            ++n;
            const double delta  = y - mean;
            mean += delta / static_cast<double>(n);
            const double delta2 = y - mean;
            M2   += delta * delta2;
        }
    }

    const double var = (n > 1) ? (M2 / static_cast<double>(n - 1)) : 0.0;
    const double se  = (n > 0) ? std::sqrt(var / static_cast<double>(n)) : 0.0;

    return {mean, se};
}
