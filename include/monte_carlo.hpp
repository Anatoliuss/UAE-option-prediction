#pragma once
#include <random>
#include <cmath>
#include <cstdint>
#include <utility>
#include <limits>
#include <algorithm>
// Monte Carlo option pricing
// Assumes stock price follows Geometric Brownian Motion (GBM) (same as Black-Scholes)
// Idea: simulate many paths of the stock price at expiration, then average the payoffs
//Difference from Black-Scholes: we simulate paths instead of using a closed-form formula
//currently the best method for pricing options with complex features (like American options, barriers, etc.)
// Parameters:
// S = current stock price  
// K = strike price
// r = risk-free interest rate
// q = continuous dividend yield
// sigma = volatility of the underlying asset
// T = time to expiration

inline std::pair<double,double> mc_euro_option_price(
    double S, double K, double r, double q, double sigma, double T,
    bool is_call,
    std::uint64_t n_paths = 100'000, // number of simulated paths
    std::uint64_t seed = 42, // random seed for reproducibility
    bool antithetic = true) // antithetic - use pairs of paths to reduce variance
{
    if (S <= 0.0 || K <= 0.0) // invalid stock or strike price
        return {std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN()};

    if (T <= 0.0) // if time to expiration is zero, we can use the intrinsic value
    {
        double payoff;
        if (is_call)
        {
            payoff = std::max(0.0, S - K); // same as binomial or Black-Scholes
        }
        else
        {
            payoff = std::max(0.0, K - S);
        }
        return std::make_pair(payoff, 0.0);
    }

    if (sigma <= 0.0) { // if volatility is zero, we can use the deterministic forward price
        // ST = S * exp((r - q) * T)
        const double ST = S * std::exp((r - q) * T);
        const double disc = std::exp(-r * T);
        double payoff;
        if (is_call) {
            payoff = std::max(0.0, ST - K);
        } else {
            payoff = std::max(0.0, K - ST);
        }
        return {disc * payoff, 0.0};
    }

    if (n_paths == 0) // no paths to simulate
        return {std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN()};

    std::mt19937_64 gen(seed); //random number generator
    // Normal distribution for standard normal random variables
    std::normal_distribution<double> norm01(0.0, 1.0);

    const double drift = (r - q - 0.5 * sigma * sigma) * T; // drift term (explanation: drift is the expected return of the stock price)
    // volT = sigma * sqrt(T) is the volatility term for the GBM model
    const double volT  = sigma * std::sqrt(T);
    const double disc  = std::exp(-r * T);

    auto payoff_disc = [&](double Z){ // Z is a standard normal random variable
        // ST = S * exp((r - q) * T + sigma * sqrt(T
        double ST = S * std::exp(drift + volT * Z);
        double payoff;
        if (is_call) {
            payoff = std::max(0.0, ST - K);
        } else {
            payoff = std::max(0.0, K - ST);
        }
        return disc * payoff;
    };

    double mean = 0.0;
    double M2   = 0.0;    
    std::uint64_t n = 0;

    const bool use_pairs = antithetic && (n_paths >= 2); // use pairs of paths to reduce variance
    std::uint64_t loops; // number of iterations to run
    if (use_pairs) {
        loops = n_paths / 2;
    } else {
        loops = n_paths;
    }

    for (std::uint64_t i = 0; i < loops; ++i) { 
        const double z  = norm01(gen);
        if (use_pairs) {
            const double y1 = payoff_disc(z); // first path
            const double y2 = payoff_disc(-z);// antithetic path
            // Average the two paths to reduce variance
            const double y  = 0.5 * (y1 + y2);

            n++;
            const double delta  = y - mean;
            mean += delta / static_cast<double>(n); // update mean
            const double delta2 = y - mean;
            M2   += delta * delta2; //M2 is the second moment (]the expected value of the squared deviation) about the mean
        } else { // single path simulation
            const double y = payoff_disc(z);
            n++;
            const double delta  = y - mean; 
            mean += delta / static_cast<double>(n);
            const double delta2 = y - mean;
            M2   += delta * delta2;
        }
    }

    double var = 0.0;
    double se = 0.0;
    if (n > 1) {
        var = M2 / (n - 1);
        se = std::sqrt(var / n);
    }

    return {mean, se}; // return mean and standard error
}
