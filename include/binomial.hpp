#pragma once
#include <cmath>
#include <vector>
#include <limits>
#include <algorithm>
//Binomial option pricing model
//Assumption: stock price follows a binomial distribution (up or down at each discrete step)
//over time, we get a binomial tree of possible stock prices at expiration
//it is basically a discrete replication of black-scholes (which is continuous time)
//Parameters:
// S = current stock price
// K = strike price
// r = risk-free interest rate
// σ = volatility of the underlying asset
// T = time to expiration
// N = number of steps in the binomial tree
inline double binomial_price(
    double S, double K, double r, double sigma, double T,
    int steps, bool is_call, bool american = false, double q = 0.0)
{
    // Check for invalid input values
    if (steps <= 0 || S <= 0.0 || K <= 0.0 || T <= 0.0 || sigma <= 0.0) {
        // If any parameter is not valid, return NaN
        return std::numeric_limits<double>::quiet_NaN();
    }
    //We want the binomial variance to match the GBM variance in delta t
    const double dt   = T / steps; // Time step
    const double u    = std::exp(sigma * std::sqrt(dt)); // Up factor
    const double d    = 1.0 / u; // Down factor
    const double disc = std::exp(-r * dt); // Discount factor per step
    const double p    = (std::exp((r - q) * dt) - d) / (u - d); // Risk-neutral probability of up move
    if (!(p > 0.0 && p < 1.0))
        return std::numeric_limits<double>::quiet_NaN(); // Probability must be in (0, 1)

    // V[j] means node with j UP moves (more standard)
    std::vector<double> V(steps + 1);
    for (int j = 0; j <= steps; ++j) { // fill the vector with option payoffs at maturity (leaf nodes)
        const double ST = S * std::pow(u, j) * std::pow(d, steps - j);
        if (is_call) {
            V[j] = std::max(ST - K, 0.0); // Call option payoff.
        } else {
            V[j] = std::max(K - ST, 0.0); // Put option payoff
        }
    }
    // Backward induction to calculate option price at time 0
    // Idea: we go from the leaf nodes (maturity) to the root node (current time)
    for (int i = steps - 1; i >= 0; --i) { //iterate backwards from maturity to time 0
        for (int j = 0; j <= i; ++j) { // j is the number of up moves at this node
            // Calculate the expected value at this node
            const double cont = disc * (p * V[j + 1] + (1.0 - p) * V[j]);
            if (american) { // If American option , we can exercise early
                double Sij = S * std::pow(u, j) * std::pow(d, i - j);
                double ex; // exercise value
                if (is_call) {
                    ex = Sij - K;
                    if (ex < 0.0) ex = 0.0; // max(0, S-K)
                } else {
                    ex = K - Sij;
                    if (ex < 0.0) ex = 0.0; // max(0, K-S)
                }
                V[j] = std::max(cont, ex);
            } else { //if European option, just take the continuation value
                V[j] = cont;
            }
        }
    }
    return V[0];
}
