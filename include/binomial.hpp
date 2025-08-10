#pragma once
#include <cmath>
#include <vector>
#include <limits>
#include <algorithm>

inline double binomial_price(
    double S, double K, double r, double sigma, double T,
    int steps, bool is_call, bool american = false, double q = 0.0)
{
    if (steps <= 0 || S <= 0.0 || K <= 0.0 || T <= 0.0 || sigma <= 0.0)
        return std::numeric_limits<double>::quiet_NaN();

    const double dt   = T / steps;
    const double u    = std::exp(sigma * std::sqrt(dt));
    const double d    = 1.0 / u;
    const double disc = std::exp(-r * dt);
    const double p    = (std::exp((r - q) * dt) - d) / (u - d);
    if (!(p > 0.0 && p < 1.0))
        return std::numeric_limits<double>::quiet_NaN();

    // V[j] means node with j UP moves (more standard)
    std::vector<double> V(steps + 1);
    for (int j = 0; j <= steps; ++j) {
        const double ST = S * std::pow(u, j) * std::pow(d, steps - j);
        V[j] = is_call ? std::max(ST - K, 0.0) : std::max(K - ST, 0.0);
    }

    for (int i = steps - 1; i >= 0; --i) {
        for (int j = 0; j <= i; ++j) {
            const double cont = disc * (p * V[j + 1] + (1.0 - p) * V[j]);
            if (american) {
                const double Sij = S * std::pow(u, j) * std::pow(d, i - j);
                const double ex  = is_call ? std::max(Sij - K, 0.0)
                                           : std::max(K - Sij, 0.0);
                V[j] = std::max(cont, ex);
            } else {
                V[j] = cont;
            }
        }
    }
    return V[0];
}
