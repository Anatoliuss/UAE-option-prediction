#pragma once                      
#include <cmath>                  
#include <utility>             
#include <corecrt_math_defines.h>
#include <vector>
double binomial_price(
    double S, double K, double r, double sigma, double T,
    int steps, bool is_call, bool american = false, double q = 0.0
){
    const double dt = T / steps;
    const double u  = std::exp(sigma * std::sqrt(dt)); 
    const double d  = 1.0 / u;                  
    const double disc = std::exp(-r * dt);      // discount per step
    const double p = (std::exp((r - q) * dt) - d) / (u - d); // risk-neutral prob

    std::vector<double> V(steps + 1);
    for (int j = 0; j <= steps; ++j) {
        const double ST = S * std::pow(u, steps - j) * std::pow(d, j);
        const double intrinsic = is_call ? std::max(ST - K, 0.0)
                                        : std::max(K - ST, 0.0);
        V[j] = intrinsic;
    }
    for (int i = steps - 1; i >= 0; --i) {
        for (int j = 0; j <= i; ++j) {
            const double cont = disc * (p * V[j] + (1.0 - p) * V[j + 1]);

            if (american) {
                const double S_ij = S * std::pow(u, i - j) * std::pow(d, j);
                const double intrinsic = is_call ? std::max(S_ij - K, 0.0)
                                                : std::max(K - S_ij, 0.0);
                V[j] = std::max(cont, intrinsic);
            } else {
                V[j] = cont;
            }
        }
        if (steps <= 0) return std::numeric_limits<double>::quiet_NaN();
        if (!(p > 0.0 && p < 1.0)) {

            return std::numeric_limits<double>::quiet_NaN();
        }
    }
    return V[0];
}
