#pragma once                      
#include <cmath>                  
#include <utility>             
#include <corecrt_math_defines.h>
inline double norm_cdf(double x) {
    return 0.5 * std::erfc(-x * M_SQRT1_2);  
}
// Computes d1 and d2 for Black-Scholes formula
// d1 = (ln(S/K) + (r + σ²/2)T) / (σ√T)
// d2 = d1 - σ√T
struct D1D2 {
    double d1;
    double d2;
};

inline D1D2 compute_d1d2(double S, double K,
                         double r, double sigma, double T) {
    const double sigma_sqrtT = sigma * std::sqrt(T);             // σ√T
    const double num = std::log(S / K) + (r + 0.5 * sigma * sigma) * T;
    const double d1 = num / sigma_sqrtT;
    const double d2 = d1 - sigma_sqrtT;                          // d2 = d1 − σ√T
    return {d1, d2};
}

// ---------- Public pricing API ----------
inline double call_price(double S, double K,
                         double r, double sigma, double T) {
    auto [d1, d2] = compute_d1d2(S, K, r, sigma, T);
    return S * norm_cdf(d1) - K * std::exp(-r * T) * norm_cdf(d2);
}

inline double put_price(double S, double K,
                        double r, double sigma, double T) {
    auto [d1, d2] = compute_d1d2(S, K, r, sigma, T);
    return K * std::exp(-r * T) * norm_cdf(-d2) - S * norm_cdf(-d1);
}
