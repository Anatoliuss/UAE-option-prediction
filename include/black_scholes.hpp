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
inline double d1(double S,double K,double r,double sigma,double T){
    return compute_d1d2(S,K,r,sigma,T).d1;
}
inline double d2(double S,double K,double r,double sigma,double T){
    return compute_d1d2(S,K,r,sigma,T).d2;
}
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
// Greeks (no dividends)
inline double delta_call(double S,double K,double r,double sigma,double T) {
    const double D1 = d1(S,K,r, sigma, T); 
    return norm_cdf(D1);  // N(d1)
}

inline double delta_put(double S,double K,double r,double sigma,double T) {
    const double D1 = d1(S,K,r,  sigma, T);
    return norm_cdf(D1) - 1.0;  // -(1 - N(d1))
}
// Standard normal PDF
inline double norm_pdf(double x) {
    static constexpr double INV_SQRT_2PI = 0.3989422804014327; // 1/sqrt(2π)
    return INV_SQRT_2PI * std::exp(-0.5 * x * x);
}

// Vega (no dividends). Same for calls & puts.
inline double vega(double S, double K, double r, double sigma, double T) {
    auto [D1, D2] = compute_d1d2(S, K, r, sigma, T);
    return S * norm_pdf(D1) * std::sqrt(T);
}
inline double implied_vol_call_bisect(double S,double K,double r,double T,
                                      double market_price,
                                      double lo = 1e-6,
                                      double hi = 5.0,
                                      double tol = 1e-8,
                                      int max_iters = 100);

inline double implied_vol_call(double S,double K,double r,double T,
                               double market_price,
                               double sigma0 = 0.2,
                               double tol = 1e-8,
                               int max_newton = 20)
{
    // Try Newton
    double sigma = sigma0;
    for (int i = 0; i < max_newton; ++i) {
        const double price = call_price(S,K,r,sigma,T);
        const double diff  = price - market_price;
        if (std::fabs(diff) < tol) return sigma;

        const double v = vega(S,K,r,sigma,T);
        if (v < 1e-12 || !std::isfinite(v)) break;

        double step = diff / v;
        sigma -= step;

        if (!std::isfinite(sigma)) break;
        if (sigma < 1e-6 || sigma > 5.0) break; 
        if (std::fabs(step) < tol) return sigma;
    }

    // bisection if fallback needed
    return implied_vol_call_bisect(S,K,r,T,market_price);
}

inline double implied_vol_call_bisect(double S,double K,double r,double T,
                                      double market_price,
                                      double lo,
                                      double hi,
                                      double tol,
                                      int max_iters)
{
    auto f = [&](double sig){ return call_price(S,K,r,sig,T) - market_price; };

    double flo = f(lo), fhi = f(hi);
    if (flo == 0.0) return lo;
    if (fhi == 0.0) return hi;
    if (flo * fhi > 0.0) {
        return (std::fabs(flo) < std::fabs(fhi)) ? lo : hi;
    }

    for (int i = 0; i < max_iters; ++i) {
        double mid = 0.5 * (lo + hi);
        double fmid = f(mid);
        if (std::fabs(fmid) < tol || (hi - lo) * 0.5 < tol) return mid;
        if (flo * fmid <= 0.0) { hi = mid; fhi = fmid; }
        else { lo = mid; flo = fmid; }
    }
    return 0.5 * (lo + hi);
}
inline double implied_vol_put_bisect(double S,double K,double r,double T,
                                     double market_price,
                                     double lo = 1e-6,
                                     double hi = 5.0,
                                     double tol = 1e-8,
                                     int max_iters = 100);

inline double implied_vol_put(double S,double K,double r,double T,
                              double market_price,
                              double sigma0 = 0.2,
                              double tol = 1e-8,
                              int max_newton = 20)
{
    double sigma = sigma0;
    for (int i = 0; i < max_newton; ++i) {
        const double price = put_price(S,K,r,sigma,T);
        const double diff  = price - market_price;   // f(sigma)
        if (std::fabs(diff) < tol) return sigma;

        const double v = vega(S,K,r,sigma,T);
        if (v < 1e-12 || !std::isfinite(v)) break;

        const double step = diff / v;
        sigma -= step;

        if (!std::isfinite(sigma)) break;
        if (sigma < 1e-6 || sigma > 5.0) break;
        if (std::fabs(step) < tol) return sigma;
    }
    return implied_vol_put_bisect(S,K,r,T,market_price);
}

inline double implied_vol_put_bisect(double S,double K,double r,double T,
                                     double market_price,
                                     double lo,
                                     double hi,
                                     double tol,
                                     int max_iters)
{
    auto f = [&](double sig){ return put_price(S,K,r,sig,T) - market_price; };

    double flo = f(lo), fhi = f(hi);
    if (flo == 0.0) return lo;
    if (fhi == 0.0) return hi;
    if (flo * fhi > 0.0) {
        return (std::fabs(flo) < std::fabs(fhi)) ? lo : hi;
    }
    for (int i = 0; i < max_iters; ++i) {
        double mid = 0.5*(lo+hi);
        double fmid = f(mid);
        if (std::fabs(fmid) < tol || (hi-lo)*0.5 < tol) return mid;
        if (flo * fmid <= 0.0) { hi = mid; fhi = fmid; }
        else { lo = mid; flo = fmid; }
    }
    return 0.5*(lo+hi);
}
