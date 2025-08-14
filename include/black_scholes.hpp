#pragma once                      
#include <cmath>                  
#include <utility>             
#include <corecrt_math_defines.h>
//Key idea: It models the price of European call and put options using the Black-Scholes model, which assumes a log-normal distribution of stock prices.
//It computes the option prices, Greeks (delta, vega), and implied volatility.
//Black-Scholes assumption: stock price follows Geometric Brownian Motion (GBM)
//Black-Scholes PDE: ∂C/∂t + 0.5σ²S²∂²C/∂S² + (r-q)S∂C/∂S - rC = 0
//Solving for european type :
//Black-Scholes formula: C(S, K, T) = S * N(d1) - K * e^(-rT) * N(d2)
//Letters meaning:
// S = current stock price
// K = strike price
// T = time to expiration
// r = risk-free interest rate
// σ = volatility of the underlying asset
inline double norm_cdf(double x) {
    return 0.5 * std::erfc(-x * (1/std::sqrt(2)));  
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
    return S * norm_cdf(d1) - K * std::exp(-r * T) * norm_cdf(d2); // using the formula
}

inline double put_price(double S, double K,
                        double r, double sigma, double T) {
    auto [d1, d2] = compute_d1d2(S, K, r, sigma, T);
    return K * std::exp(-r * T) * norm_cdf(-d2) - S * norm_cdf(-d1); //also using the formula
}
// Greeks 
inline double delta_call(double S,double K,double r,double sigma,double T) {
    const double D1 = d1(S,K,r, sigma, T); 
    return norm_cdf(D1);  // slope of the option price curve with respect to S (basically first derivative)
}

inline double delta_put(double S,double K,double r,double sigma,double T) {
    const double D1 = d1(S,K,r,  sigma, T);
    return norm_cdf(D1) - 1.0;  // -(1 - N(d1)), same but for put options

}
// Standard normal PDF
inline double norm_pdf(double x) {
    static constexpr double INV_SQRT_2PI = 0.3989422804014327; // 1/sqrt(2π)
    return INV_SQRT_2PI * std::exp(-0.5 * x * x);
}

// Vega
inline double vega(double S, double K, double r, double sigma, double T) {
    auto [D1, D2] = compute_d1d2(S, K, r, sigma, T);
    return S * norm_pdf(D1) * std::sqrt(T); //second derivative of option price wrt sigma
}
inline double implied_vol_call_bisect(double S,double K,double r,double T,
                                      double market_price,
                                      double lo = 1e-6,
                                      double hi = 5.0,
                                      double tol = 1e-8, // tolerance for convergence
                                      int max_iters = 100); //declaration

inline double implied_vol_call(double S,double K,double r,double T,
                               double market_price,
                               double sigma0 = 0.2,
                               double tol = 1e-8,// tolerance for convergence
                               int max_newton = 20)
{
    // Try Newton
    //Newton's method: iteratively improve guess for sigma until the option price matches the market price
    double sigma = sigma0;
    for (int i = 0; i < max_newton; ++i) {
        const double price = call_price(S,K,r,sigma,T);
        const double diff  = price - market_price;
        if (std::fabs(diff) < tol){
            return sigma;  // found a good sigma
        } 

        const double v = vega(S,K,r,sigma,T);
        if (v < 1e-12 || !std::isfinite(v)) {
            break;
        } // if vega is too small or not finite, stop
        double step = diff / v; 
        sigma -= step;

        if (!std::isfinite(sigma)) {
            break;
        }
        if (sigma < 1e-6 || sigma > 5.0) {
            break; 
        }
        if (std::fabs(step) < tol) {
            return sigma;
        }
    }

    // bisection if fallback needed (sigma didn't converge)
    return implied_vol_call_bisect(S,K,r,T,market_price);
}
//bisection method: finds a root of the function f(sigma) = call_price(S,K,r,sigma,T) - market_price
// It finds a sigma such that the call price matches the market price (basically a brute-force search)
inline double implied_vol_call_bisect(double S,double K,double r,double T,
                                      double market_price,
                                      double lo,
                                      double hi,
                                      double tol,
                                      int max_iters)
{
    // Define a function f(sigma) that gives the difference between the Black-Scholes call price and the market price
    auto f = [&](double sigma) {
        double price = call_price(S, K, r, sigma, T); // Calculate the call price using Black-Scholes
        return price - market_price; 
    };

    double flo = f(lo), fhi = f(hi); // Evaluate f at the bounds
    if (flo == 0.0) return lo;
    if (fhi == 0.0) return hi;
    if (flo * fhi > 0.0) {
        // If the function does not change sign over [lo, hi], return the endpoint with the smaller absolute value
        if (std::fabs(flo) < std::fabs(fhi))
            return lo;
        else
            return hi;
    }

    for (int i = 0; i < max_iters; ++i) {
        double mid = 0.5 * (lo + hi);
        double fmid = f(mid);
        if (std::fabs(fmid) < tol || (hi - lo) * 0.5 < tol) {
            return mid;
        }
        if (flo * fmid <= 0.0) { 
            hi = mid; 
            fhi = fmid;
        } else {
            lo = mid;
            flo = fmid;
        }
    }
    return 0.5 * (lo + hi);
}
// Put implied volatility (similar to call, but for put options)
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
        if (v < 1e-12 || !std::isfinite(v)) {
            break;
        } // if vega is too small or not finite, stop
        const double step = diff / v;
        sigma -= step;

        if (!std::isfinite(sigma)) {
            break;
        }
        if (sigma < 1e-6 || sigma > 5.0) {
            break;
        }
        if (std::fabs(step) < tol) {
            return sigma;
        }
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
    if (flo == 0.0) {
        return lo;
    }
    if (fhi == 0.0) {
        return hi;
    }
    if (flo * fhi > 0.0) {
        if (std::fabs(flo) < std::fabs(fhi)) {
            return lo;
        } else {
            return hi;
        }
    }
    for (int i = 0; i < max_iters; ++i) {
        double mid = 0.5*(lo+hi);
        double fmid = f(mid);
        if (std::fabs(fmid) < tol || (hi-lo)*0.5 < tol) {
            return mid;
        }
        if (flo * fmid <= 0.0) {
            hi = mid;
            fhi = fmid;
        } else {
            lo = mid;
            flo = fmid;
        }
    }
    return 0.5*(lo+hi);
}
