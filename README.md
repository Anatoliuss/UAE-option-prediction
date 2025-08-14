# **Quantitative Finance Option Pricing Toolkit**

## **Overview**

This project implements three foundational option pricing approaches:

1. **Black–Scholes–Merton Model** (closed-form, European options)
2. **Binomial Tree Model** (discrete-time lattice)
3. **Monte Carlo Simulation** (path-based stochastic pricing)

It reads option parameters from a CSV file and outputs option prices (and standard errors for Monte Carlo) using the selected method.

**Key features**:

* European call and put pricing
* Continuous dividend yield support (binomial, Monte Carlo)
* Implied volatility estimation (Newton + bisection fallback)
* Configurable steps for binomial, paths for Monte Carlo
## **Mathematical Models**

### 1. Black–Scholes–Merton (BSM) Model

Closed-form solution for European options assuming:

* Underlying follows **Geometric Brownian Motion (GBM)**
* Constant volatility and interest rate
* Continuous trading, no arbitrage

**References:**

* Black, F., & Scholes, M. (1973). *The Pricing of Options and Corporate Liabilities*. Journal of Political Economy, 81(3).
* Merton, R.C. (1973). *Theory of Rational Option Pricing*. Bell Journal of Economics and Management Science, 4(1).

###2. Binomial Options Pricing Model

Discrete-time model where the underlying price moves up ($u$) or down ($d$) each step:

The model builds a recombining price tree and works backwards from terminal payoffs to present value.

**References:**

* Cox, J.C., Ross, S.A., & Rubinstein, M. (1979). *Option Pricing: A Simplified Approach*. Journal of Financial Economics, 7(3).

---

###3. Monte Carlo Simulation

Simulates many possible future prices using GBM:
The option price is the discounted mean payoff. Variance reduction via antithetic variates is supported.

**References:**

* Boyle, P.P. (1977). *Options: A Monte Carlo Approach*. Journal of Financial Economics, 4(3).
* Glasserman, P. (2003). *Monte Carlo Methods in Financial Engineering*. Springer.

Data for the project was taken from Yahoo! Finance
