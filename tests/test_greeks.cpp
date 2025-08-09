#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "black_scholes.hpp"

TEST_CASE("Delta (no dividends)") {
    const double S=100, K=100, r=0.05, sigma=0.2, T=1.0;
    const double dC = delta_call(S,K,r,sigma,T);
    const double dP = delta_put (S,K,r,sigma,T);

    REQUIRE(dC == Catch::Approx(0.6368306511756191).epsilon(1e-12));
    REQUIRE(dP == Catch::Approx(-0.3631693488243809).epsilon(1e-12));
}

TEST_CASE("Delta OTM/ITM (no dividends)") {
    const double S=100, K=110, r=0.05, sigma=0.2, T=1.0;
    const double dC = delta_call(S,K,r,sigma,T);
    const double dP = delta_put (S,K,r,sigma,T);

    REQUIRE(dC == Catch::Approx(0.44964793063717595).epsilon(1e-12));
    REQUIRE(dP == Catch::Approx(-0.5503520693628241).epsilon(1e-12));
}
