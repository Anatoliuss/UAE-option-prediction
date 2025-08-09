#include <iostream>
#include "binomial.hpp"
int main() {
    double S=100, K=100, r=0.05, sigma=0.2, T=1.0;
    for (int N : {1, 5, 25, 100, 500}) {
        double price = binomial_price(S,K,r,sigma,T,N, /*is_call=*/true, /*american=*/false, /*q=*/0.0);
        std::cout << "N=" << N << "  price=" << price << "\n";
    }
}
