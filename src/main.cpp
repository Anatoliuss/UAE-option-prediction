#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>

#include "black_scholes.hpp"
#include "binomial.hpp"
#include "monte_carlo.hpp"

// Simple struct to hold one option's parameters
struct OptionRow {
    double S, K, r, sigma, T;
    std::string type; // "C" or "P"
};

// Parse CSV into vector<OptionRow>
// CSV format: S,K,r,sigma,T,type
std::vector<OptionRow> load_csv(const std::string& path) {
    std::ifstream file(path);
    if (!file) {
        throw std::runtime_error("Could not open file: " + path);
    }
    std::vector<OptionRow> rows;
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::istringstream ss(line);
        std::string field;
        OptionRow row{};
        std::getline(ss, field, ','); row.S = std::stod(field);
        std::getline(ss, field, ','); row.K = std::stod(field);
        std::getline(ss, field, ','); row.r = std::stod(field);
        std::getline(ss, field, ','); row.sigma = std::stod(field);
        std::getline(ss, field, ','); row.T = std::stod(field);
        std::getline(ss, field, ','); row.type = field;
        rows.push_back(row);
    }
    return rows;
}

int main(int argc, char* argv[]) {
    std::string method;
    std::string file;
    double q = 0.0;
    int steps = 512;
    std::uint64_t npaths = 100000;

    // crude flag parsing
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--method" && i+1 < argc) method = argv[++i];
        else if (arg == "--file" && i+1 < argc) file = argv[++i];
        else if (arg == "--steps" && i+1 < argc) steps = std::stoi(argv[++i]);
        else if (arg == "--npaths" && i+1 < argc) npaths = std::stoull(argv[++i]);
        else if (arg == "--q" && i+1 < argc) q = std::stod(argv[++i]);
        else if (arg == "--help") {
            std::cout << "Usage: app --method bs|bin|mc --file path "
                         "[--steps N] [--npaths M] [--q qyield]\n";
            return 0;
        }
    }

    if (method.empty() || file.empty()) {
        std::cerr << "Error: --method and --file are required. Use --help.\n";
        return 1;
    }

    std::vector<OptionRow> rows;
    try {
        rows = load_csv(file);
    } catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
        return 1;
    }

    for (const auto& opt : rows) {
        bool is_call = (opt.type == "C" || opt.type == "c" ||
                        opt.type == "Call" || opt.type == "call");

        if (method == "bs") {
            double price = is_call
                ? call_price(opt.S, opt.K, opt.r, opt.sigma, opt.T)
                : put_price (opt.S, opt.K, opt.r, opt.sigma, opt.T);
            std::cout << price << "\n";
        }
        else if (method == "bin") {
            double price = binomial_price(opt.S, opt.K, opt.r, opt.sigma,
                                          opt.T, steps, is_call, false, q);
            std::cout << price << "\n";
        }
        else if (method == "mc") {
            auto [price,se] = mc_euro_option_price(opt.S, opt.K, opt.r, q,
                                                   opt.sigma, opt.T, is_call,
                                                   npaths, 42, true);
            std::cout << price << " (SE=" << se << ")\n";
        }
        else {
            std::cerr << "Unknown method: " << method << "\n";
            return 1;
        }
    }
    return 0;
}
