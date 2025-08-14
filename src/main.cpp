#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <stdexcept>
#include <cctype>

#include "black_scholes.hpp"
#include "binomial.hpp"
#include "monte_carlo.hpp"

struct OptionRow {
    double S, K, r, sigma, T, q;
    std::string type; // "C" or "P"
    std::string date; // optional
};

static inline std::string trim(const std::string& s) {
    size_t a = 0, b = s.size();
    while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
    while (b > a && std::isspace(static_cast<unsigned char>(s[b-1]))) --b;
    return s.substr(a, b - a);
}
static inline std::string lower(std::string s) {
    for (char &c : s) c = (char)std::tolower((unsigned char)c);
    return s;
}

std::vector<OptionRow> load_csv_headered(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("Could not open CSV: " + path);

    std::string header;
    if (!std::getline(f, header)) throw std::runtime_error("Empty CSV: " + path);

    std::vector<std::string> cols;
    {
        std::istringstream hs(header);
        std::string c;
        while (std::getline(hs, c, ',')) cols.push_back(trim(c));
    }
    std::unordered_map<std::string,int> idx;
    for (int i=0;i<(int)cols.size();++i) idx[lower(cols[i])] = i;

    auto need = [&](const char* name)->int{
        auto it = idx.find(lower(name));
        if (it==idx.end()) throw std::runtime_error(std::string("Missing required column: ")+name);
        return it->second;
    };

    const int iS     = need("S");
    const int iK     = need("K");
    const int ir     = need("r");
    const int isigma = need("sigma");
    const int iT     = need("T");
    const int itype  = need("type");
    const int iq     = idx.count("q") ? idx["q"] : -1;
    const int idate  = idx.count("date") ? idx["date"] : -1;

    std::vector<OptionRow> out;
    std::string line;
    while (std::getline(f, line)) {
        if (trim(line).empty()) continue;
        std::vector<std::string> fields;
        std::istringstream ls(line);
        std::string cell;
        while (std::getline(ls, cell, ',')) fields.push_back(trim(cell));
        if ((int)fields.size() < (int)cols.size()) fields.resize(cols.size());

        OptionRow r{};
        r.S     = std::stod(fields[iS]);
        r.K     = std::stod(fields[iK]);
        r.r     = std::stod(fields[ir]);
        r.sigma = std::stod(fields[isigma]);
        r.T     = std::stod(fields[iT]);
        r.type  = fields[itype];
        r.q     = (iq>=0 && !fields[iq].empty()) ? std::stod(fields[iq]) : 0.0;
        r.date  = (idate>=0 ? fields[idate] : "");
        out.push_back(r);
    }
    return out;
}

// ----- main CLI -----
int main(int argc, char* argv[]) {
    std::string method;
    std::string file;
    double q_override = std::numeric_limits<double>::quiet_NaN();
    int steps = 512;
    std::uint64_t npaths = 100000;

    for (int i=1;i<argc;++i){
        std::string a = argv[i];
        if (a=="--method" && i+1<argc) method = argv[++i];
        else if (a=="--file" && i+1<argc) file = argv[++i];
        else if (a=="--steps" && i+1<argc) steps = std::stoi(argv[++i]);
        else if (a=="--npaths" && i+1<argc) npaths = std::stoull(argv[++i]);
        else if (a=="--q" && i+1<argc) q_override = std::stod(argv[++i]);
        else if (a=="--help") {
            std::cout << "Usage: app --method bs|bin|mc --file <csv> "
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
        rows = load_csv_headered(file);
    } catch (const std::exception& e) {
        std::cerr << "[load error] " << e.what() << "\n";
        return 1;
    }

    for (const auto& opt : rows) {
        const bool is_call = (opt.type=="C" || opt.type=="c");
        const double q = std::isnan(q_override) ? opt.q : q_override;

        if (method=="bs") {
            double price = is_call
                ? call_price(opt.S,opt.K,opt.r,opt.sigma,opt.T)
                : put_price (opt.S,opt.K,opt.r,opt.sigma,opt.T);
            std::cout << price << "\n";
        } else if (method=="bin") {
            double price = binomial_price(opt.S,opt.K,opt.r,opt.sigma,opt.T,
                                          steps,is_call,false,q);
            std::cout << price << "\n";
        } else if (method=="mc") {
            auto [price,se] = mc_euro_option_price(opt.S,opt.K,opt.r,q,
                                                   opt.sigma,opt.T,is_call,
                                                   npaths,42,true);
            std::cout << price << " (SE=" << se << ")\n";
        } else {
            std::cerr << "Unknown method: " << method << "\n";
            return 2;
        }
    }
    return 0;
}
