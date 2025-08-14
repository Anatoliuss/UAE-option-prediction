#include <iostream>      
#include <fstream>       
#include <sstream>       
#include <string>        
#include <unordered_map> 
#include <stdexcept>     
#include <cctype>        
#include <iomanip>      
#include <limits>    
#include "black_scholes.hpp"  
#include "binomial.hpp"      
#include "monte_carlo.hpp"  

// A single row from the input CSV describing one option contract.
struct OptionRow {
    double S     = 0.0; // underlying spot price
    double K     = 0.0; // strike price
    double r     = 0.0; // risk-free rate
    double sigma = 0.0; // volatility
    double T     = 0.0; // time to maturity in years
    double q     = 0.0; // dividend yield (optional; defaults to 0)
    std::string type;   // "C" for call, "P" for put
    std::string date; 
};
struct PrintOpts {
    bool pretty = false;       
    std::string out_csv = ""; 
};

// trim the string
static inline std::string trim(const std::string& s) {
    size_t a = 0, b = s.size();
    while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
    while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
    return s.substr(a, b - a);
}

// make a string lowercase
static inline std::string lower(std::string s) {
    for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

// Read a headered CSV file and parse each subsequent line into an OptionRow.
static std::vector<OptionRow> load_csv_headered(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("Could not open CSV: " + path);

    // Read the header line and split by commas to get column names
    std::string header;
    if (!std::getline(f, header)) throw std::runtime_error("Empty CSV: " + path);

    std::vector<std::string> cols;
    {
        std::istringstream hs(header);
        std::string c;
        while (std::getline(hs, c, ',')) cols.push_back(trim(c));
    }
    std::unordered_map<std::string, int> idx;
    for (int i = 0; i < static_cast<int>(cols.size()); ++i) idx[lower(cols[i])] = i;

    // Helper to get required column index (
    auto need = [&](const char* name) -> int {
        auto it = idx.find(lower(name));
        if (it == idx.end()) throw std::runtime_error(std::string("Missing required column: ") + name);
        return it->second;
    };

    //Find indices for required and optional columns
    const int iS     = need("S");
    const int iK     = need("K");
    const int ir     = need("r");
    const int isigma = need("sigma");
    const int iT     = need("T");
    const int itype  = need("type");
    const int iq     = idx.count("q") ? idx["q"] : -1;         // optional
    const int idate  = idx.count("date") ? idx["date"] : -1;   // optional

    // read each data line, split by comma
    std::vector<OptionRow> out;
    std::string line;
    while (std::getline(f, line)) {
        if (trim(line).empty()) continue; // skip blank lines

        std::vector<std::string> fields;
        std::istringstream ls(line);
        std::string cell;
        while (std::getline(ls, cell, ',')) fields.push_back(trim(cell));
        if (static_cast<int>(fields.size()) < static_cast<int>(cols.size())) fields.resize(cols.size());

        OptionRow r{};
        r.S     = std::stod(fields[iS]);
        r.K     = std::stod(fields[iK]);
        r.r     = std::stod(fields[ir]);
        r.sigma = std::stod(fields[isigma]);
        r.T     = std::stod(fields[iT]);
        r.type  = fields[itype];
        r.q     = (iq >= 0 && !fields[iq].empty()) ? std::stod(fields[iq]) : 0.0;
        r.date  = (idate >= 0 ? fields[idate] : "");
        out.push_back(r);
    }
    return out;
}

// If pretty-printing is enabled, print a CSV header once to stdout
static void print_header_if_needed(const PrintOpts& po, bool is_mc) {
    if (!po.pretty) return;
    if (is_mc) {
        std::cout << "date,S,K,T,type,price,se\n";
    } else {
        std::cout << "date,S,K,T,type,price\n";
    }
}

// Print one result row in a friendly CSV format to stdout
static void print_row_pretty(const OptionRow& opt, double price, double se, bool is_mc) {
    std::cout << std::fixed << std::setprecision(4); //for readability
    std::cout << (opt.date.empty() ? "-" : opt.date) << ","
              << opt.S << "," << opt.K << "," << opt.T << ","
              << opt.type << "," << price;
    if (is_mc) std::cout << "," << se;
    std::cout << "\n";
}

int main(int argc, char* argv[]) {
    std::string method;   // which pricing method to use: "bs", "bin", or "mc"
    std::string file;     // path to input CSV
    double q_override = std::numeric_limits<double>::quiet_NaN(); // override all rows if provided
    int steps = 512;      // binomial tree steps 
    std::uint64_t npaths = 100000; // Monte Carlo paths
    PrintOpts po;         // output formatting options

    int i = 1;
    while (i < argc) {
        std::string a = argv[i];
        if (a == "--method" && i + 1 < argc) { method = argv[++i]; }
        else if (a == "--file" && i + 1 < argc) { file = argv[++i]; }
        else if (a == "--steps" && i + 1 < argc) { steps = std::stoi(argv[++i]); }
        else if (a == "--npaths" && i + 1 < argc) { npaths = std::stoull(argv[++i]); }
        else if (a == "--q" && i + 1 < argc) { q_override = std::stod(argv[++i]); }
        else if (a == "--pretty") { po.pretty = true; }
        else if (a == "--out" && i + 1 < argc) { po.out_csv = argv[++i]; }
        else if (a == "--help") {
            std::cout << "Usage: app --method bs|bin|mc --file <csv> "
                         "[--steps N] [--npaths M] [--q qyield] "
                         "[--pretty] [--out output.csv]\n";
            return 0;
        }
        ++i;
    }

    if (method.empty() || file.empty()) {
        std::cerr << "Error: --method and --file are required. Use --help.\n";//error handling
        return 1;
    }
    std::ofstream csv_out;
    bool write_csv = !po.out_csv.empty();
    if (write_csv) {
        csv_out.open(po.out_csv);
        if (!csv_out) {
            // If folder doesn't exist or path is invalid, keep going without file output
            std::cerr << "[warn] could not open --out " << po.out_csv << " for writing; continuing without file.\n";
            write_csv = false;
        } else {
            // Write header and set consistent numeric formatting for file
            if (method == "mc") csv_out << "date,S,K,T,type,price,se\n";
            else                 csv_out << "date,S,K,T,type,price\n";
            csv_out << std::fixed << std::setprecision(6);
        }
    }
    std::vector<OptionRow> rows;
    try {
        rows = load_csv_headered(file);
    } catch (const std::exception& e) {
        std::cerr << "[load error] " << e.what() << "\n";
        return 1;
    }
    for (std::size_t rowIndex = 0; rowIndex < rows.size(); ++rowIndex) {
        const OptionRow& opt = rows[rowIndex];
        const bool is_call = (opt.type == "C" || opt.type == "c");
        const double q_used = std::isnan(q_override) ? opt.q : q_override;

        static bool header_printed = false;
        if (po.pretty) {
            if (!header_printed) {
                print_header_if_needed(po, method == "mc");
                header_printed = true;
            }
        }

        if (method == "bs") {
            // Black–Scholes: closed-form formula (q is ignored in current implementation)
            const double price = is_call
                ? call_price(opt.S, opt.K, opt.r, opt.sigma, opt.T)
                : put_price (opt.S, opt.K, opt.r, opt.sigma, opt.T);

            if (po.pretty) {
                print_row_pretty(opt, price, /*se=*/0.0, /*is_mc=*/false);
            } else {
                std::cout << std::setprecision(8) << price << "\n";
            }
            if (write_csv) {
                csv_out << (opt.date.empty() ? "-" : opt.date) << ","
                        << opt.S << "," << opt.K << "," << opt.T << ","
                        << opt.type << "," << price << "\n";
            }
        }
        else if (method == "bin") {
            // Binomial tree: uses dividend yield (q_used)
            const double price = binomial_price(
                opt.S, opt.K, opt.r, opt.sigma, opt.T,
                steps, is_call, /*american=*/false, q_used
            );

            if (po.pretty) {
                print_row_pretty(opt, price, /*se=*/0.0, /*is_mc=*/false);
            } else {
                std::cout << std::setprecision(8) << price << "\n";
            }
            if (write_csv) {
                csv_out << (opt.date.empty() ? "-" : opt.date) << ","
                        << opt.S << "," << opt.K << "," << opt.T << ","
                        << opt.type << "," << price << "\n";
            }
        }
        else if (method == "mc") {
            // Monte Carlo: returns price and an estimated standard error
            const std::pair<double,double> mc_result = mc_euro_option_price(
                opt.S, opt.K, opt.r, q_used, opt.sigma, opt.T,
                is_call, npaths, /*seed=*/42, /*antithetic=*/true
            );
            const double price = mc_result.first;
            const double se    = mc_result.second;

            if (po.pretty) {
                print_row_pretty(opt, price, se, /*is_mc=*/true);
            } else {
                std::cout << std::fixed << std::setprecision(6)
                          << price << " (SE=" << se << ")\n";
            }
            if (write_csv) {
                csv_out << (opt.date.empty() ? "-" : opt.date) << ","
                        << opt.S << "," << opt.K << "," << opt.T << ","
                        << opt.type << "," << price << "," << se << "\n";
            }
        }
    }

    return 0; 
}
