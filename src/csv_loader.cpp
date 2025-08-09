#include "CsvRow.hpp"
#include <fstream>
#include <vector>

std::vector<CsvRow> load_csv(const std::string& path) {
    std::ifstream f(path);
    std::vector<CsvRow> rows;
    std::string line;

    std::getline(f, line);            // skip header
    while (std::getline(f, line)) {
        if (!line.empty())
            rows.push_back(CsvRow::from_line(line));
    }
    return rows;
}
