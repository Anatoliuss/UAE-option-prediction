#pragma once
#include <string>
#include <sstream>

struct CsvRow {
    std::string date;
    double close = 0.0;
    // Parses a CSV line and returns a CsvRow object
    static CsvRow from_line(const std::string& line) {
        std::stringstream ss(line);
        std::string token;
        CsvRow row;

        std::getline(ss, row.date, ',');          // Date
        std::getline(ss, token, ',');             // Open (skip)
        std::getline(ss, token, ',');             // High (skip)
        std::getline(ss, token, ',');             // Low  (skip)
        std::getline(ss, token, ',');             // Close
        row.close = std::stod(token);

        return row;
    }
};
