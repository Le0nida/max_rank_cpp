#include "csvutils.h"

#include <iomanip>
#include <iostream>

std::vector<Point> readCSV(const std::string& filename, int numRecords, int dimensions) {
    std::ifstream file(filename);
    std::vector<Point> data;
    std::string line, word;

    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    // Skip the header line if present
    if (!std::getline(file, line)) {
        throw std::runtime_error("Empty or invalid file: " + filename);
    }

    data.reserve(numRecords);

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::vector<double> row;
        int id = -1;
        bool first = true;

        while (std::getline(ss, word, ',')) {
            try {
                if (first) {
                    id = std::stoi(word); // First column is the ID
                    first = false;
                } else {
                    row.push_back(std::stod(word)); // Remaining columns are coordinates
                }
            } catch (const std::invalid_argument& e) {
                std::cerr << "Invalid argument: " << word
                          << " in file " << filename << std::endl;
                throw;
            } catch (const std::out_of_range& e) {
                std::cerr << "Out of range: " << word
                          << " in file " << filename << std::endl;
                throw;
            }
        }

        // Ensure row matches expected dimensions
        if (id != -1 && (int)row.size() == dimensions) {
            data.emplace_back(row, id);
        } else {
            throw std::runtime_error("Row does not match expected dimensions: " + line);
        }
    }
    return data;
}

std::vector<int> readQuery(const std::string& filename, int numQueries) {
    std::ifstream file(filename);
    std::vector<int> query;
    std::string line;
    query.reserve(numQueries);

    while (getline(file, line)) {
        try {
            query.push_back(stoi(line));
        } catch (const std::invalid_argument& e) {
            std::cerr << "Invalid argument in query file: " << line << std::endl;
            throw;
        } catch (const std::out_of_range& e) {
            std::cerr << "Out of range in query file: " << line << std::endl;
            throw;
        }
    }

    if (query.size() != numQueries) {
        throw std::runtime_error("Query file does not contain the expected number of queries.");
    }

    return query;
}

template <typename T>
void writeCSV(const std::string& filename, const std::vector<std::vector<T>>& data, const std::vector<std::string>& headers) {
    std::ofstream file(filename);

    // Set precision for floating-point values
    if constexpr (std::is_same<T, double>::value) {
        file << std::fixed << std::setprecision(15);
    }

    // Write headers
    for (size_t i = 0; i < headers.size(); ++i) {
        file << headers[i];
        if (i < headers.size() - 1) {
            file << ",";
        }
    }
    file << "\n";

    // Write data
    for (const auto& row : data) {
        if (!row.empty()) {
            file << row[0]; // First element (ID or first column)

            for (size_t j = 1; j < row.size(); ++j) {
                file << "," << row[j]; // Remaining elements
            }
            file << "\n";
        }
    }
}

// Explicit Template Instantiations
template void writeCSV<int>(const std::string&, const std::vector<std::vector<int>>&, const std::vector<std::string>&);
template void writeCSV<double>(const std::string&, const std::vector<std::vector<double>>&, const std::vector<std::string>&);