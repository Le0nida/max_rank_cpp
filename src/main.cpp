//
// Created by leona on 06/08/2024.
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <vector>
#include <string>
#include <Eigen/Dense>
#include "geom.h"
#include "maxrank.h"
#include "qtree.h"
#include "query.h"
#include "cell.h"
#include <chrono>

using namespace std;

class Cell;
vector<Point> readCSV(const string& filename, int numRecords, int dimensions) {
    ifstream file(filename);
    vector<Point> data;
    string line, word;

    // Check if the file opened successfully
    if (!file.is_open()) {
        throw runtime_error("Could not open file: " + filename);
    }

    // Skip the header line
    if (!getline(file, line)) {
        throw runtime_error("Empty or invalid file: " + filename);
    }

    // Preallocate memory for data
    data.reserve(numRecords);

    while (getline(file, line)) {
        stringstream ss(line);
        vector<double> row;
        int id = -1;
        bool first = true;

        while (getline(ss, word, ',')) {
            try {
                if (first) {
                    id = stoi(word); // First value is the ID
                    first = false;
                } else {
                    row.push_back(stod(word)); // Rest are coordinates
                }
            } catch (const invalid_argument& e) {
                cerr << "Invalid argument: " << word << " in file " << filename << endl;
                throw;
            } catch (const out_of_range& e) {
                cerr << "Out of range: " << word << " in file " << filename << endl;
                throw;
            }
        }

        if (id != -1 && row.size() == dimensions) {
            data.emplace_back(row, id);
        } else {
            throw runtime_error("Row does not match expected dimensions: " + line);
        }
    }

    return data;
}

vector<int> readQuery(const string& filename, int numQueries) {
    ifstream file(filename);
    vector<int> query;
    string line;

    // Preallocate memory for queries
    query.reserve(numQueries);

    while (getline(file, line)) {
        try {
            query.push_back(stoi(line));
        } catch (const invalid_argument& e) {
            cerr << "Invalid argument in query file: " << line << endl;
            throw;
        } catch (const out_of_range& e) {
            cerr << "Out of range in query file: " << line << endl;
            throw;
        }
    }

    if (query.size() != numQueries) {
        throw runtime_error("Query file does not contain the expected number of queries.");
    }

    return query;
}

void writeCSV(const std::string& filename, const std::vector<std::vector<double>>& data, const std::vector<std::string>& headers) {
    std::ofstream file(filename);

    // Set the desired precision
    file << std::fixed << std::setprecision(15); // Adjust precision as needed

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
            // Write the id as an integer
            int id = static_cast<int>(row[0]);
            file << id;

            // Write the remaining values as doubles with high precision
            for (size_t j = 1; j < row.size(); ++j) {
                file << "," << row[j];
            }
            file << "\n"; // No extra comma at the end
        }
    }
}

int main(int argc, char* argv[]) {
    auto start = std::chrono::high_resolution_clock::now();

    if (argc < 6) {
        cerr << "Usage: " << argv[0] << " <datafile> <numRecords> <dimensions> <numQueries> <queryfile> <method>" << endl;
        return 1;
    }

    string datafile = argv[1];
    int numRecords;
    int dimensions;
    int numQueries;
    string queryfile = argv[5];
    string method = argv[6];


    try {
        numRecords = stoi(argv[2]);
        dimensions = stoi(argv[3]);
        numQueries = stoi(argv[4]);

        if (numRecords <= 0 || dimensions <= 0 || numQueries <= 0) {
            throw invalid_argument("Input numbers must be positive integers.");
        }
    } catch (const invalid_argument& e) {
        cerr << "Invalid input: " << e.what() << endl;
        return 1;
    } catch (const out_of_range& e) {
        cerr << "Input number out of range: " << e.what() << endl;
        return 1;
    }

    // Load data
    vector<Point> data = readCSV(datafile, numRecords, dimensions);
    cout << "Loaded " << data.size() << " records from " << datafile << endl;

    // Load query
    vector<int> query = readQuery(queryfile, numQueries);
    cout << "Loaded " << query.size() << " queries from " << queryfile << endl;

    // Main MaxRank routine
    vector<vector<double>> res;
    res.reserve(query.size());
    vector<vector<double>> cells;
    cells.reserve(query.size());

    if (data[0].dims > 2) {
        for (int q : query) {
            // Reset global variables
            halfspaceCache = nullptr;
            pointToHalfSpaceCache.clear();
            //globalNodeID = 0;
            cout << "#  Processing data point " << q << "  #" << endl;
            int idx = q - 1;  // Assuming query contains 1-based indices

            cout << "#  " << Eigen::Map<Eigen::VectorXd>(data[idx].coord.data(), data[idx].coord.size()).transpose() << "  #" << endl;

            int maxrank;
            vector<Cell> mincells;
            if (method == "BA") {
                //tie(maxrank, mincells) = ba_hd(data, data[idx]);
            } else {
                tie(maxrank, mincells) = aa_hd(data, data[idx]);
            }
            cout << "#  MaxRank: " << maxrank << "  NOfMincells: " << mincells.size() << "  #" << endl;

            res.push_back({(double)q, (double)maxrank});
            vector<double> cell_entry = { (double)q };
            cell_entry.insert(cell_entry.end(), mincells[0].feasible_pnt.coord.begin(), mincells[0].feasible_pnt.coord.end());
            cell_entry.push_back(1 - accumulate(mincells[0].feasible_pnt.coord.begin(), mincells[0].feasible_pnt.coord.end(), 0.0));
            cells.push_back(cell_entry);
        }
    } else {
        for (int q : query) {
            cout << "#  Processing data point " << q << "  #" << endl;
            int idx = q - 1;

            cout << "#  " << Eigen::Map<Eigen::VectorXd>(data[idx].coord.data(),
                    data[idx].coord.size()).transpose() << "  #" << endl;

            int maxrank;
            vector<Interval> mincells;
            tie(maxrank, mincells) = aa_2d(data, data[idx]);

            cout << "#  MaxRank: " << maxrank
                 << "  NOfMincells: " << mincells.size() << "  #" << endl;

            // Salvataggio su CSV
            res.push_back({(double)q, (double)maxrank});

            // Esempio di come salvare i range di ogni mincell
            // (o solo il primo mincell, a seconda delle tue necessità)
            vector<double> cell_entry = { (double)q };
            for (const auto &cell : mincells) {
                cell_entry.push_back(cell.range.first);
                cell_entry.push_back(cell.range.second);
            }
            cells.push_back(cell_entry);
        }
    }

    // Write results to CSV
    writeCSV(R"(C:\Users\leona\Desktop\maxrank.csv)", res, { "id", "maxrank" });
    writeCSV(R"(C:\Users\leona\Desktop\cells.csv)", cells, { "id", "query_found" });

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    cout << "Total execution time: " << elapsed.count() << " seconds." << endl;
    return 0;
}
