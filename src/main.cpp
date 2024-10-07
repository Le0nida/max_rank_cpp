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

using namespace std;

class Cell;
vector<Point> readCSV(const string& filename) {
    ifstream file(filename);
    vector<Point> data;
    string line, word;

    // Check if the file opened successfully
    if (!file.is_open()) {
        throw runtime_error("Could not open file: " + filename);
    }

    // Skip the header line
    if (getline(file, line)) {
        // Successfully read the header line, do nothing
    } else {
        // Failed to read the header line, possibly empty file
        throw runtime_error("Empty or invalid file: " + filename);
    }

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

        if (id != -1 && !row.empty()) {
            data.emplace_back(row, id);
        }
    }

    return data;
}

vector<int> readQuery(const string& filename) {
    ifstream file(filename);
    vector<int> query;
    string line;

    while (getline(file, line)) {
        query.push_back(stoi(line));
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
    if (argc < 4) {
        cerr << "Usage: " << argv[0] << " <datafile> <queryfile> <method>" << endl;
        return 1;
    }

    string datafile = argv[1];
    string queryfile = argv[2];
    string method = argv[3];

    // Load data
    vector<Point> data = readCSV(datafile);
    cout << "Loaded " << data.size() << " records from " << datafile << endl;

    // Load query
    auto query = readQuery(queryfile);
    cout << "Loaded " << query.size() << " queries from " << queryfile << endl;

    // Main MaxRank routine
    vector<vector<double>> res;
    vector<vector<double>> cells;

    if (data[0].dims > 2) {
        for (int q : query) {
            // Reset global variables
            halfspaceCache = nullptr;
            pointToHalfSpaceCache.clear();
            globalNodeID = 0;
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
        /*for (int q : query) {
            cout << "#  Processing data point " << q << "  #" << endl;
            int idx = q - 1;

            cout << "#  " << Eigen::Map<Eigen::VectorXd>(data[idx].coord.data(), data[idx].coord.size()).transpose() << "  #" << endl;

            int maxrank;
            vector<Interval> mincells;
            tie(maxrank, mincells) = aa_2d(data, data[idx]);
            cout << "#  MaxRank: " << maxrank << "  NOfMincells: " << mincells.size() << "  #" << endl;

            res.push_back({(double)q, (double)maxrank});
            vector<double> cell_entry = { (double)q };
            for (const auto& cell : mincells) {
                cell_entry.push_back(cell.range.first);
                cell_entry.push_back(cell.range.second);
            }
            cells.push_back(cell_entry);
        }*/
    }

    // Write results to CSV
    writeCSV(R"(C:\Users\leona\Desktop\maxrank.csv)", res, { "id", "maxrank" });
    writeCSV(R"(C:\Users\leona\Desktop\cells.csv)", cells, { "id", "query_found" });

    return 0;
}
