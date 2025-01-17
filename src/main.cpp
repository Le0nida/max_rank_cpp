#include <iostream>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <vector>
#include <string>
#include <Eigen/Dense>
#include "geom.h"
#include "maxrank.h"
#include "qtree.h"
#include "cell.h"
#include <chrono>
#include <csvutils.h>

using namespace std;

// Default parameters
int limitHamWeight = 999;
int maxLevelQTree = 99;
int maxCapacityQNode = 10;
int maxNoBinStringToCheck = 999999;

int main(const int argc, char* argv[]) {

    // Start execution timer
    const auto start = std::chrono::high_resolution_clock::now();

    // Check required parameters
    if (argc < 6) {
        std::cerr << "Usage: " << argv[0]
                  << " <datafile> <numRecords> <dimensions> <numQueries> <queryfile>\n"
                  << " [limitHamWeight] [maxLevelQTree] [maxCapacityQNode] [maxNoBinStringToCheck]\n\n"
                  << "   datafile             = CSV file containing data\n"
                  << "   numRecords           = Number of expected records\n"
                  << "   dimensions           = Number of dimensions in the dataset\n"
                  << "   numQueries           = Number of queries\n"
                  << "   queryfile            = File containing query list\n\n"
                  << "   limitHamWeight       (optional, default=999)\n"
                  << "   maxLevelQTree        (optional, default=8)\n"
                  << "   maxCapacityQNode     (optional, default=20)\n"
                  << "   maxNoBinStringToCheck (optional, default=999999)\n"
                  << std::endl;
        return 1;
    }

    // Read required parameters
    const string datafile    = argv[1];
    const string queryfile   = argv[5];
    int numRecords;
    int dimensions;
    int numQueries;

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

    // Read optional parameters if provided
    if (argc > 6) limitHamWeight       = stoi(argv[6]);
    if (argc > 7) maxLevelQTree    = stoi(argv[7]);
    if (argc > 8) maxCapacityQNode = stoi(argv[8]);
    if (argc > 9) maxNoBinStringToCheck  = stoi(argv[9]);

    // Validate optional parameters
    if (limitHamWeight < 0 || maxLevelQTree < 1 || maxCapacityQNode < 1 || maxNoBinStringToCheck < 1)
    {
        cerr << "One or more input parameters are invalid (<=0)." << endl;
        return 1;
    }

    // Load dataset
    vector<Point> data = readCSV(datafile, numRecords, dimensions);
    cout << "Loaded " << data.size() << " records from " << datafile << endl;

    // Load query indices
    vector<int> query = readQuery(queryfile, numQueries);
    cout << "Loaded " << query.size() << " queries from " << queryfile << endl;

    // Main MaxRank routine
    vector<vector<int>> res;
    res.reserve(query.size());
    vector<vector<double>> cells;
    cells.reserve(query.size());

    if (dimensions > 2) {
        for (const int q : query) {
            cout << "#  Processing data point " << q << "  #" << endl;
            const int idx = q - 1;

            cout << "#  " << Eigen::Map<Eigen::VectorXd>(data[idx].coord.data(), data[idx].coord.size()).transpose() << "  #" << endl;

            int maxrank;
            vector<Cell> mincells;
            tie(maxrank, mincells) = aa_hd(data, data[idx]);

            cout << "#  MaxRank: " << maxrank << "  NOfMincells: " << mincells.size() << "  #" << endl;

            // Saving results
            res.push_back({q, maxrank});
            vector cell_entry = { static_cast<double>(q) };
            for (const auto &cell : mincells)
            {
                cell_entry.insert(cell_entry.end(), cell.feasible_pnt.coord.begin(), cell.feasible_pnt.coord.end());
                cell_entry.push_back(1 - accumulate(cell.feasible_pnt.coord.begin(), cell.feasible_pnt.coord.end(), 0.0));

                break; // todo rimuovere e considerarli tutti
            }
            cells.push_back(cell_entry);
        }
    } else {
        for (const int q : query) {
            cout << "#  Processing data point " << q << "  #" << endl;
            const int idx = q - 1;

            cout << "#  " << Eigen::Map<Eigen::VectorXd>(data[idx].coord.data(), data[idx].coord.size()).transpose() << "  #" << endl;

            int maxrank;
            vector<Interval> mincells;
            tie(maxrank, mincells) = aa_2d(data, data[idx]);

            cout << "#  MaxRank: " << maxrank << "  NOfMincells: " << mincells.size() << "  #" << endl;

            // Saving results
            res.push_back({q, maxrank});
            vector<double> cell_entry = { static_cast<double>(q) };
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

    // Print execution time
    const auto end = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> elapsed = end - start;
    cout << "Total execution time: " << elapsed.count() << " seconds." << endl;
    return 0;
}
