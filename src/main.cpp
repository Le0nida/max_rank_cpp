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
#include <filesystem>

using namespace std;

/**
 * Global defaults for optional parameters.
 */
int limitHamWeight = 999;
int maxLevelQTree = 99;
int maxCapacityQNode = 10;
int maxNoBinStringToCheck = 999999;
int halfspacesLengthLimit = 21;

/**
 * \brief Parses command-line flags of the form --flag=value
 *        and updates the global variables accordingly.
 *
 * \param argc Number of CLI arguments
 * \param argv Array of argument strings
 */
void parseArgs(int argc, char* argv[]) {
    // Required parameters come first: datafile, numRecords, dimensions, numQueries, queryfile
    // We parse them by position, then parse optional flags from index 6 onward.
    // e.g. --limit-ham-weight=500
    for (int i = 6; i < argc; i++) {
        std::string arg = argv[i];
        // Identify flags of the form "--key=value"
        if (arg.rfind("--", 0) == 0) {
            // Find '=' sign
            size_t eqPos = arg.find('=');
            if (eqPos == std::string::npos) {
                std::cerr << "Ignoring invalid flag: " << arg << std::endl;
                continue;
            }
            std::string key = arg.substr(2, eqPos - 2);    // e.g. "limit-ham-weight"
            std::string val = arg.substr(eqPos + 1);       // e.g. "500"

            try {
                if (key == "limit-ham-weight") {
                    limitHamWeight = std::stoi(val);
                } else if (key == "max-level-qtree") {
                    maxLevelQTree = std::stoi(val);
                } else if (key == "max-capacity-qnode") {
                    maxCapacityQNode = std::stoi(val);
                } else if (key == "max-nobinstring-to-check") {
                    maxNoBinStringToCheck = std::stoi(val);
                } else if (key == "halfspaces-length-limit") {
                    halfspacesLengthLimit = std::stoi(val);
                } else {
                    std::cerr << "Unknown parameter: --" << key << std::endl;
                }
            } catch (const std::invalid_argument&) {
                std::cerr << "Invalid integer value in flag: " << arg << std::endl;
            } catch (const std::out_of_range&) {
                std::cerr << "Out-of-range integer value in flag: " << arg << std::endl;
            }
        } else {
            std::cerr << "Ignoring unrecognized argument: " << arg << std::endl;
        }
    }

    // Validate optional parameters
    if (limitHamWeight < 0 || maxLevelQTree < 1 ||
        maxCapacityQNode < 1 || maxNoBinStringToCheck < 1 ||
        halfspacesLengthLimit < 1)
    {
        throw std::runtime_error("One or more optional parameters are invalid (<=0).");
    }
}

/**
 * \brief (Optional) Reads a config file to set the global parameters.
 *        Format: Each line as key=value (like limitHamWeight=999).
 *
 * \param configFile Path to the configuration file.
 */
void parseConfigFile(const std::string& configFile) {
    // Convert to an absolute path if not already
    std::filesystem::path configPath(configFile);
    if (!configPath.is_absolute()) {
        configPath = std::filesystem::absolute(configPath);
    }

    std::ifstream in(configPath);
    if (!in.is_open()) {
        throw std::runtime_error("Could not open config file: " + configPath.string());
    }
    std::string line;
    while (std::getline(in, line)) {
        // skip empty lines or comments
        if (line.empty() || line[0] == '#') continue;
        size_t eqPos = line.find('=');
        if (eqPos == std::string::npos) continue;

        std::string key = line.substr(0, eqPos);
        std::string val = line.substr(eqPos + 1);

        try {
            if (key == "limitHamWeight") {
                limitHamWeight = std::stoi(val);
            } else if (key == "maxLevelQTree") {
                maxLevelQTree = std::stoi(val);
            } else if (key == "maxCapacityQNode") {
                maxCapacityQNode = std::stoi(val);
            } else if (key == "maxNoBinStringToCheck") {
                maxNoBinStringToCheck = std::stoi(val);
            } else if (key == "halfspacesLengthLimit") {
                halfspacesLengthLimit = std::stoi(val);
            } else {
                std::cerr << "Unknown config key: " << key << std::endl;
            }
        } catch (const std::invalid_argument&) {
            std::cerr << "Invalid integer value for key " << key << " in config file.\n";
        } catch (const std::out_of_range&) {
            std::cerr << "Out-of-range integer value for key " << key << " in config file.\n";
        }
    }

    in.close();

    // Validate again
    if (limitHamWeight < 0 || maxLevelQTree < 1 ||
        maxCapacityQNode < 1 || maxNoBinStringToCheck < 1 ||
        halfspacesLengthLimit < 1)
    {
        throw std::runtime_error("Invalid config file parameter (<=0).");
    }
}

int main(const int argc, char* argv[]) {

    // Start execution timer
    const auto start = std::chrono::high_resolution_clock::now();

    // Check required parameters
    if (argc < 6) {
        std::cerr << "Usage: " << argv[0]
                  << " <datafile> <numRecords> <dimensions> <numQueries> <queryfile>\n"
                  << " [OPTIONAL: 6th param is a config file or CLI flags like --limit-ham-weight=999]\n\n"
                  << "Example flags:\n"
                  << "  --limit-ham-weight=999\n"
                  << "  --max-level-qtree=8\n"
                  << "  --max-capacity-qnode=20\n"
                  << "  --max-nobinstring-to-check=999999\n"
                  << "  --halfspaces-length-limit=21\n"
                  << std::endl;
        return 1;
    }

    // Required parameters
    std::string datafile   = argv[1];
    std::string queryfile  = argv[5];
    int numRecords         = 0;
    int dimensions         = 0;
    int numQueries         = 0;

    try {
        numRecords = std::stoi(argv[2]);
        dimensions = std::stoi(argv[3]);
        numQueries = std::stoi(argv[4]);
        if (numRecords <= 0 || dimensions <= 0 || numQueries <= 0) {
            throw std::invalid_argument("Input numbers must be positive integers.");
        }
    } catch (const std::exception& e) {
        std::cerr << "Invalid input for required parameters: " << e.what() << std::endl;
        return 1;
    }

    // 6th parameter could be a config file or a series of --flag=value arguments
    if (argc >= 7) {  // Ensure there's a 6th argument
        std::string arg6 = argv[6];

        // Convert to absolute path if it looks like a file path
        std::filesystem::path configPath(arg6);
        if (std::filesystem::exists(configPath) && std::filesystem::is_regular_file(configPath)) {
            try {
                parseConfigFile(configPath.string());  // Load config file
            } catch (const std::exception& e) {
                std::cerr << "Error loading config file: " << e.what() << "\n";
            }
        } else if (arg6.rfind("--", 0) == 0) {
            // If it starts with "--", assume it's a CLI argument and parse it
            parseArgs(argc, argv);
        } else {
            std::cerr << "Ignoring unrecognized argument: " << arg6 << std::endl;
        }
    } else if (argc > 6) {
        parseArgs(argc, argv);
    }

    std::cout << "\n========================[ PARAMS SUMMARY ]========================\n";
    std::cout << "Data File:           " << datafile << "\n";
    std::cout << "Query File:          " << queryfile << "\n";
    std::cout << "Records Loaded:      " << numRecords << "\n";
    std::cout << "Dimensions:          " << dimensions << "\n";
    std::cout << "Queries:             " << numQueries << "\n";
    std::cout << "Algorithm Parameters:\n";
    std::cout << "   limitHamWeight:          " << limitHamWeight << "\n";
    std::cout << "   maxLevelQTree:           " << maxLevelQTree << "\n";
    std::cout << "   maxCapacityQNode:        " << maxCapacityQNode << "\n";
    std::cout << "   maxNoBinStringToCheck:   " << maxNoBinStringToCheck << "\n";
    std::cout << "   halfspacesLengthLimit:   " << halfspacesLengthLimit << "\n";

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
