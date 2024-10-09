//
// Created by leona on 06/08/2024.
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstdlib> // For malloc, free
#include <cstring> // For memset
#include <Eigen/Dense>
#include "geom.h"
#include "maxrank.h"
#include "qtree.h"
#include "query.h"
#include "cell.h"

using namespace std;

// Assuming Cell class is declared somewhere

Point** readCSV(const string& filename, int& numData) {
    ifstream file(filename);
    Point** data = nullptr;
    numData = 0;
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
        double* coord = nullptr;
        int dims = 0;
        int id = -1;
        bool first = true;

        char* line_cstr = strdup(line.c_str());
        char* token = strtok(line_cstr, ",");
        while (token != nullptr) {
            if (first) {
                id = atoi(token);
                first = false;
            } else {
                coord = (double*)realloc(coord, (dims + 1) * sizeof(double));
                coord[dims] = atof(token);
                dims++;
            }
            token = strtok(nullptr, ",");
        }
        free(line_cstr);

        if (id != -1 && dims > 0) {
            Point* point = new Point(coord, dims, id);
            free(coord);
            numData++;
            data = (Point**)realloc(data, numData * sizeof(Point*));
            data[numData - 1] = point;
        }
    }

    return data;
}

int* readQuery(const string& filename, int& numQuery) {
    ifstream file(filename);
    int* query = nullptr;
    numQuery = 0;
    string line;

    while (getline(file, line)) {
        int qid = stoi(line);
        numQuery++;
        query = (int*)realloc(query, numQuery * sizeof(int));
        query[numQuery - 1] = qid;
    }
    return query;
}

void writeCSV(const std::string& filename, double** data, int numRows, int numCols, const std::string* headers, int numHeaders) {
    std::ofstream file(filename);

    // Set the desired precision
    file << std::fixed << std::setprecision(15); // Adjust precision as needed

    // Write headers
    for (int i = 0; i < numHeaders; ++i) {
        file << headers[i];
        if (i < numHeaders - 1) {
            file << ",";
        }
    }
    file << "\n";

    // Write data
    for (int i = 0; i < numRows; ++i) {
        double* row = data[i];
        for (int j = 0; j < numCols; ++j) {
            if (j == 0) {
                // Write id as integer
                int id = static_cast<int>(row[j]);
                file << id;
            } else {
                file << "," << row[j];
            }
        }
        file << "\n"; // No extra comma at the end
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
    int numData = 0;
    Point** data = readCSV(datafile, numData);
    cout << "Loaded " << numData << " records from " << datafile << endl;

    // Load query
    int numQuery = 0;
    int* query = readQuery(queryfile, numQuery);
    cout << "Loaded " << numQuery << " queries from " << queryfile << endl;

    // Main MaxRank routine
    double** res = nullptr;
    int numRes = 0;

    double** cells = nullptr;
    int numCellsRes = 0;
    int cellEntrySize = 0;

    if (data[0]->dims > 2) {
        for (int q_idx = 0; q_idx < numQuery; ++q_idx) {
            int q = query[q_idx];
            // Reset global variables if any
            // Initialize or reset them here as needed

            cout << "#  Processing data point " << q << "  #" << endl;
            int idx = q - 1;  // Assuming query contains 1-based indices

            cout << "#  ";
            for (int d = 0; d < data[idx]->dims; ++d) {
                cout << data[idx]->coord[d] << " ";
            }
            cout << "  #" << endl;

            int maxrank = 0;
            Cell** mincells = nullptr;
            int numMinCells = 0;
            if (method == "BA") {
                // tie(maxrank, mincells) = ba_hd(data, numData, *data[idx]);
            } else {
                std::pair<int, Cell**> result = aa_hd(data, numData, *data[idx], numMinCells);
                maxrank = result.first;
                mincells = result.second;
                // Assume numMinCells is set appropriately within aa_hd
            }
            cout << "#  MaxRank: " << maxrank << "  NOfMincells: " << numMinCells << "  #" << endl;

            numRes++;
            res = (double**)realloc(res, numRes * sizeof(double*));
            res[numRes - 1] = (double*)malloc(2 * sizeof(double));
            res[numRes - 1][0] = (double)q;
            res[numRes - 1][1] = (double)maxrank;

            // Assuming mincells[0] is valid
            if (numMinCells > 0) {
                int dims = mincells[0]->feasible_pnt.dims;
                cellEntrySize = dims + 2; // id + coordinates + 1 (per 1 - sumCoords)

                numCellsRes++;
                cells = (double**)realloc(cells, numCellsRes * sizeof(double*));
                cells[numCellsRes - 1] = (double*)malloc(cellEntrySize * sizeof(double));

                cells[numCellsRes - 1][0] = (double)q;  // Inserisci l'ID (query)

                // Copia le coordinate e somma
                double sumCoords = 0.0;
                for (int d = 0; d < dims; ++d) {
                    double coord = mincells[0]->feasible_pnt.coord[d];
                    cells[numCellsRes - 1][d + 1] = coord;
                    sumCoords += coord;
                }



                // Calcola 1.0 - sumCoords e memorizza l'ultimo valore
                cells[numCellsRes - 1][dims + 1] = 1.0 - sumCoords;
            }


            // Clean up mincells
            if (mincells)
            {
                for (int c = 0; c < numMinCells; ++c) {
                    if (mincells[c])
                    {
                        free(mincells[c]);
                    }
                }
                free(mincells);
            }
        }
    } else {
        // Handle 2D case if needed
    }

    // Write results to CSV
    const std::string resHeaders[] = {"id", "maxrank"};
    writeCSV(R"(C:\Users\leona\Desktop\maxrank.csv)", res, numRes, 2, resHeaders, 2);

    const std::string cellHeaders[] = {"id", "query_found"};
    writeCSV(R"(C:\Users\leona\Desktop\cells.csv)", cells, numCellsRes, cellEntrySize, cellHeaders, 2);

    // Clean up
    for (int i = 0; i < numData; ++i) {
        delete data[i];
    }
    free(data);

    free(query);

    for (int i = 0; i < numRes; ++i) {
        free(res[i]);
    }
    free(res);

    for (int i = 0; i < numCellsRes; ++i) {
        free(cells[i]);
    }
    free(cells);

    return 0;
}
