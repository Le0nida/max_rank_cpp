//
// Created by leona on 06/08/2024.
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <Eigen/Dense>
#include "geom.h"
#include "maxrank.h"
#include "qtree.h"
#include "query.h"
#include "cell.h"

using namespace std;

// Funzione aggiornata per leggere il CSV e restituire un vector di Point
std::vector<std::shared_ptr<Point>> readCSV(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<std::shared_ptr<Point>> data;
    std::string line, word;

    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    // Skip the header line
    if (!getline(file, line)) {
        throw std::runtime_error("Empty or invalid file: " + filename);
    }

    while (getline(file, line)) {
        std::stringstream ss(line);
        std::vector<double> coord;
        int id = -1;
        bool first = true;

        // Leggi i dati separati da virgola
        std::string token;
        while (getline(ss, token, ',')) {
            if (first) {
                id = std::stoi(token);
                first = false;
            } else {
                coord.push_back(std::stod(token));
            }
        }

        if (id != -1 && !coord.empty()) {
            data.push_back(std::make_shared<Point>(coord, id));
        }
    }

    return data;
}

// Funzione aggiornata per leggere il file query
std::vector<int> readQuery(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<int> query;
    std::string line;

    while (getline(file, line)) {
        query.push_back(std::stoi(line));
    }

    return query;
}

// Funzione per scrivere i risultati su CSV
void writeCSV(const std::string& filename, const std::vector<std::vector<double>>& data, const std::vector<std::string>& headers) {
    std::ofstream file(filename);
    file << std::fixed << std::setprecision(15);  // Precisione per i valori floating point

    // Scrivi gli header
    for (size_t i = 0; i < headers.size(); ++i) {
        file << headers[i];
        if (i < headers.size() - 1) {
            file << ",";
        }
    }
    file << "\n";

    // Scrivi i dati
    for (const auto& row : data) {
        for (size_t i = 0; i < row.size(); ++i) {
            if (i == 0) {
                file << static_cast<int>(row[i]);  // Prima colonna è l'ID, lo scriviamo come intero
            } else {
                file << "," << row[i];
            }
        }
        file << "\n";
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

    // Carica i dati
    std::vector<std::shared_ptr<Point>> data = readCSV(datafile);
    cout << "Loaded " << data.size() << " records from " << datafile << endl;

    // Carica le query
    std::vector<int> query = readQuery(queryfile);
    cout << "Loaded " << query.size() << " queries from " << queryfile << endl;

    // Risultati
    std::vector<std::vector<double>> res;
    std::vector<std::vector<double>> cells;

    if (data[0]->dims > 2) {
        for (int q : query) {
            // Reset variabili globali se necessario

            cout << "#  Processing data point " << q << "  #" << endl;
            int idx = q - 1;  // Indice 1-based per le query

            cout << "#  ";
            for (int d = 0; d < data[idx]->dims; ++d) {
                cout << data[idx]->coord[d] << " ";
            }
            cout << "  #" << endl;

            int maxrank = 0;
            std::vector<std::shared_ptr<Cell>> mincells;
            if (method == "BA") {
                // tie(maxrank, mincells) = ba_hd(data, *data[idx]);
            } else {
                auto result = aa_hd(data, *data[idx]);
                maxrank = result.first;
                mincells = std::move(result.second);
            }
            cout << "#  MaxRank: " << maxrank << "  NOfMincells: " << mincells.size() << "  #" << endl;

            // Salva i risultati per la query attuale
            res.push_back({static_cast<double>(q), static_cast<double>(maxrank)});

            // Salva le celle
            if (!mincells.empty()) {
                int dims = mincells[0]->feasible_pnt.dims;
                std::vector<double> cellEntry(dims + 2); // id + coords + 1 - sum(coords)

                cellEntry[0] = static_cast<double>(q);  // ID della query

                double sumCoords = 0.0;
                for (int d = 0; d < dims; ++d) {
                    double coord = mincells[0]->feasible_pnt.coord[d];
                    cellEntry[d + 1] = coord;
                    sumCoords += coord;
                }

                cellEntry[dims + 1] = 1.0 - sumCoords;  // Calcola 1.0 - somma delle coordinate

                cells.push_back(cellEntry);
            }
        }
    } else {
        // Gestisci il caso per le dimensioni <= 2 se necessario
    }

    // Scrivi i risultati su CSV
    std::vector<std::string> resHeaders = {"id", "maxrank"};
    writeCSV(R"(C:\Users\leona\Desktop\maxrank.csv)", res, resHeaders);

    std::vector<std::string> cellHeaders = {"id", "query_found"};
    writeCSV(R"(C:\Users\leona\Desktop\cells.csv)", cells, cellHeaders);

    return 0;
}
