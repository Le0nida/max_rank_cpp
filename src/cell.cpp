#include "cell.h"
#include "halfspace.h"
#include "qnode.h"
#include <cmath>
#include <algorithm>
#include <cstring>   // Per memcpy, strcpy, strlen
#include <bitset>
#include <cstdlib>   // Per malloc, free
#include <src/Highs.h>
#include <fstream>
#include <utility>

// Definizione di ZEROEXTENT per la soglia del volume
#define ZEROEXTENT 1e-15

// Constructor for Interval class remains unchanged
Interval::Interval(const HalfLine& halfline, const std::pair<double, double>& range, int coversleft)
    : halfline(halfline), range(range), coversleft(coversleft)
{
}

bool Interval::issingular() const
{
    return std::all_of(covered.begin(), covered.end(),
                       [](const HalfSpace& hl) { return hl.arr == Arrangement::SINGULAR; });
}

Cell::Cell(int order, const std::string& mask, const std::vector<std::shared_ptr<HalfSpace>>& covered,
           const std::vector<std::shared_ptr<HalfSpace>>& halfspaces, const std::vector<std::pair<double, double>>& leaf_mbr,
           int dims, const Point& feasible_pnt)
    : order(order), mask(mask), covered(covered), halfspaces(halfspaces), leaf_mbr(leaf_mbr), dims(dims),
      feasible_pnt(feasible_pnt) {}

bool Cell::issingular() const {
    return std::all_of(covered.begin(), covered.end(), [](const std::shared_ptr<HalfSpace>& hs) {
        return hs->arr == Arrangement::SINGULAR;
    });
}

// Funzione per liberare il risultato della programmazione lineare
void free_linprog_result(LinprogResult* result) {
    if (result) {
        if (result->x) {
            delete[] result->x;
        }
        delete result;
    }
}

// Implementazione della funzione linprog_highs
LinprogResult* linprog_highs(const double* c, const double* A_ub, const double* b_ub,
                             const double* bounds, int num_vars, int num_constraints)
{
    auto* result = new LinprogResult;

    /*
    std::cout << "\n\nCC\n" << std::endl;
    for (int k = 0; k < 4; k++)
    {
        std::cout << c[k] << std::endl;
    }

    std::cout << "\n\nbounds\n" << std::endl;
    for (int k = 0; k < 8; k++)
    {
        std::cout << bounds[k] << std::endl;

    }

    std::cout << "\n\nb_ub\n" << std::endl;
    for (int k = 0; k < 11; k++)
    {
        std::cout << b_ub[k] << std::endl;
    }

    std::cout << "\n\nA_ub\n" << std::endl;
    int j = 0;
    for (int k = 0; k < (11) * 4; k++)
    {
        std::cout << A_ub[k] << ",";
        j++;

        if (j == 4)
        {
            std::cout << std::endl;
            j = 0;
        }

    }*/

    Highs highs;

    // Disabilita l'output a console di HiGHS
    HighsOptions options;
    options.output_flag = false;
    highs.passOptions(options);

    // Definisci le dimensioni del problema
    const int num_col = num_vars;
    const int num_row = num_constraints;

    // Coefficienti della funzione obiettivo
    std::vector<double> col_cost(c, c + num_col);

    // Preparazione di A_ub in formato CSR
    std::vector<int> A_start(num_row + 1);
    std::vector<int> A_index;
    std::vector<double> A_value;

    int count = 0;
    for (int i = 0; i < num_row; ++i) {
        A_start[i] = count;
        for (int j = 0; j < num_col; ++j) {
            double value = A_ub[i * num_col + j];
            if (value != 0.0) {
                A_index.push_back(j);
                A_value.push_back(value);
                count++;
            }
        }
    }
    A_start[num_row] = count;

    // Valori del lato destro (upper bounds)
    std::vector<double> row_upper(b_ub, b_ub + num_row);
    // Valori del lato sinistro (lower bounds)
    std::vector<double> row_lower(num_row, -kHighsInf);

    // Limiti delle variabili
    std::vector<double> col_lower(num_col);
    std::vector<double> col_upper(num_col);

    for (int i = 0; i < num_col; ++i) {
        col_lower[i] = bounds[2 * i];
        col_upper[i] = bounds[2 * i + 1];
    }

    // Aggiungi colonne e righe a HiGHS
    highs.addCols(num_col, col_cost.data(), col_lower.data(), col_upper.data(),
                  0, nullptr, nullptr, nullptr);
    highs.addRows(num_row, row_lower.data(), row_upper.data(),
                  A_value.size(), A_start.data(), A_index.data(), A_value.data());

    // Esegui HiGHS
    highs.run();

    // Ottieni la soluzione
    HighsSolution solution = highs.getSolution();
    HighsModelStatus model_status = highs.getModelStatus();

    // Prepara il risultato
    result->x = new double[num_col];
    for (int i = 0; i < num_col; ++i) {
        result->x[i] = solution.col_value[i];
    }
    result->fun = highs.getObjectiveValue();
    result->status = static_cast<int>(model_status);
    strncpy(result->message, highs.modelStatusToString(model_status).c_str(), 128);

    return result;
}

// Funzione per generare stringhe di Hamming
char** genhammingstrings(int strlen, int weight, int& numStrings) {
    bool botup;
    if (weight > (strlen / 2)) {
        weight = strlen - weight;
        botup = false;
    } else {
        botup = true;
    }

    int capacity = 1024;
    int size = 0;
    int* decstr = (int*)malloc(capacity * sizeof(int));

    if (weight == 0) {
        decstr[size++] = 0;
    } else if (weight == 1) {
        for (int b = 0; b < strlen; ++b) {
            if (size >= capacity) {
                capacity *= 2;
                decstr = (int*)realloc(decstr, capacity * sizeof(int));
            }
            decstr[size++] = 1 << b;
        }
    } else {
        int halfmax = (1 << (strlen - 1)) - 1;
        int curr_weight = 2;

        // Aggiungi i bit di base
        for (int b = 1; b < strlen; ++b) {
            if (size >= capacity) {
                capacity *= 2;
                decstr = (int*)realloc(decstr, capacity * sizeof(int));
            }
            decstr[size++] = (1 << b) + 1;
        }

        int* bases = (int*)malloc(size * sizeof(int));
        memcpy(bases, decstr, size * sizeof(int));
        int base_size = size;

        // Rimuovi gli elementi > halfmax
        int valid_size = 0;
        for (int i = 0; i < base_size; ++i) {
            if (bases[i] <= halfmax) {
                bases[valid_size++] = bases[i];
            }
        }

        base_size = valid_size;

        while (true) {
            // Genera nuovi bit shiftati
            while (base_size > 0) {
                int shifts_capacity = 1024;
                int shifts_size = 0;
                int* shifts = (int*)malloc(shifts_capacity * sizeof(int));

                for (int i = 0; i < base_size; ++i) {
                    int shifted = bases[i] << 1;
                    if (shifted <= halfmax) {
                        if (shifts_size >= shifts_capacity) {
                            shifts_capacity *= 2;
                            shifts = (int*)realloc(shifts, shifts_capacity * sizeof(int));
                        }
                        shifts[shifts_size++] = shifted;
                    }
                }

                // Aggiungi i risultati
                for (int i = 0; i < shifts_size; ++i) {
                    if (size >= capacity) {
                        capacity *= 2;
                        decstr = (int*)realloc(decstr, capacity * sizeof(int));
                    }
                    decstr[size++] = shifts[i];
                }

                free(bases);
                bases = shifts;
                base_size = shifts_size;
            }

            if (curr_weight < weight) {
                int new_bases_capacity = 1024;
                int new_bases_size = 0;
                int* new_bases = (int*)malloc(new_bases_capacity * sizeof(int));

                for (int i = 0; i < size; ++i) {
                    int new_dec = (decstr[i] << 1) + 1;
                    if (new_dec <= (1 << strlen) - 1) {
                        if (new_bases_size >= new_bases_capacity) {
                            new_bases_capacity *= 2;
                            new_bases = (int*)realloc(new_bases, new_bases_capacity * sizeof(int));
                        }
                        new_bases[new_bases_size++] = new_dec;
                    }
                }

                for (int i = 0; i < new_bases_size; ++i) {
                    if (size >= capacity) {
                        capacity *= 2;
                        decstr = (int*)realloc(decstr, capacity * sizeof(int));
                    }
                    decstr[size++] = new_bases[i];
                }

                free(bases);
                bases = new_bases;
                base_size = new_bases_size;
                curr_weight++;
            } else {
                free(bases);
                break;
            }
        }
    }

    numStrings = size;
    char** hamming_strings = (char**)malloc(numStrings * sizeof(char*));

    for (int i = 0; i < numStrings; ++i) {
        hamming_strings[i] = (char*)malloc((strlen + 1) * sizeof(char));
        for (int j = 0; j < strlen; ++j) {
            hamming_strings[i][j] = ((decstr[i] >> (strlen - j - 1)) & 1) + '0';
        }
        hamming_strings[i][strlen] = '\0';
    }

    free(decstr);

    if (!botup) {
        int decmax = (1 << strlen) - 1;
        for (int i = 0; i < numStrings; ++i) {
            int dec_val = 0;
            for (int j = 0; j < strlen; ++j) {
                dec_val = (dec_val << 1) | (hamming_strings[i][j] - '0');
            }
            int inverted = decmax - dec_val;
            for (int j = 0; j < strlen; ++j) {
                hamming_strings[i][j] = ((inverted >> (strlen - j - 1)) & 1) + '0';
            }
        }
    }

    return hamming_strings;
}

// Funzione per generare il file dei halfspace basato sulla stringa di Hamming
bool generateHalfspaceFile(const char* fileName, const QNode& leaf, const char* hamstr, int dims) {
    FILE* fout = fopen(fileName, "w");
    if (fout == NULL) {
        std::cout << "Errore nell'apertura del file " << fileName << std::endl;
        return false;
    }

    // Generazione di un punto fattibile (ad esempio, il centro dell'MBR)
    std::vector<std::pair<double, double>> mbr = leaf.mbr;
    std::vector<double> feasible_coords(dims);
    for (int i = 0; i < dims; ++i) {
        feasible_coords[i] = (mbr[i].first + mbr[i].second) / 2.0;
    }

    // Scrittura dell'intestazione del file
    fprintf(fout, "%d 1\n", dims); // Dimensionalità
    for (int i = 0; i < dims; ++i)
        fprintf(fout, "%f ", feasible_coords[i]); // Punto fattibile
    fprintf(fout, "\n");
    fprintf(fout, "%d\n", dims + 1);

    // Numero totale di iperpiani: 2*dims (limiti dell'MBR) + numero di halfspace
    int numHyperplanes = 2 * dims + leaf.halfspaces.size();
    fprintf(fout, "%d\n", numHyperplanes);

    // Scrittura degli iperpiani limitanti (l'MBR del nodo foglia)
    for (int i = 0; i < dims; ++i) {
        // Limite inferiore: x_i >= mbr[i].first
        for (int j = 0; j < dims; ++j)
            fprintf(fout, "%d ", (i == j) ? 1 : 0);
        fprintf(fout, "%f\n", -mbr[i].first);

        // Limite superiore: x_i <= mbr[i].second
        for (int j = 0; j < dims; ++j)
            fprintf(fout, "%d ", (i == j) ? -1 : 0);
        fprintf(fout, "%f\n", mbr[i].second);
    }

    // Scrittura dei halfspace basati sulla stringa di Hamming
    int index = 0;
    for (const auto& hsPtr : leaf.halfspaces) {
        auto hs = hsPtr;
        if (hamstr[index] == '0') {
            // Utilizza il halfspace così com'è
            for (int i = 0; i < dims; ++i)
                fprintf(fout, "%f ", hs->coeff[i]);
            fprintf(fout, "%f\n", -hs->known);
        } else if (hamstr[index] == '1') {
            // Utilizza il halfspace opposto
            for (int i = 0; i < dims; ++i)
                fprintf(fout, "%f ", -hs->coeff[i]);
            fprintf(fout, "%f\n", hs->known);
        }
        index++;
    }

    fclose(fout);
    return true;
}

// Funzione per leggere il volume dal file di output generato da qconvex
double readVolumeFromFile(const char* volumeFileName) {
    std::ifstream inFile(volumeFileName);
    if (!inFile.is_open()) {
        std::cout << "Errore nell'apertura del file del volume " << volumeFileName << std::endl;
        return 0.0;
    }

    std::string line;
    double volume = 0.0;
    while (std::getline(inFile, line)) {
        size_t pos = line.find("volume:");
        if (pos != std::string::npos) {
            std::string volStr = line.substr(pos + 7);
            volume = atof(volStr.c_str());
            break;
        }
    }
    inFile.close();
    return volume;
}

// Funzione per cercare le min-cells utilizzando qhalf e qconvex
std::vector<std::shared_ptr<Cell>> searchmincells_lp(const QNode& leaf, char** hamstrings, int numHamstrings) {
    int dims = leaf.dims;
    std::vector<std::shared_ptr<Cell>> cells;

    // Debug: stampa le informazioni generali sul nodo foglia
    std::cout << "Inizio ricerca min-cells per nodo foglia con " << dims << " dimensioni e " << leaf.halfspaces.size() << " halfspaces." << std::endl;

    auto leaf_covered = leaf.getTotalCovered();
    size_t numHalfspaces = leaf.halfspaces.size();

    if (numHalfspaces == 0) {
        // Nessun halfspace da considerare, crea una cella che copre l'MBR del nodo
        std::vector<std::pair<double, double>> mbr = leaf.mbr;
        std::vector<double> feasible_coords(dims);
        for (int i = 0; i < dims; ++i) {
            feasible_coords[i] = (mbr[i].first + mbr[i].second) / 2.0;
        }
        Point feasible_pnt(feasible_coords, dims);
        std::cout << "Nessun halfspace, crea cella che copre l'MBR del nodo." << std::endl;
        cells.push_back(std::make_shared<Cell>(0, "", leaf_covered, std::vector<std::shared_ptr<HalfSpace>>{}, leaf.mbr, dims, feasible_pnt));
        return cells;
    }

    // Percorsi per i file temporanei
    char halfspaceFileName[256];
    char volumeFileName[256];

    // Per ogni stringa di Hamming
    for (int h = 0; h < numHamstrings; ++h) {
        char* hamstr = hamstrings[h];

        // Genera nomi di file univoci
        sprintf(halfspaceFileName, "halfspace_%d_%d.txt", rand(), h);
        sprintf(volumeFileName, "volume_%d_%d.txt", rand(), h);

        // Debug: stampa la stringa di Hamming corrente
        std::cout << "Processing Hamming string #" << h + 1 << ": " << hamstr << std::endl;

        // Genera il file dei halfspace basato sulla stringa di Hamming
        bool success = generateHalfspaceFile(halfspaceFileName, leaf, hamstr, dims);
        if (!success) {
            // Impossibile generare il file dei halfspace, passa alla prossima stringa di Hamming
            std::cout << "Errore nella generazione del file dei halfspaces per la stringa di Hamming #" << h + 1 << std::endl;
            continue;
        }

        // Debug: stampa il nome del file halfspace generato
        std::cout << "Generato file halfspace: " << halfspaceFileName << std::endl;

        // Esegue qhalf e qconvex per calcolare il volume
        char command[512];
        sprintf(command, R"(C:\Users\leona\Documents\Projects\max_rank_cpp\bin\qhalf.exe Fp < %s | C:\Users\leona\Documents\Projects\max_rank_cpp\bin\qconvex FA > %s)", halfspaceFileName, volumeFileName);
        int systemResult = system(command);

        // Debug: stampa il comando eseguito
        std::cout << "Eseguito comando: " << command << std::endl;

        if (systemResult != 0) {
            // Errore nell'esecuzione di qhalf/qconvex
            std::cout << "Errore nell'esecuzione di qhalf/qconvex per la stringa di Hamming #" << h + 1 << std::endl;
            // Rimuove i file temporanei
            remove(halfspaceFileName);
            remove(volumeFileName);
            continue;
        }

        // Legge il volume dal file di output
        double volume = readVolumeFromFile(volumeFileName);

        // Debug: stampa il volume calcolato
        std::cout << "Volume calcolato per la stringa di Hamming #" << h + 1 << ": " << volume << std::endl;

        if (volume > ZEROEXTENT) {
            // Trovata una min-cell valida
            std::cout << "Min-cell valida trovata per la stringa di Hamming #" << h + 1 << std::endl;

            // Crea un punto fattibile (ad esempio, il centro dell'MBR del nodo)
            std::vector<std::pair<double, double>> mbr = leaf.mbr;
            std::vector<double> feasible_coords(dims);
            for (int i = 0; i < dims; ++i) {
                feasible_coords[i] = (mbr[i].first + mbr[i].second) / 2.0;
            }
            Point feasible_pnt(feasible_coords, dims);

            // Crea una nuova Cell e aggiungila alla lista
            cells.push_back(std::make_shared<Cell>(0, hamstr, leaf_covered, leaf.halfspaces, leaf.mbr, dims, feasible_pnt));

            // Rimuove i file temporanei
            remove(halfspaceFileName);
            remove(volumeFileName);

            // Debug: conferma la min-cell trovata e il loop che si interrompe
            std::cout << "Min-cell trovata, interruzione del ciclo." << std::endl;
            break; // Esci dal loop poiché hai trovato una min-cell valida
        }

        // Debug: stampa se la min-cell trovata ha un volume nullo o è stata scartata
        std::cout << "Nessuna min-cell valida trovata per la stringa di Hamming #" << h + 1 << " (volume nullo)." << std::endl;

        // Rimuove i file temporanei
        remove(halfspaceFileName);
        remove(volumeFileName);
    }

    return cells;
}