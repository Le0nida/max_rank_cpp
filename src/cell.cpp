#include "cell.h"
#include "halfspace.h"
#include "qnode.h"
#include <cmath>
#include <algorithm>
#include <cstring>   // Per memcpy, strcpy, strlen
#include <bitset>
#include <cstdlib>   // Per malloc, free
#include <src/Highs.h>
#include <limits>

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

// Implementazione del costruttore della classe Cell
Cell::Cell(int order, const char* mask, HalfSpace** covered, int numCovered,
           HalfSpace** halfspaces, int numHalfspaces, double** leaf_mbr, int dims,
           const Point& feasible_pnt)
    : order(order), numCovered(numCovered), numHalfspaces(numHalfspaces), dims(dims), feasible_pnt(feasible_pnt)
{
    // Copia della maschera
    int mask_len = strlen(mask);
    this->mask = (char*)malloc((mask_len + 1) * sizeof(char));
    strcpy(this->mask, mask);

    // Copia degli halfspaces coperti
    this->covered = (HalfSpace**)malloc(numCovered * sizeof(HalfSpace*));
    memcpy(this->covered, covered, numCovered * sizeof(HalfSpace*));

    // Copia degli halfspaces del nodo
    this->halfspaces = (HalfSpace**)malloc(numHalfspaces * sizeof(HalfSpace*));
    memcpy(this->halfspaces, halfspaces, numHalfspaces * sizeof(HalfSpace*));

    // Copia del leaf_mbr
    this->leaf_mbr = (double**)malloc(dims * sizeof(double*));
    for (int i = 0; i < dims; ++i) {
        this->leaf_mbr[i] = (double*)malloc(2 * sizeof(double));
        this->leaf_mbr[i][0] = leaf_mbr[i][0];
        this->leaf_mbr[i][1] = leaf_mbr[i][1];
    }
}

// Implementazione del distruttore della classe Cell
Cell::~Cell() {
    if (mask) {
        free(mask);
        mask = nullptr;
    }
    if (covered) {
        free(covered);
        covered = nullptr;
    }
    if (halfspaces) {
        free(halfspaces);
        halfspaces = nullptr;
    }
    if (leaf_mbr) {
        for (int i = 0; i < dims; ++i) {
            free(leaf_mbr[i]);
        }
        free(leaf_mbr);
        leaf_mbr = nullptr;
    }
}

// Implementazione della funzione issingular
bool Cell::issingular() const {
    for (int i = 0; i < numCovered; ++i) {
        if (covered[i]->arr != SINGULAR) {
            return false;
        }
    }
    return true;
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

// Funzione per cercare i minimi cell tramite programmazione lineare
Cell** searchmincells_lp(const QNode& leaf, char** hamstrings, int numHamstrings, int& numCells) {
    numCells = 0;
    Cell** cells = nullptr;

    int dims = leaf.dims;
    HalfSpace** leaf_covered = leaf.covered;
    int numLeafCovered = leaf.numCovered;

    HalfSpace** halfspaces = leaf.halfspaces;
    int numHalfspaces = leaf.numHalfspaces;

    if (numHalfspaces == 0) {
        // Caso in cui non ci sono halfspaces nel nodo foglia
        double** mbr = leaf.mbr;
        double* feasible_coords = (double*)malloc(dims * sizeof(double));
        for (int i = 0; i < dims; ++i) {
            feasible_coords[i] = (mbr[i][0] + mbr[i][1]) / 2.0;
        }
        Point feasible_pnt(feasible_coords, dims);
        free(feasible_coords);

        numCells = 1;
        cells = (Cell**)malloc(sizeof(Cell*));
        cells[0] = new Cell(0, "", leaf_covered, numLeafCovered, nullptr, 0, leaf.mbr, dims, feasible_pnt);
        return cells;
    }

    // Configurazione della funzione obiettivo
    int num_vars = dims + 1;
    double* c = (double*)calloc(num_vars, sizeof(double));
    c[dims] = -1.0;

    // Configurazione dei bounds
    double* bounds = (double*)malloc(2 * num_vars * sizeof(double));
    for (int d = 0; d < dims; ++d) {
        bounds[2 * d] = leaf.mbr[d][0];
        bounds[2 * d + 1] = leaf.mbr[d][1];
    }
    bounds[2 * dims] = 0.0; // Limite inferiore per la variabile slack
    bounds[2 * dims + 1] = kHighsInf; // Limite superiore per la variabile slack

    // Definizione del numero di vincoli
    int num_constraints = numHalfspaces + 1;

    // Alloca memoria per A_ub e b_ub (ma non inizializzarli qui)
    double* A_ub = (double*)malloc(num_constraints * num_vars * sizeof(double));
    double* b_ub = (double*)malloc(num_constraints * sizeof(double));

    for (int h = 0; h < numHamstrings; ++h) {
        char* hamstr = hamstrings[h];

        // **Re-inizializza A_ub e b_ub per ogni hamstr**
        // Inizializzazione di A_ub e b_ub
        for (int i = 0; i < num_constraints; ++i) {
            for (int j = 0; j < num_vars; ++j) {
                A_ub[i * num_vars + j] = 0.0;
            }
            b_ub[i] = 0.0; // Assicurati di inizializzare b_ub
        }
        // Ultima riga di A_ub (somma delle variabili <= 1)
        for (int j = 0; j < dims; ++j) {
            A_ub[numHalfspaces * num_vars + j] = 1.0;
        }
        A_ub[numHalfspaces * num_vars + dims] = 0.0;
        b_ub[numHalfspaces] = 1.0;

        for (int b = 0; b < numHalfspaces; ++b) {
            HalfSpace* hs = halfspaces[b];
            if (hamstr[b] == '0') {
                for (int i = 0; i < dims; ++i) {
                    A_ub[b * num_vars + i] = -hs->coeff[i];
                }
                A_ub[b * num_vars + dims] = 1.0; // **Corretto il segno**
                b_ub[b] = -hs->known;
            } else {
                for (int i = 0; i < dims; ++i) {
                    A_ub[b * num_vars + i] = hs->coeff[i];
                }
                A_ub[b * num_vars + dims] = 1.0; // **Corretto il segno**
                b_ub[b] = hs->known;
            }
        }

        // Risolvi il problema di programmazione lineare
        LinprogResult* result = linprog_highs(c, A_ub, b_ub, bounds, num_vars, num_constraints);

        if (result->status == static_cast<int>(HighsModelStatus::kOptimal)) {
            // Soluzione trovata
            double* solution = result->x;
            double* feasible_coords = (double*)malloc(dims * sizeof(double));
            for (int i = 0; i < dims; ++i) {
                feasible_coords[i] = solution[i];
            }
            Point feasible_pnt(feasible_coords, dims);
            free(feasible_coords);

            // Crea un nuovo Cell e aggiungilo alla lista
            Cell* cell = new Cell(0, hamstr, leaf_covered, numLeafCovered, halfspaces, numHalfspaces, leaf.mbr, dims, feasible_pnt);

            numCells++;
            cells = (Cell**)realloc(cells, numCells * sizeof(Cell*));
            cells[numCells - 1] = cell;

            free_linprog_result(result);
            break; // Esci dal loop poiché hai trovato una soluzione
        }

        free_linprog_result(result);
    }

    // Dealloca la memoria
    free(c);
    free(bounds);
    free(A_ub);
    free(b_ub);

    return cells;
}
