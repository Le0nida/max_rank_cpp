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
#include <set>
#include <utility>
#include <fstream>
#include <random>

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

// Funzione per cercare i minimi cell tramite programmazione lineare
std::vector<std::shared_ptr<Cell>> searchmincells_lp(const QNode& leaf, char** hamstrings, int numHamstrings) {
    int dims = leaf.dims;
    std::vector<std::shared_ptr<Cell>> cells;

    auto leaf_covered = leaf.getTotalCovered();

    size_t numHalfspaces = leaf.halfspaces.size();
    if (numHalfspaces == 0) {
        std::vector<std::pair<double, double>> mbr = leaf.mbr;
        std::vector<double> feasible_coords(dims);
        for (int i = 0; i < dims; ++i) {
            feasible_coords[i] = (mbr[i].first + mbr[i].second) / 2.0;
        }
        Point feasible_pnt(feasible_coords, dims);
        cells.push_back(std::make_shared<Cell>(0, "", leaf_covered, std::vector<std::shared_ptr<HalfSpace>>{}, leaf.mbr, dims, feasible_pnt));
        return cells;
    }

    // Configurazione della funzione obiettivo
    int num_vars = dims + 1;
    double* c = (double*)calloc(num_vars, sizeof(double));
    c[dims] = -1.0;

    // Configurazione dei bounds
    double* bounds = (double*)malloc(2 * num_vars * sizeof(double));
    for (int d = 0; d < dims; ++d) {
        bounds[2 * d] = leaf.mbr[d].first;
        bounds[2 * d + 1] = leaf.mbr[d].second;
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
            auto hs = leaf.halfspaces[b];
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
            std::vector<double> feasible_coords;
            feasible_coords.reserve(dims);
            for (int i = 0; i < dims; ++i) {
                feasible_coords[i] = solution[i];
            }
            Point feasible_pnt(feasible_coords, dims);

            // Crea un nuovo Cell e aggiungilo alla lista
            cells.push_back(std::make_shared<Cell>(0, hamstr, leaf_covered, leaf.halfspaces, leaf.mbr, dims, feasible_pnt));

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


int Dimen = 4 - 1;
#define MAXLOOP 1000
#define MAXDIM 10
#define ZEROEXTENT 1e-15
#define MAXNOBINSTRINGTOCHECK 60000

void myitoa(unsigned long val, char *buf, unsigned radix) {
    char *p;
    char *firstdig;
    char temp;
    unsigned digval;

    p = buf;
    firstdig = p;

    do {
        digval = (unsigned) (val % radix);
        val /= radix;

        if (digval > 9)
            *p++ = (char) (digval - 10 + 'a');
        else
            *p++ = (char) (digval + '0');

    } while (val > 0);

    *p-- = '\0';
    do {
        temp = *p;
        *p = *firstdig;
        *firstdig = temp;
        --p;
        ++firstdig;
    } while (firstdig < p);
}

void
GenString(long int stringLen, long int HammingDistance, long int currentLen, long int start, vector<char> &hammingStr,
          std::multimap<int, vector<char> > &binString) {
    if (currentLen < 0) return;

    for (long int i = start; i < stringLen; i++) {
        for (long int j = start; j < i; j++)
            hammingStr.push_back('0');
        hammingStr.push_back('1');
        GenString(stringLen, HammingDistance, currentLen - 1, i + 1, hammingStr, binString);
        if (currentLen == 0) {
            for (long int j = i + 1; j < stringLen; j++)
                hammingStr.push_back('0');

            vector<char> tmpHammingStr = hammingStr;
            typedef std::multimap<int, vector<char> >::value_type VT;
            binString.insert(VT(HammingDistance, tmpHammingStr));

            //cout << "Generated Hamming string: " << string(tmpHammingStr.begin(),tmpHammingStr.end()) << endl;

            for (long int j = i + 1; j < stringLen; j++)
                hammingStr.pop_back();
        }
        hammingStr.pop_back();
        for (long int j = start; j < i; j++)
            hammingStr.pop_back();
    }
}

void GenLenNBinaryString(long int len1, long int HammingDistance, std::multimap<int, vector<char> > &binString) {
    vector<char> hammingStr;

    if (HammingDistance == 0) {
        for (long int i = 0; i < len1; i++) hammingStr.push_back('0');
        typedef std::multimap<int, vector<char> >::value_type VT;
        binString.insert(VT(0, hammingStr));
        return;
    }

    if (HammingDistance == 1) {
        for (long int i = 0; i < len1; i++) {
            hammingStr.clear();
            for (long int j = 0; j < i; j++)
                hammingStr.push_back('0');
            hammingStr.push_back('1');
            for (long int j = i + 1; j < len1; j++)
                hammingStr.push_back('0');

            typedef std::multimap<int, vector<char> >::value_type VT;
            binString.insert(VT(1, hammingStr));
        }
        return;
    }


    GenString(len1, HammingDistance, HammingDistance - 1, 0, hammingStr, binString);

    return;
}


void GenBinaryString(long int len1, long int Max, std::multimap<int, vector<char> > &binString) {
    //ofstream fOut;
    //fOut.open("CombSubspace.txt", ios::app);

    typedef std::multimap<int, vector<char> >::value_type VT;

    binString.clear();

    char str[128], *b;
    //int len1=int(log(Max)/log(2));
    for (int i = 0; i < Max; i++) {
        myitoa(i, str, 2);
        int len = len1 - strlen(str);
        b = new char[Max + 1];
        strcpy(b, "");
        for (int j = 1; j <= len; j++) strcat(b, "0");
        strcat(b, str);
        std::cout << "i= " << i << ", binary string: " << b << std::endl;
        //fOut << b << endl;

        //count the number of '1's
        vector<char> tmpVec;
        int NoOfOnes = 0;
        for (int j = 0; j < strlen(b); j++) {
            if (b[j] == '1') NoOfOnes++;
            tmpVec.push_back(b[j]);
        }
        binString.insert(VT(NoOfOnes, tmpVec));

        delete b;
    }
    //fOut.close();

    return;
}

bool Cell::testHalfspacePair(std::shared_ptr<HalfSpace> HS1, long int IdxHS1, std::shared_ptr<HalfSpace> HS2, long int IdxHS2, const std::vector<std::pair<double, double>>& subDataSpace,
                                 std::multimap<int, string> &InValidHammingStr) //test whether two halfspaces are compatible w.r.t Hamming distance 00,01,10,11
{

    std::map<long int, long int>::iterator IntInt_mItr;
    typedef std::map<long int, long int>::value_type IntIntVT;

    typedef std::multimap<int, string>::value_type msVT;

    char HammingStr[4][3] = {"00", "01", "10", "11"};

    bool interiorPtExists;
    float InteriorPt[MAXDIM];

    //randomization process to generate an interior point
    long int loops = 0;
    long int count = 0;
    int HammingDistance;

    for (int i = 0; i < 4; i++) {
        loops = 0;
        interiorPtExists = true;
        while (true) {
            loops++;
            if (loops >= MAXLOOP + 200) //low probability that an interior point exists for current halfspaces
            {
                //cout << "exceed the maximal loop limits, none interior point exists!" << endl;
                interiorPtExists = false;
                break;
            }
            for (int j = 0; j < Dimen; j++) {
                double lower = subDataSpace[j].first;
                double upper = subDataSpace[j].second;
                InteriorPt[j] = lower + (upper - lower) * (double(rand()) / RAND_MAX);
            }

            count = 0;
            if (strcmp(HammingStr[i], "00") == 0)    //for case: ax_1+bx_2+... > d, lx_1+mx_2+... > t
            {
                HammingDistance = 0;
                float sum = 0;
                for (int j = 0; j < Dimen; j++)
                    sum = sum + HS1.get()->coeff[j] * InteriorPt[j];
                if (sum > HS1.get()->coeff[Dimen]) count++;
                sum = 0;
                for (int j = 0; j < Dimen; j++)
                    sum = sum + HS2.get()->coeff[j] * InteriorPt[j];
                if (sum > HS2.get()->coeff[Dimen]) count++;
            } else if (strcmp(HammingStr[i], "01") == 0)   //for case: ax_1+bx_2+... > d, lx_1+mx_2+... <= t
            {
                HammingDistance = 1;
                float sum = 0;
                for (int j = 0; j < Dimen; j++)
                    sum = sum + HS1.get()->coeff[j] * InteriorPt[j];
                if (sum > HS1.get()->coeff[Dimen]) count++;
                sum = 0;
                for (int j = 0; j < Dimen; j++)
                    sum = sum + HS2.get()->coeff[j] * InteriorPt[j];
                if (sum <= HS2.get()->coeff[Dimen]) count++;
            } else if (strcmp(HammingStr[i], "10") == 0)   //for case: ax_1+bx_2+... <= d, lx_1+mx_2+... > t
            {
                HammingDistance = 1;
                float sum = 0;
                for (int j = 0; j < Dimen; j++)
                    sum = sum + HS1.get()->coeff[j] * InteriorPt[j];
                if (sum <= HS1.get()->coeff[Dimen]) count++;
                sum = 0;
                for (int j = 0; j < Dimen; j++)
                    sum = sum + HS2.get()->coeff[j] * InteriorPt[j];
                if (sum > HS2.get()->coeff[Dimen]) count++;
            } else if (strcmp(HammingStr[i], "11") == 0)   //for case: ax_1+bx_2+... <= d, lx_1+mx_2+... <= t
            {
                HammingDistance = 2;
                float sum = 0;
                for (int j = 0; j < Dimen; j++)
                    sum = sum + HS1.get()->coeff[j] * InteriorPt[j];
                if (sum <= HS1.get()->coeff[Dimen]) count++;
                sum = 0;
                for (int j = 0; j < Dimen; j++)
                    sum = sum + HS2.get()->coeff[j] * InteriorPt[j];
                if (sum <= HS2.get()->coeff[Dimen]) count++;
            }
            if (count == 2) {
                //cout << "We have found an interior point!" << endl;
                break;
            }
        }
        if (!interiorPtExists) {
            //long int Val1=HS1,Val2=HS2;        //form incompatible pairs by IDs of the halfspaces
            long int Val1 = IdxHS1, Val2 = IdxHS2;    //form incompatible pairs by indices of the halfspaces in set 'intersectedHalfSpace'

            if (Val1 > Val2) {
                long int tmpInt;
                tmpInt = Val1;
                Val1 = Val2;
                Val2 = tmpInt;
            }
            char m_Buf[1024], m_Buf1[256];
            myitoa(Val1, m_Buf1, 10);
            strcpy(m_Buf, m_Buf1);
            strcat(m_Buf, "|");
            myitoa(Val2, m_Buf1, 10);
            strcat(m_Buf, m_Buf1);
            strcat(m_Buf, "|");
            strcat(m_Buf, HammingStr[i]);
            string tmpString = m_Buf;

            InValidHammingStr.insert(msVT(HammingDistance, tmpString));

            //cout << "Hamming Distance: " << HammingDistance << ", String=" << tmpString << endl;

        }
    }
    return true;
}


bool Cell::GenHammingHalfSpaces(char *OutFileName, const int Dimen, vector<char> &HammingString,
                                    std::vector<std::shared_ptr<HalfSpace>> HalfSpaceIDs, const std::vector<std::pair<double, double>>& subDataSpace) {
    int NoOfHyperplanes = 0;
    long int NoOfHalfSpaces = HalfSpaceIDs.size();

    bool interiorPtExists = true;
    float InteriorPt[MAXDIM];

    NoOfHyperplanes = 2 * Dimen + NoOfHalfSpaces;   //total number of halfspaces

    //randomization process to generate an interior point
    long int loops = 0;

    //for (int i=0;i<Dimen;i++)
    //     cout << "[" << subDataSpace[i] << "," << subDataSpace[Dimen+i] << "]"<< endl;

    while (true) {
        loops++;
        //cout << "loop " << loops << endl;
        if (loops >=
            MAXLOOP) //there is a very low probability that an interior point exists for current set of halfspaces
        {
            //cout << "exceed the maximal loop limits, none interior point exists!" << endl;
            interiorPtExists = false;
            break;
        }

        // The original was using rand(), which has a bad randomness
        // Also directly generate numbers between 0 and 1
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0, 1);

        //Discard query points over y = 1 - x (in 3D), as the last parameter (yield by 1 - x - y) would be negative
        //By doing so eventual mincells in such regions will not be found
        while(true){
            for (int i = 0; i < Dimen; i++) {
                double lower = subDataSpace[i].first;
                double upper = subDataSpace[i].second;
                InteriorPt[i] = lower + (upper - lower) * dis(gen);
            }

            float sum = 0;
            for (int i = 0; i < Dimen - 1; i++)
                sum += InteriorPt[i];

            if (InteriorPt[Dimen - 1] < 1 - sum)
                break;
        }

        /*
        for (int i = 0; i < Dimen; i++)
            InteriorPt[i] = subDataSpace[i] + (subDataSpace[Dimen + i] - subDataSpace[i]) * (float(rand()) / RAND_MAX);
        */

        int index = 0;
        long int count = 0;
        for (auto sItr = HalfSpaceIDs.begin(); sItr != HalfSpaceIDs.end(); sItr++) {
            if (HammingString[index] == '0')    //for the case where ax_1+bx_2+... <= d    1
            {
                auto hsID = *sItr;
                float sum = 0;
                for (int i = 0; i < Dimen; i++)
                    sum = sum + hsID.get()->coeff[i] * InteriorPt[i];
                if (sum <= hsID.get()->coeff[Dimen]) count++;
            } else if (HammingString[index] == '1')   //for the case where ax_1+bx_2+... > d   0
            {
                auto hsID = *sItr;
                float sum = 0;
                for (int i = 0; i < Dimen; i++)
                    sum = sum + hsID.get()->coeff[i] * InteriorPt[i];
                if (sum >= hsID.get()->coeff[Dimen]) count++;
            }
            index++;
        }
        if (count == HalfSpaceIDs.size()) {
            //cout << "We have found an interior point!" << endl;
            //for (int i = 0; i < Dimen; i++) cout << InteriorPt[i] << " ";
            //cout << endl;
            break;
        }
    }
    if (!interiorPtExists) {
        //cout << "Oops, it's unlikely that an interior point exists!!" << endl;
        return false;
    }
    //

    FILE *fout1 = fopen(OutFileName, "w");
    if (fout1 == NULL) {
        std::cout << "Error in opening file " << OutFileName << std::endl;
        getchar();
        exit(0);
    }

    fprintf(fout1, "%d 1\n", Dimen);   //the dimensionality
    for (int i = 0; i < Dimen; i++)
        fprintf(fout1, "%f ",
                InteriorPt[i]); //the feasible point found by Monte Carlo process above   *****************
    fprintf(fout1, "\n");
    fprintf(fout1, "%d\n", Dimen + 1);   //dimensionality + 1
    fprintf(fout1, "%d\n", NoOfHyperplanes);  //the total number of hyperplanes for intersection

    //output the 2*dimen bounding facets of the hypercube of the sub-dataspace [x_min,y_min,z_min,...] [x_max,y_max,z_max,...]
    int Cooef = -1;
    for (int i = 1; i <= 2; i++) {
        if (i == 2)
            Cooef = -Cooef;
        for (int j = 0; j < Dimen; j++) {
            for (int m = 0; m < Dimen; m++)
                fprintf(fout1, "%d  ", (m == j) ? Cooef : 0);
            double bound = (i == 2) ? -subDataSpace[j].second : subDataSpace[j].first;
            fprintf(fout1, "%f\n", bound);
        }
    }
    ///

    int index = 0;
    for (auto sItr = HalfSpaceIDs.begin(); sItr != HalfSpaceIDs.end(); sItr++) {
        if (HammingString[index] == '0')    //for the case where ax_1+bx_2+... <= d  1
        {
            auto hsID = *sItr;
            for (int i = 0; i < Dimen; i++)
                fprintf(fout1, "%f ", hsID.get()->coeff[i]);
            fprintf(fout1, "%f\n", -hsID.get()->coeff[Dimen]);   //the offset
        } else if (HammingString[index] == '1')   //for the case where ax_1+bx_2+... > d   0
        {
            auto hsID = *sItr;
            for (int i = 0; i < Dimen; i++)
                fprintf(fout1, "%f ", -hsID.get()->coeff[i]);
            fprintf(fout1, "%f\n", hsID.get()->coeff[Dimen]);   //the offset
        }
        index++;
    }
    fclose(fout1);
    //cout << "output data for halfspace intersection finished!" << endl;

    return true;
}

bool MbrIsValid(const int &Dimen, const float hs[], const std::vector<std::pair<double, double>>& mbr,
                vector<string> &Comb) {   //position of an MBR to a halfspace: is the point above, below, or intersected by the halfspace?

    int numAbove = 0;
    int numBelow = 0;
    int numOn = 0;

    long int numOfVertices = 0;
    numOfVertices = Comb.size();

    for (int i = 0; i < numOfVertices; i++) {
        std::vector<double> coord;

        long int numOfDimen = Comb[i].size();
        for (int j = 0; j < Dimen; j++) {
            if (Comb[i][j] == '0')
                coord[j] = mbr[j].first;
            else if (Comb[i][j] == '1')
                coord[j] = mbr[j].second;
        }
        float sum = 0;
        for (int k = 0; k < numOfDimen; k++) sum = sum + coord[k];
        if (sum > hs[Dimen]) numAbove++;
        if (sum < hs[Dimen]) numBelow++;
    }

    if (numAbove == numOfVertices) return false;
    if (numBelow == numOfVertices) return true;

}

std::vector<string> readCombinations() {

    std::vector<string> Comb;
    FILE *fp;
    char *token;
    char m_separator[] = " \n\t";
    char buf[512];
    string FileName[] = {"bin/Comb2D.txt", "bin/Comb3D.txt", "bin/Comb4D.txt", "bin/Comb5D.txt", "bin/Comb6D.txt", "bin/Comb7D.txt",
                         "bin/Comb8D.txt", "bin/Comb9D.txt"};
    string strComb;

    fp = fopen(FileName[Dimen - 2].c_str(), "r");
    if (fp == NULL) {
        std::cout << "error in fileopen!" << std::endl;
        exit(0);
    }

    fgets(buf, 512, fp);
    if (atoi(buf) != Dimen) {
        std::cout << "Error! Dimensions are not equal!" << std::endl;
        exit(0);
    }
    while (fgets(buf, 512, fp) != NULL) {
        token = strtok(buf, m_separator);
        //while (token != NULL){
        //	token = strtok(NULL,m_separator);
        //}
        string strComb = token;
        Comb.push_back(strComb);
    }
    fclose(fp);

    if (false)
        std::cout << "QuadTree Dimen=" << Dimen << ", reading combination finished!" << std::endl;

    return Comb;
}

long int Cell::optimizedInNodeIntersection(vector<std::pair<long, QNode *> > &Leaves,
                                               vector<std::set<std::shared_ptr<HalfSpace>> > &minCellHalfSpaces,
                                               vector<vector<char> > &binaryString) //optimization of within node intersection
{

    bool verbose = false;

    std::multimap<long int, QNode *> nodesToIntersect;   //store the nodes (in ascending order) by using the number of their intersected halfspaces
    //multimap<long int, QuadNode *>::iterator itr, itr1;
    vector<std::pair<long, QNode *> >::iterator itr, itr1;

    std::multimap<int, string>::iterator msItr;

    std::map<long int, long int>::iterator IntInt_mItr;

    std::set<string>::iterator ssItr, ssItr1;
    std::set<long int>::iterator sItr, sItr1;

    vector<string> FilesToRemove;

    if (Leaves.empty()) {
        std::cout << "There is no leaf nodes to perform intersection!" << std::endl;
        return -1;
    }

    FILE *fp_tmpIn;

    char Buf[1024];
    char *token;
    char m_seperator[] = " :\n\t";

    char volumeFilename[2048] = "Vol";
    myitoa(Dimen, Buf, 10);
    strcat(volumeFilename, Buf);
    myitoa(rand(), Buf, 10);
    strcat(volumeFilename, Buf);
    strcat(volumeFilename, "D.txt");

    char namePrefix[] = "./tmp/HalfSpaces";
    char nameSuffix[] = ".txt";
    char halfspaceFileName[1024];

    long int minOrder = INT_MAX;

    long int NoOfCoveredHS, NoOfIntersectedHS;
    bool NotFoundAllMinCells = true;

    long int NoOfInvalidLeaves = 0;
    long int NoOfNodesLeft = Leaves.size();

    long int NoOfBitStringsProcessed = 0;
    long int NoOfHalfSpacesInNode;
    long int NoOfPrunedBitStrings = 0;
    long int NoOfZeroExtentBinStrings = 0;
    long int NoOfDiscardedCells = 0;

    std::vector<string> Comb = readCombinations();
    for (itr = Leaves.begin(); itr != Leaves.end();)    // && NotFoundAllMinCells;)
    {

        //prune away leaf nodes that lie about hyperplane q_1+q2+...+q_d < 1;
        float queryPlane[MAXDIM];
        for (int i = 0; i < Dimen + 1; i++) queryPlane[i] = 1;
        bool isValid = MbrIsValid(Dimen, queryPlane, (*itr).second->mbr, Comb);
        if (!isValid) {
            //cout << "Leaf node " << (*itr).second->NodeID << " is pruned!" << endl;
            NoOfInvalidLeaves++;
            ++itr;
            continue;
        }

        NoOfCoveredHS = (*itr).first;
        if (NoOfCoveredHS > minOrder)
            break; //terminate searching min-cells, cause no cell with smaller order exists any more

        NoOfIntersectedHS = 0;
        NoOfIntersectedHS = ((*itr).second)->halfspaces.size();
        if (NoOfIntersectedHS > 0)
            nodesToIntersect.insert(std::pair<long int, QNode *>(NoOfIntersectedHS, (*itr).second));
        else {
            itr++;
            continue;
        }
        itr1 = itr;
        itr1++;
        nodesToIntersect.clear();
        while (true) {
            if (itr1 == Leaves.end() || NoOfCoveredHS != (*itr1).first) break;

            NoOfIntersectedHS = 0;
            NoOfIntersectedHS = ((*itr1).second)->halfspaces.size();
            if (NoOfIntersectedHS > 0)
                nodesToIntersect.insert(std::pair<long int, QNode *>(NoOfIntersectedHS, (*itr1).second));
            itr1++;
        }
        if (itr1 == Leaves.end()) break;
        itr = itr1;

        //cout << "Processing node " << itr1->second->NodeID << endl;

        //perform in-node halfspace intersection for nodes sorting in ascending order according to the number of intersected halfspaces
        //cout << "Size of nodesToIntersect: " << nodesToIntersect.size() << endl;
        for (std::multimap<long, QNode *>::iterator itr1 = nodesToIntersect.begin();
             itr1 != nodesToIntersect.end(); itr1++) {


            NoOfHalfSpacesInNode = itr1->first;

            NoOfNodesLeft--;
            if (NoOfHalfSpacesInNode <= 1) continue;
            //cout << "Number of nodes left to process :" << NoOfNodesLeft << endl;
            //cout <<"examining node " << itr1->second->NodeID << ", #inter.HS:" << itr1->second->intersectedHalfspace->size() << endl;

            //test compatibility of Hamming distance for each pair of halfspaces
            std::multimap<int, string> InValidHammingStr;
            long int idx1 = 0, idx2;
            for (auto sItr = (itr1->second->halfspaces).begin();
                 sItr != (itr1->second->halfspaces).end(); sItr++) {
                auto hs1 = (*sItr);
                vector<std::shared_ptr<HalfSpace>>::iterator sItr1;
                sItr1 = sItr;
                sItr1++;
                idx2 = idx1 + 1;
                for (; sItr1 != (itr1->second->halfspaces).end(); sItr1++) {
                    auto hs2 = (*sItr1);
                    testHalfspacePair(hs1, idx1, hs2, idx2, itr1->second->mbr, InValidHammingStr);
                    idx2++;
                }
                idx1++;
            }
            /*    test output of the incompatible Halfspace pair
            size_t substrPos1,substrPos2;
            for (msItr=InValidHammingStr.begin();msItr!=InValidHammingStr.end();msItr++)
            {
                 cout << "Hamming Distance: " << msItr->first << ", String=" << msItr->second << endl;
                 substrPos1=msItr->second.find("|");
                 cout << msItr->second.substr(0,substrPos1) << " ";
                 substrPos2=msItr->second.find("|",substrPos1+1);
                 cout << msItr->second.substr(substrPos1+1,substrPos2-substrPos1-1) << " ";
                 cout << msItr->second.substr(substrPos2+1) << endl;
            }
            //*/
            //end of testing compatibility

            //intersect the halfspaces in each leaf node
            long int NoOfCombinations = (int) pow(2.0, NoOfHalfSpacesInNode);
            std::multimap<int, vector<char> > binString;
            long int HammingDistance = 0;
            long int LoopCounter = 0;
            bool stopIncrHammingDist = false;
            long int HSidx1, HSidx2;
            size_t substrPos1, substrPos2;
            while ((HammingDistance <= NoOfHalfSpacesInNode) && !stopIncrHammingDist) {
                //if (LoopCounter >= MAXNOBINSTRINGTOCHECK) break;
                if (minOrder < INT_MAX && (HammingDistance + NoOfCoveredHS) > minOrder)
                    break;

                binString.clear();
                GenLenNBinaryString(NoOfHalfSpacesInNode, HammingDistance, binString);  //generate all the combinations
                std::multimap<int, vector<char> >::iterator itrHamming;
                for (itrHamming = binString.begin(); itrHamming != binString.end(); itrHamming++) {
                    LoopCounter++;
                    if (LoopCounter >= MAXNOBINSTRINGTOCHECK) {
                        std::cout << "Maximum loop limit exceeds!" << std::endl;
                        stopIncrHammingDist = true;
                        break;
                    }

                    NoOfBitStringsProcessed++;
                    if (verbose)
                        std::cout << "Testing Hamming string: "
                             << string(itrHamming->second.begin(), itrHamming->second.end()) << std::endl;
                    //Optimzed part Task 4: prune away Hamming String that contain incompatible halfspace pairs
                    bool isValid = true;
                    for (msItr = InValidHammingStr.begin(); msItr != InValidHammingStr.end(); msItr++) {
                        //cout << "Hamming Distance: " << msItr->first << ", String=" << msItr->second << endl;
                        substrPos1 = msItr->second.find("|");
                        HSidx1 = atoi((msItr->second.substr(0, substrPos1)).c_str());
                        substrPos2 = msItr->second.find("|", substrPos1 + 1);
                        HSidx2 = atoi((msItr->second.substr(substrPos1 + 1, substrPos2 - substrPos1 - 1)).c_str());
                        string tmpString = msItr->second.substr(substrPos2 + 1);
                        if (itrHamming->second[HSidx1] == tmpString[0] && itrHamming->second[HSidx2] == tmpString[1]) {
                            //cout << "String " << msItr->second << " matched pattern in " << msItr->second << endl;
                            isValid = false;
                            break;
                        }
                    }
                    if (!isValid) {
                        NoOfPrunedBitStrings++;
                        if (verbose)
                            std::cout << "Current testing Hamming string is invalid!" << std::endl;
                        if (verbose)
                            std::cout << "#pruned bit-strings: " << NoOfPrunedBitStrings << std::endl;

                        continue;
                    }
                    //end of Optimized Task 4

                    if (verbose)
                        std::cout << "Mininmal Order=" << minOrder << std::endl;

                    //prune current Hamming string if there exist two halfspaces that do not compatible
                    if (HammingDistance != itrHamming->first) {
                        if (verbose)
                            std::cout << "Intersecting halfspaces with hamming distance = " << itrHamming->first << std::endl;
                        HammingDistance = itrHamming->first;
                    }

                    //generate halfspaces for intersection based on hamming distance
                    char tmpBuf[64];
                    myitoa(itr1->second->nodeID, tmpBuf, 10);
                    strcpy(halfspaceFileName, namePrefix);
                    strcat(halfspaceFileName, tmpBuf);
                    strcat(halfspaceFileName, "_");
                    //strcat(halfspaceFileName,itrHamming->second);
                    //string tmpStr=string(itrHamming->second.begin(),itrHamming->second.end());
                    myitoa(rand(), tmpBuf, 10);
                    strcat(halfspaceFileName, tmpBuf);
                    strcat(halfspaceFileName, nameSuffix);
                    //cout << "Hamming string will be written to " << halfspaceFileName << endl;

                    bool nonZeroExtent = GenHammingHalfSpaces(halfspaceFileName, Dimen, itrHamming->second,
                                                              itr1->second->halfspaces, itr1->second->mbr);
                    if (nonZeroExtent == false) {
                        //cout << "Discard Hamming binstring " << itrHamming->second << ", for zero-extent! " << endl;
                        NoOfZeroExtentBinStrings++;
                        continue;
                    }

                    //perform intersection for all the halfspaces inside the set 'intersected halfspaces' of a leaf node
                    char sys_string[4096] = "bin\\qhalf Fp < ";
                    strcat(sys_string, halfspaceFileName);
                    strcat(sys_string, " | bin\\qconvex FA > ");
                    strcat(sys_string, volumeFilename);
                    system(sys_string);
                    if (verbose)
                        std::cout << "Command: " << sys_string << " performed..." << std::endl;

                    //open the file to read the volumes of the intersection of the halfspaces
                    std::ifstream fp_in(volumeFilename, std::ios::in);
                    string fileLine, volText;
                    size_t substrPos;
                    if (verbose)
                        std::cout << "Read and output cell volume:" << std::endl;
                    while (getline(fp_in, fileLine)) {
                        substrPos = fileLine.find("volume:");
                        if (substrPos != std::string::npos)
                            volText = fileLine.substr(substrPos + 7);
                    }
                    fp_in.close();

                    float volume = atof(volText.c_str());
                    if (verbose)
                        std::cout << "Volume=" << volume << std::endl;

                    string tmpString = halfspaceFileName;
                    FilesToRemove.push_back(tmpString);  //collect the files for removal later

                    if (volume > ZEROEXTENT)  //ZEROEXTENT: 1e-10
                    {
                        //store the min-cell found
                        if (minOrder >= (HammingDistance + NoOfCoveredHS)) {
                            minOrder = HammingDistance + NoOfCoveredHS;
                            std::set<std::shared_ptr<HalfSpace>> tmpSet;
                            std::copy((itr1->second->halfspaces).begin(),
                                      (itr1->second->halfspaces).end(),
                                      std::inserter(tmpSet, tmpSet.begin()));
                            minCellHalfSpaces.push_back(tmpSet);
                            vector<char> tmpVec = itrHamming->second;
                            binaryString.push_back(tmpVec);
                            if (verbose)
                                std::cout << "Node: " << itr1->second->nodeID << ", found a min-cell with binstring:"
                                     << string(itrHamming->second.begin(), itrHamming->second.end()) << std::endl;
                        } else {
                            NoOfDiscardedCells++;
                            stopIncrHammingDist = true;
                            string tmpString = halfspaceFileName;
                            FilesToRemove.push_back(tmpString);//collect the files for removal later
                            break;
                        }
                    } else {
                        NoOfZeroExtentBinStrings++;
                        if (verbose)
                            std::cout << "Node: " << itr1->second->nodeID << ", HammingDist: " << HammingDistance
                                 << " is empty!" << std::endl;
                    }
                    if (verbose)
                        std::cout << "#min-cells found so far:" << minCellHalfSpaces.size() << std::endl;
                }
                if (!stopIncrHammingDist) HammingDistance++;
                if (minOrder < INT_MAX && minOrder < (HammingDistance + NoOfCoveredHS))
                    break;
                else {
                    if (verbose)
                        std::cout << "Next, process bit-string with Hamming dist: " << HammingDistance << std::endl;
                    //getchar();
                }
            }

            if (FilesToRemove.size() >= 100) {
                vector<string>::iterator vItr = FilesToRemove.begin();
                while (!FilesToRemove.empty()) {
                    remove((*vItr).c_str());  //delete the halfspace data file
                    FilesToRemove.erase(vItr);
                    vItr = FilesToRemove.begin();
                }
            }
        }
        //if (minCellHalfSpaces.size()>0) NotFoundAllMinCells=false;
    }
    vector<string>::iterator vItr = FilesToRemove.begin();
    while (!FilesToRemove.empty()) {
        remove((*vItr).c_str());  //delete the halfspace data file
        FilesToRemove.erase(vItr);
        vItr = FilesToRemove.begin();
    }

    if (verbose) std::cout << "Number of invalid leaf nodes = " << NoOfInvalidLeaves << std::endl;
    if (verbose) std::cout << "#total bit-strings processed: " << NoOfBitStringsProcessed << std::endl;
    if (verbose) std::cout << "#pruned bit-strings: " << NoOfPrunedBitStrings << std::endl;
    if (verbose)
        std::cout << "#min-Cells found: " << minCellHalfSpaces.size() << ", minOrder=" << minOrder + NoOfCoveredHS << std::endl;


    if (minOrder == INT_MAX) return NoOfCoveredHS;
    return minOrder + NoOfCoveredHS;
}

