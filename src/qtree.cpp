#include "qtree.h"
#include "qnode.h"
#include "geom.h"
#include "halfspace.h"
#include <iostream>
#include <queue>
#include <future>
#include <thread>    // std::thread::hardware_concurrency
#include <algorithm> // std::min

extern int numOfSubdivisions;

QTree::QTree(int dims, int maxhsnode)
    : dims(dims), maxhsnode(maxhsnode), root(nullptr)
{
    // Calcoliamo 2^dims
    numOfSubdivisions = (1 << dims);

    // Creiamo la root "classica"
    root = createroot();

    // Pre-calcoliamo i subMBR [0,1]^dims una volta sola
    std::vector<std::array<double,2>> globalMBR(dims, {0.0, 1.0});
    precomputedSubMBRs = macroSplitMBR(globalMBR);

    // Inizializziamo macroRoots
    macroRoots.resize(numOfSubdivisions, nullptr);
}

QTree::~QTree() {
    destroyAllNodes();
    root = nullptr;
}

// Crea una root con MBR unitaria [0,1]^dims
QNode* QTree::createroot() {
    std::vector<std::array<double,2>> mbr(dims, {0.0, 1.0});
    QNode* r = new QNode(this, nullptr, mbr);
    return r;
}

// Distrugge root e tutti i macroRoots
void QTree::destroyAllNodes() {
    if (root) {
        // BFS su root
        std::queue<QNode*> Q;
        Q.push(root);
        while (!Q.empty()) {
            QNode* current = Q.front();
            Q.pop();
            for (auto child : current->children) {
                if (child) Q.push(child);
            }
            delete current;
        }
        leaves.clear();
        root = nullptr;
    }
    // Distruggi i subroot se esistono
    for (auto sr : macroRoots) {
        if (!sr) continue;
        std::queue<QNode*> Q;
        Q.push(sr);
        while (!Q.empty()) {
            QNode* current = Q.front();
            Q.pop();
            for (auto child : current->children) {
                if (child) Q.push(child);
            }
            delete current;
        }
    }
    macroRoots.clear();
}

// Inserimento incrementale classico
void QTree::inserthalfspaces(const std::vector<long int>& halfspaces) {
    if (halfspaces.empty() || !root) return;
    root->insertHalfspaces(halfspaces);
}

/**
 * Inserimento parallelo via macro-split:
 *  1) Distribuzione parallela dei nuovi halfspaces nei subMBR
 *  2) Creazione/aggiornamento parallelo dei relativi sotto-alberi
 */
void QTree::inserthalfspacesMacroSplit(const std::vector<long int>& halfspaces) {
    if (halfspaces.empty()) return;

    // Preleviamo i subMBR pre-calcolati
    const auto& subMBRs = precomputedSubMBRs;
    const int nSub = static_cast<int>(subMBRs.size());

    // Vettore di vettori su cui accumulare gli halfspaces per ogni subMBR
    // Inizializziamo già alla dimensione nSub per evitare push_back su subHS
    std::vector< std::vector<long> > subHS(nSub);

    // 1) PARALLELIZZIAMO la distribuzione: suddividiamo i halfspaces in "chunk"
    //    e in ciascun chunk distribuiamo i halfspaces nei subMBR (in parallelo).
    const size_t totalHS = halfspaces.size();
    const unsigned int hwThreads = std::max(1U, std::thread::hardware_concurrency());
    const size_t chunkSize = (totalHS + hwThreads - 1) / hwThreads;
        // divisione "ceil" tra i thread

    // Vettore "parziale" di dimensione #thread, ognuno contiene un subHS "locale"
    // da combinare alla fine
    std::vector< std::vector< std::vector<long> > > partialRes(hwThreads,
                                   std::vector<std::vector<long>>(nSub));

    std::vector<std::future<void>> distributionFutures;
    distributionFutures.reserve(hwThreads);

    for (unsigned int t = 0; t < hwThreads; t++) {
        const size_t start = t * chunkSize;
        if (start >= totalHS) break; // niente più chunk
        const size_t end = std::min(start + chunkSize, totalHS);

        // Avviamo un task asincrono per l'intervallo [start, end)
        distributionFutures.push_back(std::async(std::launch::async,
            [this, &halfspaces, &subMBRs, start, end, &partialRes, t]()
        {
            for (size_t idx = start; idx < end; ++idx) {
                long hsID = halfspaces[idx];
                auto hsPtr = halfspaceCache->get(hsID);
                if (!hsPtr) continue;

                const double known = hsPtr->known;
                // Ciclo su i subMBR
                const int localNSub = (int)subMBRs.size();
                for (int i = 0; i < localNSub; i++) {
                    double minVal = 0, maxVal = 0;
                    // Sfruttiamo i coefficienti
                    const auto& cVec = hsPtr->coeff;
                    for (size_t d = 0; d < cVec.size(); d++) {
                        double c = cVec[d];
                        // Più veloce di if-else? Forse un if qui è meglio di branching inevitabile
                        if (c >= 0) {
                            minVal += c * subMBRs[i][d][0];
                            maxVal += c * subMBRs[i][d][1];
                        } else {
                            minVal += c * subMBRs[i][d][1];
                            maxVal += c * subMBRs[i][d][0];
                        }
                    }
                    // Classica logica di Overlapped/Below/Above
                    if (maxVal < known) {
                        // BELOW
                        partialRes[t][i].push_back(hsID);
                    } else if (minVal <= known) {
                        // Overlapped
                        partialRes[t][i].push_back(hsID);
                    }
                    // Se minVal > known => Above => ignoro
                }
            }
        }));
    }

    // Attendo la fine di tutti i task di distribuzione
    for (auto &f : distributionFutures) {
        f.get();
    }

    // 1b) Combiniamo i risultati parziali in subHS[i]
    for (unsigned int t = 0; t < hwThreads; t++) {
        for (int i = 0; i < nSub; i++) {
            auto &localVec = partialRes[t][i];
            if (!localVec.empty()) {
                // Evitiamo troppi insert uno a uno
                subHS[i].insert(subHS[i].end(), localVec.begin(), localVec.end());
            }
        }
    }

    // 2) PARALLELIZZIAMO la costruzione/aggiornamento dei sub-root
    //    (logica identica a prima, ma a chunk di subMBR)
    std::vector<std::future<void>> buildFutures;
    buildFutures.reserve(nSub);

    for (int i = 0; i < nSub; i++) {
        if (subHS[i].empty()) {
            // Nessun halfspace per questa partizione => niente da fare
            continue;
        }
        // Avviamo in parallelo la build o l'aggiornamento
        buildFutures.push_back(std::async(std::launch::async, [this, i, &subHS, &subMBRs]() {
            if (!macroRoots[i]) {
                // Non esiste => costruiamo
                macroRoots[i] = buildSubtree(subMBRs[i], subHS[i]);
            } else {
                // Esiste => aggiorniamo con i nuovi halfspaces
                macroRoots[i]->insertHalfspaces(subHS[i]);
            }
        }));
    }

    // Aspettiamo tutti i sub-root
    for (auto &f : buildFutures) {
        f.get();
    }
}

// Pre-calcola le 2^dims suddivisioni di [0,1]^dims
std::vector< std::vector<std::array<double,2>> >
QTree::macroSplitMBR(const std::vector<std::array<double,2>>& globalMBR) const
{
    std::vector< std::vector<std::array<double,2>> > ret;
    ret.resize( (1 << dims) );

    for (int mask = 0; mask < (1 << dims); mask++) {
        ret[mask].resize(dims);
        for (int d = 0; d < dims; d++) {
            double mid = 0.5 * (globalMBR[d][0] + globalMBR[d][1]);
            if (mask & (1 << d)) {
                ret[mask][d] = { mid, globalMBR[d][1] };
            } else {
                ret[mask][d] = { globalMBR[d][0], mid };
            }
        }
    }
    return ret;
}

// Crea un sotto-albero con la logica incrementale
QNode* QTree::buildSubtree(const std::vector<std::array<double,2>>& subMBR,
                           const std::vector<long>& subHS) const
{
    // Creiamo un QNode radice di sub-albero
    QNode* rootNode = new QNode(const_cast<QTree*>(this), nullptr, subMBR);
    // Inserimento incrementale
    rootNode->insertHalfspaces(subHS);
    return rootNode;
}

// Ritorna tutte le leaves (macro-roots)
std::vector<QNode*> QTree::getAllLeaves() const {
    std::vector<QNode*> all;
    all.reserve(128); // micro-ottimizzazione: riserva un po' di spazio

    for (auto sr : macroRoots) {
        if (!sr) continue;
        std::queue<QNode*> Q;
        Q.push(sr);
        while (!Q.empty()) {
            QNode* curr = Q.front();
            Q.pop();
            if (curr->isLeaf()) {
                all.push_back(curr);
            } else {
                for (auto c : curr->children) {
                    if (c) Q.push(c);
                }
            }
        }
    }
    return all;
}

// Esegue updateAllOrders su root + su tutti i subroot
void QTree::updateAllOrders() {
    // Se vogliamo aggiornare anche la root "classica"
    if (root) {
        root->setOrder(root->covered.size());
        std::queue<QNode*> QQ;
        QQ.push(root);
        while (!QQ.empty()) {
            QNode* curr = QQ.front();
            QQ.pop();
            for (auto c: curr->children) {
                if (c) {
                    size_t newOrder = curr->getOrder() + c->covered.size();
                    c->setOrder(newOrder);
                    QQ.push(c);
                }
            }
        }
    }

    // Poi aggiorniamo i macroRoots
    for (auto sr : macroRoots) {
        if (!sr) continue;
        sr->setOrder(sr->covered.size());
        std::queue<QNode*> QQ;
        QQ.push(sr);
        while (!QQ.empty()) {
            QNode* curr = QQ.front();
            QQ.pop();
            for (auto c: curr->children) {
                if (c) {
                    size_t newOrder = curr->getOrder() + c->covered.size();
                    c->setOrder(newOrder);
                    QQ.push(c);
                }
            }
        }
    }
}
