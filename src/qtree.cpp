#include "qtree.h"
#include "qnode.h"
#include "geom.h"
#include "halfspace.h"
#include <iostream>
#include <queue>
#include <future>   // per std::async, se vuoi

extern int numOfSubdivisions;

QTree::QTree(int dims, int maxhsnode)
    : dims(dims), maxhsnode(maxhsnode), root(nullptr)
{
    numOfSubdivisions = (1 << dims); // 2^dims
    root = createroot();
}

QTree::~QTree() {
    destroyAllNodes();
    root = nullptr;
}

QNode* QTree::createroot() {
    // MBR [0,1]^dims
    std::vector<std::array<double,2>> mbr(dims, {0.0, 1.0});
    QNode* r = new QNode(this, nullptr, mbr);
    return r;
}

void QTree::destroyAllNodes() {
    if (!root) return;
    // BFS su root
    {
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
    if (halfspaces.empty()) return;
    if (!root) return;
    root->insertHalfspaces(halfspaces);
}

// *** MODIFICA PRINCIPALE ***
// Ora *non* distruggiamo più l'albero esistente: riutilizziamo o creiamo i sub-root necessari.
void QTree::inserthalfspacesMacroSplit(const std::vector<long int>& halfspaces) {
    if (halfspaces.empty()) return;

    // NON distruggiamo l'albero esistente:
    // (rimosso) //destroyAllNodes();

    // Definiamo l'MBR globale
    std::vector<std::array<double,2>> globalMBR(dims, {0.0, 1.0});

    // Generiamo subMBR
    auto subMBRs = macroSplitMBR(globalMBR);
    int nSub = subMBRs.size(); // di solito = 2^dims

    // Prepara i vettori di halfspaces per ogni subMBR
    std::vector< std::vector<long> > subHS(nSub);

    // Distribuiamo i nuovi halfspaces tra i vari subMBR
    for (auto hsID : halfspaces) {
        auto hsPtr = halfspaceCache->get(hsID);
        if (!hsPtr) continue;

        for (int i = 0; i < nSub; i++) {
            double minVal = 0, maxVal = 0;
            for (int d = 0; d < dims; d++) {
                double c = hsPtr->coeff[d];
                if (c >= 0) {
                    minVal += c * subMBRs[i][d][0];
                    maxVal += c * subMBRs[i][d][1];
                } else {
                    minVal += c * subMBRs[i][d][1];
                    maxVal += c * subMBRs[i][d][0];
                }
            }
            if (maxVal < hsPtr->known) {
                // BELOW => potremmo ignorare
                subHS[i].push_back(hsID);
            } else if (minVal > hsPtr->known) {
                // ABOVE => ignoriamo
            } else {
                // Overlapped
                subHS[i].push_back(hsID);
            }
        }
    }

    // Assicuriamoci che macroRoots abbia la taglia giusta
    if ((int)macroRoots.size() < nSub) {
        macroRoots.resize(nSub, nullptr);
    }

    // Costruzione/aggiornamento parallelo dei sotto-alberi
    std::vector<std::future<void>> futures;
    futures.reserve(nSub);

    for (int i = 0; i < nSub; i++) {
        if (subHS[i].empty()) {
            // Nessun halfspace per questa partizione: non c'è nulla da aggiungere
            continue;
        }
        futures.push_back(std::async(std::launch::async, [this, i, &subHS, &subMBRs]() {
            if (!macroRoots[i]) {
                // Creiamo un nuovo subtree se non esiste
                macroRoots[i] = buildSubtree(subMBRs[i], subHS[i]);
            } else {
                // Altrimenti inseriamo semplicemente i nuovi halfspaces
                macroRoots[i]->insertHalfspaces(subHS[i]);
            }
        }));
    }

    // Aspettiamo che tutti i thread abbiano terminato
    for (auto &f : futures) {
        f.get();
    }
}

// Suddivide [0,1]^dims in 2^dims subMBR
std::vector< std::vector<std::array<double,2>> >
QTree::macroSplitMBR(const std::vector<std::array<double,2>>& globalMBR) const
{
    std::vector< std::vector<std::array<double,2>> > ret;
    ret.resize( (1<<dims) );
    for (int mask = 0; mask < (1 << dims); mask++) {
        ret[mask].resize(dims);
        for (int d = 0; d < dims; d++) {
            double mid = 0.5 * (globalMBR[d][0] + globalMBR[d][1]);
            if (mask & (1 << d)) {
                ret[mask][d] = {mid, globalMBR[d][1]};
            } else {
                ret[mask][d] = {globalMBR[d][0], mid};
            }
        }
    }
    return ret;
}

// Costruisce un subtree con la logica incrementale
QNode* QTree::buildSubtree(const std::vector<std::array<double,2>>& subMBR,
                           const std::vector<long>& subHS) const
{
    QNode* rootNode = new QNode(const_cast<QTree*>(this), nullptr, subMBR);
    rootNode->insertHalfspaces(subHS);
    return rootNode;
}

// Ritorna tutte le leaves dai macroRoots
std::vector<QNode*> QTree::getAllLeaves() const {
    std::vector<QNode*> all;
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
    // se usi ancora root
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

    // e poi su macroRoots
    for (auto sr : macroRoots) {
        if (!sr) continue;
        sr->setOrder(sr->covered.size());
        std::queue<QNode*> QQ;
        QQ.push(sr);
        while (!QQ.empty()) {
            QNode* curr = QQ.front();
            QQ.pop();
            for (auto c: curr->children) {
                if (!c) continue;
                c->setOrder(curr->getOrder() + c->covered.size());
                QQ.push(c);
            }
        }
    }
}
