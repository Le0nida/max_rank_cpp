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

// *** NUOVO ***
// 1) Cancella albero esistente (root, macroRoots)
// 2) Macro-splitta MBR globale in subregioni
// 3) Partiziona i halfspaces (duplicandoli se servono) su tali subregioni
// 4) In parallelo, costruisce i sub-root con la solita insertHalfspaces
void QTree::inserthalfspacesMacroSplit(const std::vector<long int>& halfspaces) {
    if (halfspaces.empty()) return;

    // Distruggiamo l'albero old
    destroyAllNodes();

    // Definiamo l'MBR globale
    std::vector<std::array<double,2>> globalMBR(dims, {0.0, 1.0});

    // Generiamo subMBR
    auto subMBRs = macroSplitMBR(globalMBR);
    int nSub = subMBRs.size(); // di solito = 2^dims

    // Prepara i vettori di halfspaces per ogni subMBR
    std::vector< std::vector<long> > subHS(nSub);

    // Un unico pass per distribuire i halfspaces
    // Se un halfspace OverLAPPED con più subMBR, lo duplico in piu subHS[..]
    for (auto hsID : halfspaces) {
        auto hsPtr = halfspaceCache->get(hsID);
        if (!hsPtr) continue;

        // Controlla overlap con ciascun subMBR
        for (int i=0; i<nSub; i++) {
            // Esegui la logica "MbrVersusHalfSpace" su subMBR[i]
            // -> se OVERLAPPED => subHS[i].push_back(hsID)
            double minVal=0, maxVal=0;
            for (int d=0; d<dims; d++) {
                double c = hsPtr->coeff[d];
                if (c>=0) {
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

    // Ora costruiamo in parallelo i sotto-alberi
    macroRoots.resize(nSub, nullptr);

    // Esempio: std::async
    std::vector< std::future<QNode*> > futures;
    futures.reserve(nSub);

    for (int i=0; i<nSub; i++) {
        if (subHS[i].empty()) {
            // subMBR vuoto => nessun halfspace => potremmo lasciare null
            macroRoots[i] = nullptr;
            continue;
        }
        // lancio un task
        futures.push_back( std::async(std::launch::async, [this, i, &subHS, &subMBRs]() {
            // Costruisco un subtree => buildSubtree
            return buildSubtree(subMBRs[i], subHS[i]);
        }));
    }
    // Raccogliamo i risultati
    // NB: l'ordine dei future in "futures" deve corrispondere al subHS
    // Un semplice trucco: store (i, future) in un array
    int idxf=0;
    for (int i=0; i<nSub; i++) {
        if (subHS[i].empty()) continue;
        QNode* sr = futures[idxf++].get();
        macroRoots[i] = sr;
    }
}

// Suddivide [0,1]^dims in 2^dims subMBR
std::vector< std::vector<std::array<double,2>> >
QTree::macroSplitMBR(const std::vector<std::array<double,2>>& globalMBR) const
{
    std::vector< std::vector<std::array<double,2>> > ret;
    ret.resize( (1<<dims) );
    for (int mask=0; mask<(1<<dims); mask++) {
        ret[mask].resize(dims);
        for (int d=0; d<dims; d++) {
            double mid = 0.5*(globalMBR[d][0] + globalMBR[d][1]);
            if (mask & (1<<d)) {
                ret[mask][d] = {mid, globalMBR[d][1]};
            } else {
                ret[mask][d] = {globalMBR[d][0], mid};
            }
        }
    }
    return ret;
}

// Costruisce un subtree interno usando la logica del QNode
// E' un "mini-quadtree" con radice => new QNode(...) e poi node->insertHalfspaces(...).
// Se preferisci un "bulk build" più evoluto, puoi implementarlo qui.
// Per restare coerenti con la tua implementazione, usiamo insertHalfspaces.
QNode* QTree::buildSubtree(const std::vector<std::array<double,2>>& subMBR,
                           const std::vector<long>& subHS) const
{
    // Costruiamo un QNode con parent=null, e "owner=this" se vogliamo,
    // ma se vogliamo che each subroot gestisca le sue leaves in "this->leaves",
    // potremmo pass owner=this => cosi' "registerLeaf" funziona.
    QNode* rootNode = new QNode(const_cast<QTree*>(this), nullptr, subMBR);
    // Adesso, “inseriamo” tutti i halfspaces
    rootNode->insertHalfspaces(subHS);
    return rootNode;
}

// Ritorna tutte le leaves: in un “macro-split”, prendiamo leaves di macroRoots + eventuale root
std::vector<QNode*> QTree::getAllLeaves() const {
    std::vector<QNode*> all;
    // se vuoi includere pure la "old root"
    // BFS su root (se non nullo e non usato) => a volte potresti ignorarlo.

    // BFS su ciascun macroRoots
    for (auto sr : macroRoots) {
        if (!sr) continue;
        std::queue<QNode*> Q;
        Q.push(sr);
        while (!Q.empty()) {
            QNode* curr = Q.front(); Q.pop();
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
            QNode* curr = QQ.front(); QQ.pop();
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
            QNode* curr = QQ.front(); QQ.pop();
            for (auto c: curr->children) {
                if (!c) continue;
                c->setOrder(curr->getOrder() + c->covered.size());
                QQ.push(c);
            }
        }
    }
}
