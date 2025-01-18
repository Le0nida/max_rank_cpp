#include "qtree.h"

QTree::QTree(const int dims, const int maxhsnode, const int maxLevel)
    : dims(dims),
      maxhsnode(maxhsnode),
      maxLevel(maxLevel),
      root(nullptr)
{

    // Create the classical root covering [0,1]^dims
    root = createroot();

    // Precompute sub-MBRs for macro-split
    const std::vector<std::array<double,2>> globalMBR(dims, {0.0, 1.0});
    precomputedSubMBRs = macroSplitMBR(globalMBR);

    // Prepare macroRoots, one for each sub-MBR
    macroRoots.resize(numOfSubdivisions, nullptr);
}

QTree::~QTree() {
    destroyAllNodes();
}

QNode* QTree::createroot() {
    // Root node with MBR = [0,1]^dims
    const std::vector<std::array<double,2>> mbr(dims, {0.0, 1.0});
    return new QNode(this, nullptr, mbr, 0);
}

void QTree::destroyAllNodes() {
    // 1) Destroy the classical root subtree
    if (root) {
        std::queue<QNode*> q;
        q.push(root);
        while (!q.empty()) {
            QNode* current = q.front();
            q.pop();
            for (auto child : current->children) {
                if (child) q.push(child);
            }
            delete current;
        }
        root = nullptr;
    }

    // 2) Destroy each macro-root subtree
    for (auto sr : macroRoots) {
        if (!sr) continue;
        std::queue<QNode*> q;
        q.push(sr);
        while (!q.empty()) {
            QNode* current = q.front();
            q.pop();
            for (auto child : current->children) {
                if (child) q.push(child);
            }
            delete current;
        }
    }
    macroRoots.clear();
}

std::vector< std::vector<std::array<double,2>> >
QTree::macroSplitMBR(const std::vector<std::array<double,2>>& globalMBR) const {
    // Each bit in 'mask' picks lower or upper half for each dimension
    std::vector< std::vector<std::array<double,2>> > result(1 << dims);

    for (int mask = 0; mask < (1 << dims); mask++) {
        result[mask].resize(dims);
        double mid;
        for (int d = 0; d < dims; d++) {
            mid = 0.5 * (globalMBR[d][0] + globalMBR[d][1]);
            if (mask & (1 << d)) {
                result[mask][d] = { mid, globalMBR[d][1] };
            } else {
                result[mask][d] = { globalMBR[d][0], mid };
            }
        }
    }
    return result;
}

QNode* QTree::buildSubtree(const std::vector<std::array<double,2>>& subMBR,
                           const std::vector<long>& subHS) const {
    // Create a new root node for this subMBR
    auto* rootNode = new QNode(const_cast<QTree*>(this), nullptr, subMBR, 1);
    // Insert halfspaces incrementally
    rootNode->insertHalfspaces(subHS);
    return rootNode;
}

void QTree::inserthalfspacesMacroSplit(const std::vector<long int>& halfspaces) {
    if (halfspaces.empty()) return;

    // Prepare containers
    int nSub = static_cast<int>(precomputedSubMBRs.size());
    std::vector< std::vector<long> > subHS(nSub);

    // Multi-threaded distribution
    size_t totalHS = halfspaces.size();
    unsigned int hwThreads = std::max(1U, std::thread::hardware_concurrency());
    size_t chunkSize = (totalHS + hwThreads - 1) / hwThreads;

    // Partial results per thread
    std::vector<std::vector<std::vector<long>>> partialRes(
        hwThreads, std::vector<std::vector<long>>(nSub)
    );

    // 1) Distribute halfspaces among subMBRs in parallel
    std::vector<std::future<void>> distributionFutures;
    distributionFutures.reserve(hwThreads);

    for (unsigned int t = 0; t < hwThreads; t++) {
        size_t start = t * chunkSize;
        if (start >= totalHS) break;
        size_t end = std::min(start + chunkSize, totalHS);

        // Parallel task: assign halfspaces to subMBRs
        distributionFutures.push_back(std::async(std::launch::async,
            [this, &halfspaces, start, end, &partialRes, t, nSub]()
        {
            for (size_t idx = start; idx < end; ++idx) {
                long hsID = halfspaces[idx];
                const auto hsPtr = halfspaceCache->get(hsID);
                if (!hsPtr) continue;

                double known = hsPtr->known;
                const auto& cVec = hsPtr->coeff;
                // Evaluate each precomputed subMBR
                for (int i = 0; i < nSub; i++) {
                    double minVal = 0, maxVal = 0;
                    for (size_t d = 0; d < cVec.size(); d++) {
                        double c = cVec[d];
                        if (c >= 0) {
                            minVal += c * precomputedSubMBRs[i][d][0];
                            maxVal += c * precomputedSubMBRs[i][d][1];
                        } else {
                            minVal += c * precomputedSubMBRs[i][d][1];
                            maxVal += c * precomputedSubMBRs[i][d][0];
                        }
                    }
                    // Classical Overlapped/Below/Above logic
                    if (maxVal < known || minVal <= known) {
                        // Entirely below or Overlapped => add it
                        partialRes[t][i].push_back(hsID);
                    }
                }
            }
        }));
    }

    // Wait for distribution tasks to finish
    for (auto &f : distributionFutures) {
        f.get();
    }

    // Combine partial results
    for (unsigned int t = 0; t < hwThreads; t++) {
        for (int i = 0; i < nSub; i++) {
            if (!partialRes[t][i].empty()) {
                subHS[i].insert(subHS[i].end(),
                                partialRes[t][i].begin(),
                                partialRes[t][i].end());
            }
        }
    }

    // 2) Build/update sub-roots in parallel
    std::vector<std::future<void>> buildFutures;
    buildFutures.reserve(nSub);

    for (int i = 0; i < nSub; i++) {
        if (subHS[i].empty()) continue;

        buildFutures.push_back(std::async(std::launch::async,
            [this, i, &subHS]()
        {
            if (!macroRoots[i]) {
                // Create a new sub-root
                macroRoots[i] = buildSubtree(precomputedSubMBRs[i], subHS[i]);
            } else {
                // Update existing sub-root
                macroRoots[i]->insertHalfspaces(subHS[i]);
            }
        }));
    }

    // Wait for all sub-root builds
    for (auto &f : buildFutures) {
        f.get();
    }
}

std::vector<QNode*> QTree::getAllLeaves() const {
    std::vector<QNode*> result;
    result.reserve(128);

    // Collect leaves from each macro-root
    for (auto sr : macroRoots) {
        if (!sr) continue;

        std::queue<QNode*> q;
        q.push(sr);

        while (!q.empty()) {
            QNode* curr = q.front();
            q.pop();
            if (curr->leaf) {
                result.push_back(curr);
            } else {
                for (auto c : curr->children) {
                    if (c) q.push(c);
                }
            }
        }
    }
    return result;
}

void QTree::updateAllOrders() {
    // 1) Update root subtree if it exists
    if (root) {
        root->order = root->covered.size();
        std::queue<QNode*> queueNodes;
        queueNodes.push(root);
        while (!queueNodes.empty()) {
            QNode* curr = queueNodes.front();
            queueNodes.pop();
            for (auto c : curr->children) {
                if (c) {
                    size_t newOrder = curr->order + c->covered.size();
                    c->order = newOrder;
                    queueNodes.push(c);
                }
            }
        }
    }

    // 2) Update each macro-root
    for (auto sr : macroRoots) {
        if (!sr) continue;
        sr->order = sr->covered.size();
        std::queue<QNode*> queueNodes;
        queueNodes.push(sr);
        while (!queueNodes.empty()) {
            QNode* curr = queueNodes.front();
            queueNodes.pop();
            for (auto c : curr->children) {
                if (c) {
                    size_t newOrder = curr->order + c->covered.size();
                    c->order = newOrder;
                    queueNodes.push(c);
                }
            }
        }
    }
}
