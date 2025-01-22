#include "qtree.h"

QTree::QTree(const int dims, const int maxhsnode, const int maxLevel)
    : dims(dims),
      maxhsnode(maxhsnode),
      maxLevel(maxLevel),
      root(nullptr)
{
    // Create the classical root covering [0,1]^dims
    root = createroot();

    // Precompute sub-MBRs for macro-split (using float)
    const std::vector<std::array<float,2>> globalMBR(dims, {0.0f, 1.0f});
    precomputedSubMBRs = macroSplitMBR(globalMBR);

    // Prepare macroRoots, one for each sub-MBR
    macroRoots.resize(numOfSubdivisions, nullptr);
}

QTree::~QTree() {
    destroyAllNodes();
}

QNode* QTree::createroot() {
    // Root node with MBR = [0,1]^dims (float)
    const std::vector<std::array<float,2>> mbr(dims, {0.0f, 1.0f});
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

std::vector< std::vector<std::array<float,2>> >
QTree::macroSplitMBR(const std::vector<std::array<float,2>>& globalMBR) const
{
    // Each bit in 'mask' picks lower or upper half for each dimension
    std::vector< std::vector<std::array<float,2>> > result(1 << dims);

    for (int mask = 0; mask < (1 << dims); mask++) {
        result[mask].resize(dims);
        for (int d = 0; d < dims; d++) {
            float minVal = globalMBR[d][0];
            float maxVal = globalMBR[d][1];
            float mid    = 0.5f * (minVal + maxVal);
            if (mask & (1 << d)) {
                result[mask][d] = { mid, maxVal };
            } else {
                result[mask][d] = { minVal, mid };
            }
        }
    }
    return result;
}

QNode* QTree::buildSubtree(const std::vector<std::array<float,2>>& subMBR,
                           const std::vector<long>& subHS) const
{
    // Create a new root node for this subMBR
    auto* rootNode = new QNode(const_cast<QTree*>(this), nullptr, subMBR, 1);
    // Insert halfspaces incrementally
    rootNode->insertHalfspaces(subHS);
    return rootNode;
}

void QTree::inserthalfspacesMacroSplit(const std::vector<long int>& halfspaces) {
    if (halfspaces.empty()) return;

    const int nSub = (int) precomputedSubMBRs.size();
    // For each subMBR, store two lists: fullyCovered, partialOverlapped
    std::vector<std::vector<long>> fullyCovered(nSub);
    std::vector<std::vector<long>> partialOverlapped(nSub);

    const size_t totalHS = halfspaces.size();
    const unsigned int hwThreads = std::max(1U, std::thread::hardware_concurrency());
    const size_t chunkSize = (totalHS + hwThreads - 1) / hwThreads;

    // partial results per thread
    struct HSBatch {
        std::vector<long> fully;
        std::vector<long> partial;
    };
    std::vector<std::vector<HSBatch>> partialRes(hwThreads, std::vector<HSBatch>(nSub));

    // 1) Distribute halfspaces in parallel
    std::vector<std::future<void>> distributionFutures;
    distributionFutures.reserve(hwThreads);

    for (unsigned int t = 0; t < hwThreads; t++) {
        size_t start = t * chunkSize;
        if (start >= totalHS) break;
        size_t end = std::min(start + chunkSize, totalHS);

        distributionFutures.push_back(std::async(std::launch::async,
            [this, &halfspaces, start, end, &partialRes, t, nSub]()
        {
            for (size_t idx = start; idx < end; ++idx) {
                long hsID = halfspaces[idx];
                auto hsPtr = halfspaceCache->get(hsID);
                if (!hsPtr) continue;

                double known = hsPtr->known;           // We keep halfspace in double
                const auto& cVec = hsPtr->coeff;       // also double

                // For each subMBR
                for (int i = 0; i < nSub; i++) {
                    // We'll do minVal, maxVal in double for numeric stability
                    double minVal = 0, maxVal = 0;
                    for (size_t d = 0; d < cVec.size(); d++) {
                        double c = cVec[d];
                        float low = precomputedSubMBRs[i][d][0];
                        float high= precomputedSubMBRs[i][d][1];
                        if (c >= 0) {
                            minVal += c * low;
                            maxVal += c * high;
                        } else {
                            minVal += c * high;
                            maxVal += c * low;
                        }
                    }

                    // "Fully covered" => store once in fully list
                    if (maxVal < known) {
                        partialRes[t][i].fully.push_back(hsID);
                    }
                    // partial overlap => store in partial list
                    else if (minVal <= known) {
                        partialRes[t][i].partial.push_back(hsID);
                    }
                    // else => skip
                }
            }
        }));
    }

    // Wait distribution
    for (auto &f : distributionFutures) {
        f.get();
    }

    // Combine partial results
    for (unsigned int t = 0; t < hwThreads; t++) {
        for (int i = 0; i < nSub; i++) {
            if (!partialRes[t][i].fully.empty()) {
                fullyCovered[i].insert(
                    fullyCovered[i].end(),
                    partialRes[t][i].fully.begin(),
                    partialRes[t][i].fully.end()
                );
            }
            if (!partialRes[t][i].partial.empty()) {
                partialOverlapped[i].insert(
                    partialOverlapped[i].end(),
                    partialRes[t][i].partial.begin(),
                    partialRes[t][i].partial.end()
                );
            }
        }
    }

    // 2) Build/Update sub-roots
    std::vector<std::future<void>> buildFutures;
    buildFutures.reserve(nSub);

    for (int i = 0; i < nSub; i++) {
        // if no halfspace at all => skip
        if (fullyCovered[i].empty() && partialOverlapped[i].empty()) continue;

        // **Memory saver**:
        // If partialOverlapped[i] is empty => entire subMBR is fully covered,
        // no need for a macro-root subtree.
        if (partialOverlapped[i].empty()) {
            // Optionally store these "fully covered" halfspaces in a global structure
            // if you need them for later queries. Or just skip creation:
            continue;
        }

        // Otherwise, partial coverage => we do need a macro-root
        buildFutures.push_back(std::async(std::launch::async, [this, i, &fullyCovered, &partialOverlapped]() {
            if (!macroRoots[i]) {
                // Create a sub-root
                QNode* subRoot = new QNode(const_cast<QTree*>(this),
                                           nullptr,
                                           precomputedSubMBRs[i],
                                           1);

                // "fully covered" => put in subRoot->covered
                // (Delta coverage approach: these are newly discovered coverage at this root)
                subRoot->covered.insert(subRoot->covered.end(),
                                        fullyCovered[i].begin(),
                                        fullyCovered[i].end());

                // partial => insert halfspaces
                subRoot->insertHalfspaces(partialOverlapped[i]);

                macroRoots[i] = subRoot;
            } else {
                // existing subRoot =>
                // We add "fully covered" to subRoot->covered (delta coverage at this level)
                macroRoots[i]->covered.insert(
                    macroRoots[i]->covered.end(),
                    fullyCovered[i].begin(),
                    fullyCovered[i].end()
                );
                // partial => insert
                macroRoots[i]->insertHalfspaces(partialOverlapped[i]);
            }
        }));
    }

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
