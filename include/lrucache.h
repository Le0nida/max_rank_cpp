//
// Created by leona on 09/09/2024.
//

#ifndef LRUCACHE_H
#define LRUCACHE_H

#include <unordered_map>
#include <list>
#include <memory>  // Include per std::unique_ptr
#include <unordered_set>

#include "qnode.h"

class LRUCache {
private:
    // Map of nodeID to std::unique_ptr of QNode
    std::unordered_map<int, std::list<std::pair<int, std::shared_ptr<QNode>>>::iterator> cache;

    // LRU list to track node usage
    std::list<std::pair<int, std::shared_ptr<QNode>>> lruList;

    // Maximum cache size
    int cacheSize;

    // Contiene gli ID dei nodi "bloccati"
    std::unordered_set<int> lockedNodes;

public:
    explicit LRUCache(int size) : cacheSize(size) {}

    // Get a node by ID, load from disk if not present in cache
    std::shared_ptr<QNode> get(int nodeID);

    void add(std::shared_ptr<QNode> qnode);
    void invalidate(int nodeID);

    void lockNode(int nodeID);     // Blocca il nodo in cache
    void unlockNode(int nodeID);   // Sblocca il nodo in cache

    // Helper function to generate file path for node serialization
    std::string getFilePath(int nodeID);

    // Save the current cache state to disk
    void saveCacheToDisk();
};

// Dichiarazione della cache globale
extern LRUCache globalCache;

#endif // LRUCACHE_H
