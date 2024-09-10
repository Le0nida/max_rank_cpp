//
// Created by leona on 09/09/2024.
//

#ifndef LRUCACHE_H
#define LRUCACHE_H

#include <unordered_map>
#include <list>
#include <memory>  // Include per std::unique_ptr
#include "qnode.h"

class LRUCache {
private:
    // Map of nodeID to std::unique_ptr of QNode
    std::unordered_map<int, std::unique_ptr<QNode>> cache;

    // LRU list to track node usage
    std::list<int> lruList;

    // Maximum cache size
    int cacheSize;

public:
    LRUCache(int size) : cacheSize(size) {}

    // Get a node by ID, load from disk if not present in cache
    QNode* get(int nodeID);

    // Helper function to generate file path for node serialization
    std::string getFilePath(int nodeID);

    // Save the current cache state to disk
    void saveCacheToDisk();
};

#endif // LRUCACHE_H
