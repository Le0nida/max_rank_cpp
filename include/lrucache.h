//
// Created by leona on 09/09/2024.
//

#ifndef LRUCACHE_H
#define LRUCACHE_H

#include <iostream>
#include <unordered_map>
#include <list>
#include <memory>
#include <mutex>
#include <unordered_set>
#include "qnode.h"
#include "utils.h"
#include <leveldb/db.h>

class LRUCache {
private:
    // Map of nodeID to iterator in lruList
    std::unordered_map<int, std::list<std::pair<int, std::shared_ptr<QNode>>>::iterator> cache;

    // LRU list to track node usage
    std::list<std::pair<int, std::shared_ptr<QNode>>> lruList;

    // Maximum cache size
    int cacheSize;

    // Contains IDs of "locked" nodes
    std::unordered_set<int> lockedNodes;

    // LevelDB database pointer
    leveldb::DB* db;

    // Mutex for thread safety
    std::mutex cacheMutex;

    // Helper methods
    void saveNodeToDB(int nodeID, const std::shared_ptr<QNode>& node);
    std::shared_ptr<QNode> loadNodeFromDB(int nodeID);

public:
    explicit LRUCache();
    ~LRUCache();

    // Get a node by ID, load from DB if not present in cache
    std::shared_ptr<QNode> get(int nodeID);

    void add(std::shared_ptr<QNode> qnode);
    void invalidate(int nodeID);

    void lockNode(int nodeID);     // Lock the node in cache
    void unlockNode(int nodeID);   // Unlock the node in cache

    void cleanup(); // Not needed with LevelDB
};

// Declaration of global cache
extern LRUCache globalCache;

#endif // LRUCACHE_H
