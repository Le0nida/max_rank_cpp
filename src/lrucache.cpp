#include "lrucache.h"
#include <leveldb/write_batch.h>
#include <sstream>
#include <iostream>
#include <thread>

// Definition of global cache
LRUCache globalCache;

LRUCache::LRUCache() {
    const int MIN_CACHE_SIZE = 1000;
    const int MAX_CACHE_SIZE = 1000000;

    // Calculate cache size based on available memory
    size_t availableMemory = getAvailableMemory();
    if (availableMemory == 0) {
        std::cerr << "Failed to retrieve available memory. Using default cache size of 10,000." << std::endl;
        cacheSize = 10000;
    } else {
        // Use 10% of available memory for cache
        size_t cacheMemoryUsage = availableMemory / 10;
        size_t estimatedNodeSize = sizeof(QNode) + 1024; // Adjust as necessary
        cacheSize = static_cast<int>(cacheMemoryUsage / estimatedNodeSize);

        if (cacheSize < MIN_CACHE_SIZE) {
            cacheSize = MIN_CACHE_SIZE;
        } else if (cacheSize > MAX_CACHE_SIZE) {
            cacheSize = MAX_CACHE_SIZE;
        }
        std::cout << "Cache size set to " << cacheSize << " based on available memory." << std::endl;
    }

    cacheSize = 8000;

    // Open LevelDB database
    leveldb::Options options;
    options.create_if_missing = true;
    leveldb::Status status = leveldb::DB::Open(options, "nodesdb", &db);
    if (!status.ok()) {
        std::cerr << "Unable to open/create LevelDB database: " << status.ToString() << std::endl;
        exit(EXIT_FAILURE);
    }
}

LRUCache::~LRUCache() {
    for (auto& [nodeID, nodeIter] : cache) {
        auto& node = nodeIter->second;
        saveNodeToDB(nodeID, node);
    }

    delete db;
}

void LRUCache::saveNodeToDB(int nodeID, const std::shared_ptr<QNode>& node) {
    std::stringstream ss;
    node->serialize(ss);
    leveldb::Status status = db->Put(leveldb::WriteOptions(), std::to_string(nodeID), ss.str());
    if (!status.ok()) {
        std::cerr << "Failed to write node " << nodeID << " to LevelDB: " << status.ToString() << std::endl;
    }
}

std::shared_ptr<QNode> LRUCache::loadNodeFromDB(int nodeID) {
    std::string data;
    leveldb::Status status = db->Get(leveldb::ReadOptions(), std::to_string(nodeID), &data);
    if (!status.ok()) {
        std::cerr << "Node with ID " << nodeID << " not found in LevelDB." << std::endl;
        return nullptr;
    }
    std::stringstream ss(data);
    auto node = std::make_shared<QNode>();
    node->deserialize(ss);
    return node;
}

std::shared_ptr<QNode> LRUCache::get(int nodeID) {
    {
        std::lock_guard<std::mutex> lock(cacheMutex); // Acquire lock only for cache access
        if (cache.find(nodeID) != cache.end()) {
            lruList.splice(lruList.begin(), lruList, cache[nodeID]);
            return cache[nodeID]->second;
        }
    }

    // Load node from DB outside of the locked section
    auto node = loadNodeFromDB(nodeID);
    if (node) {
        add(node); // Add back to cache
        return node;
    } else {
        std::cerr << "Node with ID " << nodeID << " not found in LevelDB." << std::endl;
        return nullptr;
    }
}

void LRUCache::add(std::shared_ptr<QNode> qnode) {
    int nodeID = qnode->getNodeID();

    // Invalidate the node before locking to avoid deadlock
    invalidate(nodeID); // No lock inside invalidate()

    {
        std::lock_guard<std::mutex> lock(cacheMutex); // Acquire lock for cache modification
        if (cache.size() >= cacheSize) {
            int evictID = lruList.back().first;
            auto evictNode = lruList.back().second;
            lruList.pop_back();

            // Save the evicted node to LevelDB asynchronously
            std::thread([this, evictID, evictNode]() {
                saveNodeToDB(evictID, evictNode);
            }).detach();

            cache.erase(evictID);
        }

        lruList.emplace_front(nodeID, qnode);
        cache[nodeID] = lruList.begin();
    }
}

void LRUCache::invalidate(int nodeID) {
    std::lock_guard<std::mutex> lock(cacheMutex); // Lock only during cache modification
    if (lockedNodes.find(nodeID) != lockedNodes.end()) {
        std::cerr << "Warning: Attempted to evict a locked node with ID " << nodeID << std::endl;
        return;
    }

    if (cache.find(nodeID) != cache.end()) {
        auto it = cache[nodeID];
        lruList.erase(it);
        cache.erase(nodeID);
    }
}

void LRUCache::lockNode(int nodeID) {
    std::lock_guard<std::mutex> lock(cacheMutex);
    lockedNodes.insert(nodeID);
}

void LRUCache::unlockNode(int nodeID) {
    std::lock_guard<std::mutex> lock(cacheMutex);
    lockedNodes.erase(nodeID);
}
