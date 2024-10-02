#ifndef LRUCACHE_H
#define LRUCACHE_H

#include <future>
#include <iostream>
#include <unordered_map>
#include <list>
#include <map>
#include <memory>
#include <mutex>
#include <unordered_set>
#include <thread>
#include <windows.h> // For memory-mapped files on Windows

#include "qnode.h"
#include "utils.h"

class LRUCache {
private:
    struct CacheNode {
        int nodeID;
        std::shared_ptr<QNode> qnode;
        int frequency;
    };

    // Map frequency to a list of CacheNodes
    std::map<int, std::list<CacheNode>> freqMap;

    // Map nodeID to iterator pointing to the CacheNode in freqMap
    std::unordered_map<int, std::list<CacheNode>::iterator> nodeMap;

    // Map of nodeID to iterator in lruList
    std::unordered_map<int, std::list<std::pair<int, std::shared_ptr<QNode>>>::iterator> cache;

    // LRU list to track node usage
    std::list<std::pair<int, std::shared_ptr<QNode>>> lruList;

    // Maximum cache size
    int cacheSize;

    // Contains IDs of "locked" nodes
    std::unordered_set<int> lockedNodes;

    // For asynchronous I/O
    std::mutex ioMutex;
    std::condition_variable ioCondVar;
    std::unordered_map<int, std::future<void>> pendingOperations;

    // Index mapping nodeID to file offset
    struct IndexEntry {
        size_t offset;
    };
    std::unordered_map<int, IndexEntry> index;

    // Path to data file
    std::string dataFilePath = "nodes.dat";

    // Memory-mapped file variables
    HANDLE hFile = INVALID_HANDLE_VALUE;
    HANDLE hMapFile = NULL;
    size_t fileSize = 0;
    char* pFileData = nullptr;
    size_t currentOffset = 0;

    // Helper methods for memory-mapped file
    void openMemoryMappedFile();
    void closeMemoryMappedFile();

public:
    explicit LRUCache();
    ~LRUCache();

    // Get a node by ID, load from disk if not present in cache
    std::shared_ptr<QNode> get(int nodeID);

    void add(std::shared_ptr<QNode> qnode);
    void invalidate(int nodeID);

    void lockNode(int nodeID);     // Lock the node in cache
    void unlockNode(int nodeID);   // Unlock the node in cache

    void cleanup();
};

// Declaration of the global cache
extern LRUCache globalCache;

#endif // LRUCACHE_H
