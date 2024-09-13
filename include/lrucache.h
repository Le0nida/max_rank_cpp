//
// Created by leona on 09/09/2024.
//

#ifndef LRUCACHE_H
#define LRUCACHE_H

#include <future>
#include <unordered_map>
#include <list>
#include <memory>  // Include per std::unique_ptr
#include <mutex>
#include <unordered_set>
#include <thread>

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

    // For asynchronous I/O
    std::mutex ioMutex;
    std::condition_variable ioCondVar;
    std::unordered_map<int, std::future<void>> pendingOperations;

    // Index mapping nodeID to file offset and size
    struct IndexEntry {
        std::streampos offset;
    };
    std::unordered_map<int, IndexEntry> index;

    // File stream for nodes.dat
    std::fstream dataFile;

    // Path to data and index files
    std::string dataFilePath = "nodes.dat";
    std::string indexFilePath = "index.dat";

    // Helper methods
    void loadIndex();
    void saveIndex();

public:
    explicit LRUCache(int size) : cacheSize(size) {
        // Open the data file in read-write mode, create if doesn't exist
        dataFile.open(dataFilePath, std::ios::in | std::ios::out | std::ios::binary);
        if (!dataFile.is_open()) {
            // File doesn't exist, create it
            dataFile.open(dataFilePath, std::ios::out | std::ios::binary);
            dataFile.close();
            dataFile.open(dataFilePath, std::ios::in | std::ios::out | std::ios::binary);
        }

        // Load the index from disk
        loadIndex();
    }

    ~LRUCache() {
        // Ensure all nodes in cache are saved
        for (auto& [nodeID, nodeIter] : cache) {
            auto& node = nodeIter->second;
            std::streampos offset = dataFile.tellp();
            size_t size = node->saveToDisk(dataFile);
            index[nodeID] = {offset};
        }

        // Save index to disk
        saveIndex();

        // Close data file
        if (dataFile.is_open()) {
            dataFile.close();
        }

        // Perform cleanup if necessary
        cleanup();
    }

    // Get a node by ID, load from disk if not present in cache
    std::shared_ptr<QNode> get(int nodeID);

    void add(std::shared_ptr<QNode> qnode);
    void invalidate(int nodeID);

    void lockNode(int nodeID);     // Blocca il nodo in cache
    void unlockNode(int nodeID);   // Sblocca il nodo in cache

    void cleanup();


};

// Dichiarazione della cache globale
extern LRUCache globalCache;

#endif // LRUCACHE_H
