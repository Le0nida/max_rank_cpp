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

class QNode; // Forward declaration
class Context;

class LRUCache {
public:
    explicit LRUCache(Context& ctx);
    ~LRUCache();

    std::shared_ptr<QNode> get(int nodeID);
    void add(std::shared_ptr<QNode> qnode);
    void invalidate(int nodeID);

    void lockNode(int nodeID);
    void unlockNode(int nodeID);

    void cleanup();
    void clear();

private:
    Context& ctx;

    struct CacheNode {
        int nodeID;
        std::shared_ptr<QNode> qnode;
        int frequency;
    };

    std::map<int, std::list<CacheNode>> freqMap;
    std::unordered_map<int, std::list<CacheNode>::iterator> nodeMap;

    std::unordered_map<int, std::list<std::pair<int, std::shared_ptr<QNode>>>::iterator> cache;
    std::list<std::pair<int, std::shared_ptr<QNode>>> lruList;

    int cacheSize;

    std::unordered_set<int> lockedNodes;

    std::mutex ioMutex;
    std::condition_variable ioCondVar;
    std::unordered_map<int, std::future<void>> pendingOperations;

    struct IndexEntry {
        size_t offset;
    };
    std::unordered_map<int, IndexEntry> index;

    std::string dataFilePath = "nodes.dat";

    HANDLE hFile = INVALID_HANDLE_VALUE;
    HANDLE hMapFile = NULL;
    size_t fileSize = 0;
    char* pFileData = nullptr;
    size_t currentOffset = 0;

    void openMemoryMappedFile();
    void closeMemoryMappedFile();
};

#endif // LRUCACHE_H

