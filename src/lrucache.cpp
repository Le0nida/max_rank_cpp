#include "lrucache.h"
#include <fstream>
#include <iostream>
#include <filesystem>
#include <future>

// Definition of the global cache
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
        // Estimate the size of a QNode in memory (this is a rough estimate)
        size_t estimatedNodeSize = sizeof(QNode) + 1024; // Adjust as necessary
        cacheSize = static_cast<int>(cacheMemoryUsage / estimatedNodeSize);

        if (cacheSize < MIN_CACHE_SIZE) {
            cacheSize = MIN_CACHE_SIZE;
        } else if (cacheSize > MAX_CACHE_SIZE) {
            cacheSize = MAX_CACHE_SIZE;
        }
        std::cout << "Cache size set to " << cacheSize << " based on available memory." << std::endl;
    }

    // Open the memory-mapped file
    openMemoryMappedFile();
}

LRUCache::~LRUCache() {
    // Ensure all nodes in cache are saved
    for (auto& [nodeID, nodeIter] : cache) {
        auto& node = nodeIter->second;

        // Save the node to the memory-mapped file
        size_t offset = currentOffset;
        size_t dataSize = node->saveToDisk(pFileData, offset);

        // Update currentOffset
        currentOffset += dataSize;

        // Update the index
        index[nodeID] = {offset};
    }

    // Close the memory-mapped file
    closeMemoryMappedFile();

    // Perform cleanup if necessary
    cleanup();
}

void LRUCache::openMemoryMappedFile() {
    const size_t MAX_FILE_SIZE = 6ULL * 1024ULL * 1024ULL * 1024ULL; // 6 GB

    // Open or create the file
    hFile = CreateFileA(dataFilePath.c_str(), GENERIC_READ | GENERIC_WRITE, 0, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
    if (hFile == INVALID_HANDLE_VALUE) {
        std::cerr << "Failed to open or create data file: " << dataFilePath << std::endl;
        exit(EXIT_FAILURE);
    }

    // Set the file size to MAX_FILE_SIZE
    LARGE_INTEGER liFileSize;
    liFileSize.QuadPart = MAX_FILE_SIZE;
    if (!SetFilePointerEx(hFile, liFileSize, NULL, FILE_BEGIN) || !SetEndOfFile(hFile)) {
        std::cerr << "Failed to set file size." << std::endl;
        CloseHandle(hFile);
        exit(EXIT_FAILURE);
    }
    fileSize = MAX_FILE_SIZE;

    // Create a file mapping object
    hMapFile = CreateFileMapping(hFile, NULL, PAGE_READWRITE, liFileSize.HighPart, liFileSize.LowPart, NULL);
    if (hMapFile == NULL) {
        std::cerr << "Failed to create file mapping object." << std::endl;
        CloseHandle(hFile);
        exit(EXIT_FAILURE);
    }

    // Map the entire file into the address space
    pFileData = static_cast<char*>(MapViewOfFile(hMapFile, FILE_MAP_ALL_ACCESS, 0, 0, 0));
    if (pFileData == NULL) {
        std::cerr << "Failed to map view of file." << std::endl;
        CloseHandle(hMapFile);
        CloseHandle(hFile);
        exit(EXIT_FAILURE);
    }

    // Initialize currentOffset to 0
    currentOffset = 0;
}

void LRUCache::closeMemoryMappedFile() {
    if (pFileData) {
        UnmapViewOfFile(pFileData);
        pFileData = nullptr;
    }
    if (hMapFile != NULL) {
        CloseHandle(hMapFile);
        hMapFile = NULL;
    }
    if (hFile != INVALID_HANDLE_VALUE) {
        CloseHandle(hFile);
        hFile = INVALID_HANDLE_VALUE;
    }
}

std::shared_ptr<QNode> LRUCache::get(int nodeID) {
    // Check if node is in cache
    if (cache.find(nodeID) != cache.end()) {
        lruList.splice(lruList.begin(), lruList, cache[nodeID]);
        return cache[nodeID]->second;
    }

    // Not in cache, load from disk
    auto it = index.find(nodeID);
    if (it != index.end()) {
        auto node = std::make_shared<QNode>();

        size_t offset = it->second.offset;

        // Load the node from memory-mapped file
        node->loadFromDisk(pFileData, offset);

        add(node);
        return node;
    } else {
        std::cerr << "Node with ID " << nodeID << " not found in index." << std::endl;
        return nullptr;
    }
}

void LRUCache::add(std::shared_ptr<QNode> qnode) {
    int nodeID = qnode->getNodeID();
    invalidate(nodeID);

    if (cache.size() >= cacheSize) {
        // Remove the least recently used node
        int evictID = lruList.back().first;
        auto evictNode = lruList.back().second;
        lruList.pop_back();

        // Save the node to the memory-mapped file
        size_t offset = currentOffset;

        // Ensure we have enough space
        size_t estimatedSize = evictNode->estimatedSize() + sizeof(size_t); // Include size field
        if (currentOffset + estimatedSize > fileSize) {
            std::cerr << "Error: Exceeded maximum file size for memory-mapped file." << std::endl;
            // Handle error, perhaps expand the file or exit
            // For now, we'll exit
            exit(EXIT_FAILURE);
        }

        size_t dataSize = evictNode->saveToDisk(pFileData, offset);

        // Update currentOffset
        currentOffset += dataSize;

        // Update the index with the offset
        index[evictID] = {offset};

        cache.erase(evictID);
    }

    lruList.emplace_front(nodeID, qnode);
    cache[nodeID] = lruList.begin();
}

void LRUCache::invalidate(int nodeID) {
    // If the node is locked, do nothing
    if (lockedNodes.find(nodeID) != lockedNodes.end()) {
        std::cerr << "Warning: Attempted to evict a locked node with ID " << nodeID << std::endl;
        return;
    }

    // If the node is present in the cache, remove it
    if (cache.find(nodeID) != cache.end()) {
        auto it = cache[nodeID];
        lruList.erase(it);  // Remove the node from the LRU list
        cache.erase(nodeID);  // Remove the node from the cache
    }
}

void LRUCache::lockNode(int nodeID) {
    // Add the node ID to the set of locked nodes
    lockedNodes.insert(nodeID);
}

void LRUCache::unlockNode(int nodeID) {
    // Remove the node ID from the set of locked nodes
    lockedNodes.erase(nodeID);
}

// Method to delete .dat files in batch
void batchDeleteFiles(const std::vector<std::string>& filepaths) {
    for (const auto& path : filepaths) {
        std::filesystem::remove(path);
    }
}

// Utility function to retrieve all .dat files from a directory
std::vector<std::string> getAllDatFiles(const std::string& directoryPath) {
    std::vector<std::string> files;
    for (const auto& entry : std::filesystem::directory_iterator(directoryPath)) {
        if (entry.path().extension() == ".dat") {
            files.push_back(entry.path().string());
        }
    }
    return files;
}

// Improved version of cleanup with parallel deletion
void LRUCache::cleanup() {
    std::vector<std::string> files = getAllDatFiles(".");

    // Divide files into batches for parallel deletion
    const int batchSize = 6000;
    std::vector<std::future<void>> futures;

    for (size_t i = 0; i < files.size(); i += batchSize) {
        std::vector<std::string> batch(files.begin() + i,
                                       files.begin() + std::min(files.size(), i + batchSize));

        // Create a separate thread for each batch
        futures.push_back(std::async(std::launch::async, batchDeleteFiles, batch));
    }

    // Wait for all threads to finish
    for (auto& future : futures) {
        future.get();
    }
}
