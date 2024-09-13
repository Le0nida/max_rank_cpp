//
// Created by leona on 09/09/2024.
//

#include "lrucache.h"
#include <fstream>
#include <iostream>
#include <filesystem>
#include <future>

// Definizione della cache globale
LRUCache globalCache(10000);

// lrucache.cpp
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
        node->loadFromDisk(dataFile, it->second.offset);
        add(node);
        return node;
    } else {
        std::cerr << "Node with ID " << nodeID << " not found in index." << std::endl;
        return nullptr;
    }
}

void LRUCache::lockNode(int nodeID) {
    // Aggiungi l'ID del nodo all'insieme dei nodi bloccati
    lockedNodes.insert(nodeID);
}

void LRUCache::unlockNode(int nodeID) {
    // Rimuovi l'ID del nodo dall'insieme dei nodi bloccati
    lockedNodes.erase(nodeID);
}

// lrucache.cpp
void LRUCache::add(std::shared_ptr<QNode> qnode) {
    int nodeID = qnode->getNodeID();
    invalidate(nodeID);

    if (cache.size() >= cacheSize) {
        // Remove the least recently used node
        int evictID = lruList.back().first;
        auto evictNode = lruList.back().second;
        lruList.pop_back();

        // Ensure we're at the end of the file
        dataFile.clear(); // Clear any error flags
        dataFile.seekp(0, std::ios::end);

        // Get the current offset
        std::streampos offset = dataFile.tellp();

        // Save the node to the data file
        evictNode->saveToDisk(dataFile);

        // Flush the output to synchronize the stream
        dataFile.flush();

        // Update the index with the offset only
        index[evictID] = {offset};

        cache.erase(evictID);
    }

    lruList.emplace_front(nodeID, qnode);
    cache[nodeID] = lruList.begin();
}


void LRUCache::invalidate(int nodeID) {
    // Se il nodo è bloccato, non fare nulla
    if (lockedNodes.find(nodeID) != lockedNodes.end()) {
        std::cerr << "Warning: Attempted to evict a locked node with ID " << nodeID << std::endl;
        return;
    }

    // Se il nodo è presente nella cache, rimuovilo
    if (cache.find(nodeID) != cache.end()) {
        auto it = cache[nodeID];
        lruList.erase(it);  // Rimuovi il nodo dalla lista LRU
        cache.erase(nodeID);  // Rimuovi il nodo dalla cache
    }
}

// Metodo per cancellare i file .dat in batch
void batchDeleteFiles(const std::vector<std::string>& filepaths) {
    for (const auto& path : filepaths) {
        std::filesystem::remove(path);
    }
}

// Funzione di utility per recuperare tutti i file .dat da una directory
std::vector<std::string> getAllDatFiles(const std::string& directoryPath) {
    std::vector<std::string> files;
    for (const auto& entry : std::filesystem::directory_iterator(directoryPath)) {
        if (entry.path().extension() == ".dat") {
            files.push_back(entry.path().string());
        }
    }
    return files;
}

// Versione migliorata di cleanup, con cancellazione parallela
void LRUCache::cleanup() {
    std::vector<std::string> files = getAllDatFiles(".");

    // Dividi i file in batch per cancellazioni parallele
    const int batchSize = 6000;
    std::vector<std::future<void>> futures;

    for (size_t i = 0; i < files.size(); i += batchSize) {
        std::vector<std::string> batch(files.begin() + i,
                                       files.begin() + std::min(files.size(), i + batchSize));

        // Crea un thread separato per ogni batch
        futures.push_back(std::async(std::launch::async, batchDeleteFiles, batch));
    }

    // Attendi che tutti i thread abbiano terminato
    for (auto& future : futures) {
        future.get();
    }
}

void LRUCache::loadIndex() {
    std::ifstream in(indexFilePath, std::ios::binary);
    if (!in.is_open()) {
        // Index file doesn't exist; start with empty index
        return;
    }

    size_t indexSize;
    in.read(reinterpret_cast<char*>(&indexSize), sizeof(indexSize));

    for (size_t i = 0; i < indexSize; ++i) {
        int nodeID;
        IndexEntry entry;
        in.read(reinterpret_cast<char*>(&nodeID), sizeof(nodeID));
        in.read(reinterpret_cast<char*>(&entry.offset), sizeof(entry.offset));
        index[nodeID] = entry;
    }

    in.close();
}

void LRUCache::saveIndex() {
    std::ofstream out(indexFilePath, std::ios::binary | std::ios::trunc);
    if (!out.is_open()) {
        std::cerr << "Failed to open index file for writing: " << indexFilePath << std::endl;
        return;
    }

    size_t indexSize = index.size();
    out.write(reinterpret_cast<const char*>(&indexSize), sizeof(indexSize));

    for (const auto& [nodeID, entry] : index) {
        out.write(reinterpret_cast<const char*>(&nodeID), sizeof(nodeID));
        out.write(reinterpret_cast<const char*>(&entry.offset), sizeof(entry.offset));
    }

    out.close();
}