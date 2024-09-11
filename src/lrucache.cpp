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

std::shared_ptr<QNode> LRUCache::get(int nodeID) {
    if (cache.find(nodeID) != cache.end()) {
        lruList.splice(lruList.begin(), lruList, cache[nodeID]);
        return cache[nodeID]->second;
    }

    // Load from disk if not in cache
    auto node = std::make_shared<QNode>();
    std::string filePath = getFilePath(nodeID);
    node->loadFromDisk(filePath);
    add(node);
    return node;
}

void LRUCache::lockNode(int nodeID) {
    // Aggiungi l'ID del nodo all'insieme dei nodi bloccati
    lockedNodes.insert(nodeID);
}

void LRUCache::unlockNode(int nodeID) {
    // Rimuovi l'ID del nodo dall'insieme dei nodi bloccati
    lockedNodes.erase(nodeID);
}

void LRUCache::add(std::shared_ptr<QNode> qnode) {

    int nodeID = qnode->getNodeID();
    invalidate(nodeID);
    if (cache.size() >= cacheSize) {
        // Rimuovi il nodo meno recentemente usato
        int evictID = lruList.back().first;  // Ottieni l'ID del nodo da rimuovere
        auto evictNode = lruList.back().second;  // Ottieni il nodo da rimuovere
        lruList.pop_back();

        // Salva il nodo su disco prima di rimuoverlo dalla cache
        evictNode->saveToDisk(getFilePath(evictID));
        cache.erase(evictID);  // Rimuovi il nodo dalla cache
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

std::string LRUCache::getFilePath(int nodeID) {
    return "node_" + std::to_string(nodeID) + ".dat";
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