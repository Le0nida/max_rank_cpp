//
// Created by leona on 09/09/2024.
//

#include "lrucache.h"
#include <fstream>
#include <iostream>

QNode* LRUCache::get(int nodeID) {
    // Controlla se il nodo è già in cache
    if (cache.find(nodeID) != cache.end()) {
        // Muovi il nodo in cima alla lista LRU
        lruList.remove(nodeID);
        lruList.push_front(nodeID);
        return cache[nodeID].get();  // Restituisce il puntatore al nodo
    }

    // Se non è in cache, carica il nodo dal disco
    std::unique_ptr<QNode> node = std::make_unique<QNode>();
    node->loadFromDisk(getFilePath(nodeID));

    // Se la cache è piena, rimuovi il nodo meno recentemente usato
    if (cache.size() >= cacheSize) {
        int evictID = lruList.back();
        lruList.pop_back();
        cache[evictID]->saveToDisk(getFilePath(evictID));
        cache.erase(evictID);
    }

    // Inserisci il nuovo nodo nella cache e sposta la proprietà
    cache[nodeID] = std::move(node);
    lruList.push_front(nodeID);
    return cache[nodeID].get();  // Restituisce il puntatore al nodo
}

// Salva lo stato della cache su disco
void LRUCache::saveCacheToDisk() {
    for (auto& pair : cache) {
        int nodeID = pair.first;
        pair.second->saveToDisk(getFilePath(nodeID));
    }
}

std::string LRUCache::getFilePath(int nodeID) {
    return "node_" + std::to_string(nodeID) + ".dat";
}