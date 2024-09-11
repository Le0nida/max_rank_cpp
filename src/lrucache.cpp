//
// Created by leona on 09/09/2024.
//

#include "lrucache.h"
#include <fstream>
#include <iostream>

// Definizione della cache globale
LRUCache globalCache(10000);

std::shared_ptr<QNode> LRUCache::get(int nodeID) {
    if (cache.find(nodeID) != cache.end()) {
        lruList.splice(lruList.begin(), lruList, cache[nodeID]);
        return cache[nodeID]->second;
    }

    // Load from disk if not in cache
    auto node = std::make_shared<QNode>();
    node->loadFromDisk(getFilePath(nodeID));
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
