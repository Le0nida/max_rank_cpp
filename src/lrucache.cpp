//
// Created by leona on 09/09/2024.
//

#include "lrucache.h"
#include <fstream>
#include <iostream>

// Definizione della cache globale
LRUCache globalCache(100);

std::shared_ptr<QNode> LRUCache::get(int nodeID) {
    if (cache.find(nodeID) != cache.end()) {
        lruList.remove(nodeID);
        lruList.push_front(nodeID);
        return cache[nodeID];
    }

    // Load from disk if not in cache
    auto node = std::make_shared<QNode>();
    node->loadFromDisk(getFilePath(nodeID));
    add(node);
    return cache[nodeID];
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
        int evictID = lruList.back();
        lruList.pop_back();
        cache[evictID]->saveToDisk(getFilePath(evictID));  // Salva su disco
        cache.erase(evictID);
    }
    cache[nodeID] = std::move(qnode);
    lruList.push_front(nodeID);
}

void LRUCache::invalidate(int nodeID) {
    if (lockedNodes.find(nodeID) != lockedNodes.end()) {
        // Se il nodo è bloccato, non fare nulla
        std::cerr << "Warning: Attempted to evict a locked node with ID " << nodeID << std::endl;
        return;
    }
    if (cache.find(nodeID) != cache.end()) {
        cache.erase(nodeID);
        lruList.remove(nodeID);
    }
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
