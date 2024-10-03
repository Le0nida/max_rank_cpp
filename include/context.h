//
// Created by leona on 03/10/2024.
//

#ifndef CONTEXT_H
#define CONTEXT_H

#include "geom.h" // For Point and PointHash
#include <unordered_map>
#include <memory>

class HalfSpaceCache; // Forward declaration
class LRUCache;       // Forward declaration

class Context {
public:
    // Member variables
    std::shared_ptr<HalfSpaceCache> halfspaceCache;
    std::unordered_map<Point, long, PointHash> pointToHalfSpaceCache;
    std::shared_ptr<LRUCache> cache;
    int globalNodeID;

    // Constructor declaration
    explicit Context(size_t cacheSize);
};

#endif // CONTEXT_H