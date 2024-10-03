//
// Created by leona on 03/10/2024.
//

#include "context.h"
#include "lrucache.h"
#include "halfspace.h"

Context::Context(size_t cacheSize)
    : halfspaceCache(std::make_shared<HalfSpaceCache>(cacheSize)),
      cache(std::make_shared<LRUCache>(*this)),
      globalNodeID(0) {}