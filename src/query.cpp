#include "query.h"
#include <algorithm>
#include <cmath>

// ------------------------------------------------------------------------
// Existing helper functions (unchanged)
// ------------------------------------------------------------------------
std::vector<Point> getdominators(const std::vector<Point>& data, const Point& p)
{
    std::vector<Point> dominators;
    dominators.reserve(data.size() / 2);
    for (const auto& r : data) {
        bool less_equal    = true;
        bool strictly_less = false;
        for (int i = 0; i < p.dims; ++i) {
            if (r.coord[i] > p.coord[i]) {
                less_equal = false;
                break;
            }
            if (r.coord[i] < p.coord[i]) {
                strictly_less = true;
            }
        }
        if (less_equal && strictly_less) {
            dominators.push_back(r);
        }
    }
    return dominators;
}

std::vector<Point> getdominees(const std::vector<Point>& data, const Point& p)
{
    std::vector<Point> dominees;
    dominees.reserve(data.size() / 2);
    for (const auto& r : data) {
        bool greater_equal = true;
        bool strictly_greater = false;
        for (int i = 0; i < p.dims; ++i) {
            if (r.coord[i] < p.coord[i]) {
                greater_equal = false;
                break;
            }
            if (r.coord[i] > p.coord[i]) {
                strictly_greater = true;
            }
        }
        if (greater_equal && strictly_greater) {
            dominees.push_back(r);
        }
    }
    return dominees;
}

std::vector<Point> getincomparables(const std::vector<Point>& data, const Point& p)
{
    std::vector<Point> incomp;
    incomp.reserve(data.size() / 2);
    for (const auto& r : data) {
        bool less = false;
        bool greater = false;
        for (int i = 0; i < p.dims; ++i) {
            if (r.coord[i] < p.coord[i]) less = true;
            if (r.coord[i] > p.coord[i]) greater = true;
        }
        if (less && greater) {
            incomp.push_back(r);
        }
    }
    return incomp;
}

// ------------------------------------------------------------------------
// Production-ready getskyline using libspatialindex (official API)
// with bulk loading + incremental dominance
// ------------------------------------------------------------------------
#include <spatialindex/SpatialIndex.h>
#include <memory>
#include <iostream>
#include <vector>
#include <limits>


/**
 * \class SkylineDataStream
 * \brief Feeds data to R*-tree bulk-load. Must implement Tools::IObjectStream.
 */
class SkylineDataStream : public SpatialIndex::IDataStream
                       , public Tools::IObjectStream // some versions require this
{
public:
    SkylineDataStream(const std::vector<Point>& pts)
        : m_points(pts), m_dims(pts.empty()?0:pts[0].dims),
          m_index(0), m_size((uint32_t)pts.size())
    {}

    // ========== Methods from IDataStream ==========

    bool hasNext() /*no override*/ {
        return (m_index < m_points.size());
    }

    SpatialIndex::IData* getNext() /*no override*/ {
        if (!hasNext()) return nullptr;
        const auto &p = m_points[m_index++];
        double* low  = new double[m_dims];
        double* high = new double[m_dims];
        for (int i=0; i<m_dims; i++){
            low[i]  = p.coord[i];
            high[i] = p.coord[i];
        }
        SpatialIndex::Region r(low, high, m_dims);

        auto id = static_cast<SpatialIndex::id_type>(p.id);

        return new SpatialIndex::RTree::Data(0, nullptr, r, id);
    }

    double getCurrentRadius() /*no override*/ {
        // Some older versions require this;
        // if your version has no getCurrentRadius,
        // just remove it.
        return 0.0;
    }

    // ========== Methods from Tools::IObjectStream ==========
    // Some versions are pure virtual:
    uint32_t size() /*no override*/ {
        return m_size;
    }

    void rewind() /*no override*/ {
        m_index=0;
    }

    // Possibly Tools::IObject::getIdentifier() is also pure virtual
    // but typically not. If so, define a dummy.

private:
    const std::vector<Point>& m_points;
    int m_dims;
    size_t m_index;
    uint32_t m_size;
};

/**
 * \class SkylineVisitor
 * \brief A visitor to retrieve all data from the tree in a single pass.
 *        Some versions require visitData(std::vector<const IData*>&).
 */
class SkylineVisitor : public SpatialIndex::IVisitor
{
public:
    std::vector<Point> retrieved;

    void visitData(std::vector<const SpatialIndex::IData*>& v) override /* no override */ {
        for (auto* ptr : v) {
            if (ptr) visitData(*ptr);
        }
    }

    // Single-data call
    void visitData(const SpatialIndex::IData& d) override /*no override*/ {
        SpatialIndex::IShape* shape = nullptr;
        d.getShape(&shape);
        if (!shape) return;

        auto region = dynamic_cast<SpatialIndex::Region*>(shape);
        if (region) {
            int dims = (int) region->getDimension();
            std::vector<double> c(dims);
            for (int i=0; i<dims; i++){
                c[i] = region->getLow(i);
            }
            retrieved.emplace_back(c, (int)d.getIdentifier());
        }
        delete shape;
    }

    // Some versions also have a pure virtual
    // "visitData(const std::vector<const IData*>& v)"
    // We define it as well, just calling visitData(...) in a loop.
    void visitData(const std::vector<const SpatialIndex::IData*>& v) /*no override*/ {
        for (auto &ptr : v) {
            if (ptr) visitData(*ptr);
        }
    }

    // Optionally visitNode
    void visitNode(const SpatialIndex::INode& n) /*no override*/ {
        // We do nothing special
    }
};

/**
 * \brief Checks if s strictly dominates p in all dims, s != p.
 */
static bool sDominates(const Point& p, const Point& s)
{
    bool strictly=false;
    for (int i=0; i<p.dims; i++){
        if (s.coord[i] > p.coord[i]) return false;
        if (s.coord[i] < p.coord[i]) strictly=true;
    }
    return strictly;
}

/**
 * \brief An incremental approach to build the final skyline, like your naive method but on a smaller set.
 */
static std::vector<Point> incrementalSkyline(std::vector<Point> data)
{
    std::vector<Point> window;
    window.reserve(data.size()/2);
    auto dominates = [&](const Point& p, const Point& r){
        bool less_equal=true, strictly_less=false;
        for (int i=0; i<p.dims; i++){
            if (p.coord[i] > r.coord[i]){
                less_equal=false;
                break;
            }
            if (p.coord[i] < r.coord[i]){
                strictly_less=true;
            }
        }
        return (less_equal && strictly_less);
    };

    for (auto &pnt : data){
        bool isDom=false;
        for (auto &w: window){
            if (dominates(w, pnt)) {
                isDom=true;
                break;
            }
        }
        if (!isDom) {
            window.erase(std::remove_if(window.begin(), window.end(),
                            [&](const Point &w2){return dominates(pnt, w2);}),
                         window.end());
            window.push_back(pnt);
        }
    }
    return window;
}

/**
 * \brief Final getskyline:
 *        1) Bulk-load the dataset into an R*-tree
 *        2) Single "infinite region" query to retrieve all data
 *        3) final incremental pass
 */
std::vector<Point> getskyline(const std::vector<Point>& data)
{
    if (data.empty()) {
        return {};
    }

    // 1) Bulk load
    std::unique_ptr<SpatialIndex::IStorageManager> store(
        SpatialIndex::StorageManager::createNewMemoryStorageManager()
    );
    SkylineDataStream ds(data);

    int dims = data[0].dims;
    SpatialIndex::id_type indexID=1;
    double fillFactor=0.7;
    uint32_t indexCap=100, leafCap=100;

    // Some versions do: SpatialIndex::RTree::createAndBulkLoadNewRTree
    // Others do: SpatialIndex::RTree::createAndBulkLoadNewRTree(...
    // The function is the same in official 1.8, just ensure the arguments match.
    std::unique_ptr<SpatialIndex::ISpatialIndex> tree(
        SpatialIndex::RTree::createAndBulkLoadNewRTree(
            SpatialIndex::RTree::BLM_STR,
            ds,
            *store,
            fillFactor,
            indexCap,
            leafCap,
            (uint32_t)dims,
            SpatialIndex::RTree::RV_RSTAR,
            indexID
        )
    );

    // 2) Single pass retrieval => infinite region
    SkylineVisitor visitor;
    std::vector<double> low(dims, -std::numeric_limits<double>::max());
    std::vector<double> high(dims, std::numeric_limits<double>::max());
    SpatialIndex::Region fullRegion(&low[0], &high[0], (uint32_t)dims);

    // intersectsWithQuery => calls visitor.visitData(...) for each record
    tree->intersectsWithQuery(fullRegion, visitor);

    // 3) final incremental pass
    auto sky = incrementalSkyline(std::move(visitor.retrieved));
    return sky;
}
