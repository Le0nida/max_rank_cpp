#include "maxrank.h"

int numOfSubdivisions = 0;

// Definisci e inizializza le variabili globali
HalfSpaceCache* halfspaceCache = nullptr;
std::unordered_map<Point, long, PointHash> pointToHalfSpaceCache;

std::pair<int, std::vector<Cell>> aa_hd(const std::vector<Point>& data, const Point& p) {

    // Reset global variables
    halfspaceCache = nullptr;
    pointToHalfSpaceCache.clear();

    int dims = static_cast<int>(p.dims - 1);
    numOfSubdivisions = (int) pow(2.0, dims);
    float queryPlane[10];
    // Imposta queryPlane in modo fisso a 1
    for (int i = 0; i < dims + 1; i++) {
        queryPlane[i] = 1.0f;
    }
    std::vector<std::string> Comb = readCombinations(dims);
    for (int i = 0; i < dims + 1; i++) queryPlane[i] = 1;

    QTree qt(dims, maxCapacityQNode, maxLevelQTree);
    std::vector<Point> dominators = getdominators(data, p);
    std::vector<Point> incomp = getincomparables(data, p);

    // Inizializzo la cache per gli halfspaces
    initializeCache(data.size());

    auto updateqt = [&](const std::vector<Point>& old_sky) {
        std::cout << "> getting skyline ... " << std::endl;
        auto start = std::chrono::high_resolution_clock::now();
        std::vector<Point> new_sky = getskyline(incomp);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout << "> skyline time: " << elapsed.count() << " seconds.\n" << std::endl;


        start = std::chrono::high_resolution_clock::now();
        std::cout << "> building halfspaces ... " << std::endl;
        std::vector<long> new_halfspaces = genhalfspaces(p, new_sky);
        std::vector<long> unique_new_halfspaces;
        for (const auto& hs : new_halfspaces) {
            bool found = false;
            for (const auto& os : old_sky) {
                if (os.id == hs) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                unique_new_halfspaces.push_back(hs);
            }
        }

        new_halfspaces = std::move(unique_new_halfspaces);
        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        std::cout << "> building halfspaces time: " << elapsed.count() << " seconds.\n" << std::endl;

        start = std::chrono::high_resolution_clock::now();
        std::cout << "> " << new_halfspaces.size() << " halfspace(s) to insert" << std::endl;
        if (!new_halfspaces.empty()) {
            //qt.inserthalfspaces(new_halfspaces);
            qt.inserthalfspacesMacroSplit(new_halfspaces);
            std::cout << "> " << new_halfspaces.size() << " halfspace(s) have been inserted" << std::endl;
        }
        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        std::cout << "> inserting halfspaces time: " << elapsed.count() << " seconds.\n" << std::endl;

        auto new_leaves = qt.getAllLeaves();//qt.getLeaves();
        //std::cout << "> " << new_leaves.size() << " total leaves" << std::endl;
        qt.updateAllOrders();

        std::sort(new_leaves.begin(), new_leaves.end(), [](QNode* a, QNode* b) { return a->order < b->order; });

        return std::make_pair(new_sky, new_leaves);
    };


    auto [sky, leaves] = updateqt({});

    int minorder_singular = std::numeric_limits<int>::max();
    std::vector<Cell> mincells_singular;
    int n_exp = 0;

    while (true) {
        std::cout << "Cycle number " << n_exp << std::endl;
        int minorder = std::numeric_limits<int>::max();
        std::vector<Cell> mincells;

        auto start = std::chrono::high_resolution_clock::now();
        for (auto leaf : leaves) {
            int leaf_order = static_cast<int>(leaf->order);
            if (leaf_order > minorder || leaf_order > minorder_singular) {
                break;
            }
            //prune away leaf nodes that lie about hyperplane q_1+q2+...+q_d < 1;
            if (!MbrIsValid(leaf->mbr, Comb, dims, queryPlane)) {
                continue;
            }

            int hamweight = 0;
            while (hamweight <= leaf->halfspaces.size() && leaf_order + hamweight <= minorder && leaf_order + hamweight <= minorder_singular && hamweight <= limitHamWeight) {
                // std::cout << "Hamweight " << hamweight << ", numero hs: " << leaf->halfspaces.size();
                std::vector<std::string> hamstrings = genhammingstrings(static_cast<int>(leaf->halfspaces.size()), hamweight);
                // std::cout << ", Hamstring " << hamstrings.size();
                std::vector<Cell> cells = searchmincells_lp(*leaf, hamstrings);
                // std::cout << ", Celle " << cells.size() << std::endl;
                if (!cells.empty()) {
                    for (auto& cell : cells) {
                        cell.order = leaf_order + hamweight;
                    }

                    if (minorder > leaf_order + hamweight) {
                        minorder = leaf_order + hamweight;
                        mincells = std::move(cells);
                    } else {
                        mincells.insert(mincells.end(), cells.begin(), cells.end());
                    }
                    break;
                }
                if (hamstrings.size() > maxNoBinStringToCheck) break;
                hamweight++;
            }
        }
        std::cout << "> Expansion " << n_exp << ": Found " << mincells.size() << " mincell(s)" << std::endl;
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout << "> expansion time: " << elapsed.count() << " seconds.\n" << std::endl;

        int new_singulars = 0;
        std::vector<std::shared_ptr<HalfSpace>> to_expand;
        for (auto& cell : mincells) {
            if (cell.issingular()) {
                minorder_singular = cell.order;
                mincells_singular.push_back(cell);
                new_singulars++;
            } else {
                for (const auto k : cell.covered) {
                    const auto hs = halfspaceCache->get(k);
                    if (hs->arr == Arrangement::AUGMENTED && std::find(to_expand.begin(), to_expand.end(), hs) == to_expand.end()) {
                        to_expand.push_back(hs);
                    }
                }
            }
        }
        if (new_singulars > 0) {
            std::cout << "> Expansion " << n_exp << ": Found " << new_singulars << " singular mincell(s) with a minorder of " << minorder_singular << std::endl;
        }

        if (to_expand.empty()) {
            return {static_cast<int>(dominators.size()) + minorder_singular + 1, mincells_singular};
        }

        n_exp++;
        std::cout << "> Expansion " << n_exp << ": " << to_expand.size() << " halfspace(s) will be expanded" << std::endl;
        for (const auto& hs : to_expand) {
            hs->arr = Arrangement::SINGULAR;
            auto it = std::find_if(incomp.begin(), incomp.end(), [&](const Point& pt) { return hs->pntID == pt.id; });
            if (it != incomp.end()) {
                incomp.erase(it);
            }

        }
        std::tie(sky, leaves) = updateqt(sky);
    }
}


std::pair<int, std::vector<Interval>> aa_2d(const std::vector<Point>& data, const Point& p) {
    // 1) Troviamo dominatori e incomparabili
    std::vector<Point> dominators = getdominators(data, p);
    std::vector<Point> incomp     = getincomparables(data, p);

    // 2) Skyline dei soli incomparabili
    std::vector<Point> sky = getskyline(incomp);

    // 3) Creiamo la halfline relativa al punto p
    //    (in 2D la score function di p si rappresenta come retta y = m*x + q)
    auto p_line = std::make_shared<HalfLine>(p);

    // 4) Costruiamo un vettore di Intervals che rappresentano le "suddivisioni"
    //    dell'asse x in cui l'ordine non cambia
    std::vector<Interval> intervals;
    intervals.reserve(sky.size() + 1);

    // Per ogni punto in sky:
    // - Creiamo la halfline corrispondente
    // - Troviamo l'intersezione con la halfline di p
    // - Creiamo un intervallo [NaN, x_intersect], se x_intersect esiste,
    //   con coversleft = (sp_line->q < p_line->q) (come nel Python).
    for (const auto& sp : sky) {
        auto sp_line = std::make_shared<HalfLine>(sp);
        Point isect = find_halflines_intersection(*p_line, *sp_line);

        // Se r.m == s.m, in Python ottenevamo None: qui isect avrà Infinity
        // Gestiamo quindi il caso parallelo
        if (std::isinf(isect.coord[0])) {
            // Niente intersezione "finita": a tua scelta se ignorare o gestire diversamente
            continue;
        }

        double x_intersect = isect.coord[0];
        bool coversleft    = (sp_line->q < p_line->q);

        intervals.emplace_back(
            sp_line,
            std::make_pair(std::numeric_limits<double>::quiet_NaN(), x_intersect),
            coversleft
        );
    }

    // 5) Aggiungiamo anche un intervallo "fittizio" che arriva fino a x=1
    //    (come fa la versione Python, Interval(None, [NaN, 1], false))
    intervals.emplace_back(
        nullptr,
        std::make_pair(std::numeric_limits<double>::quiet_NaN(), 1.0),
        false
    );

    std::cout << "> " << sky.size() << " halfline(s) have been inserted" << std::endl;

    // 6) Avviamo il ciclo di espansione
    int n_exp = 0;
    std::vector<Interval> mincells_singular;  // Per collezionare i mincells singolari finali

    while (true) {
        // 6a) Ordiniamo gli intervalli in base all'estremo destro (range.second),
        //     esattamente come in Python sorted(key=lambda cell: cell.range[1])
        std::sort(intervals.begin(), intervals.end(), [](const Interval& a, const Interval& b){
            return a.range.second < b.range.second;
        });

        double last_end = 0.0;
        int minorder = std::numeric_limits<int>::max();
        std::vector<Interval> mincells;

        // Costruiamo la "copertura" iniziale: in Python si faceva:
        //   covering = [cell.halfline for cell in intervals if cell.coversleft]
        //   ma lì veniva gestito in modo un po’ diverso.
        //   Qui lo facciamo on-the-fly nel loop.
        std::vector<std::shared_ptr<HalfLine>> covering;
        covering.reserve(intervals.size());

        // Prima passata: raccogliamo le halflines che hanno coversleft==true
        // e che hanno range.second < 0?
        // Oppure seguiamo pedissequamente la logica Python:
        for (auto &cell : intervals) {
            if (cell.coversleft) {
                // L’idea in Python era: se coversleft è True, la halfline 'esce'
                // al confine di cell.range.second. Quindi prima di quell’ x,
                // la halfline è "attiva".
                // Se volessimo avere lo stesso effetto, potremmo aggiungere
                // la halfline a covering.
                covering.push_back(cell.halfline);
            }
        }

        // 6b) Cicliamo sugli intervalli nell’ordine (già sortato):
        for (auto &cell : intervals) {
            // L'order è semplicemente la dimensione della covering
            cell.order   = (int) covering.size();
            cell.covered = covering;

            // Aggiorniamo l’estremo sinistro dell’intervallo
            cell.range.first = last_end;
            // Aggiorniamo last_end
            last_end = cell.range.second;

            // Aggiorniamo minorder e mincells
            if (cell.order < minorder) {
                minorder = cell.order;
                mincells.clear();
                mincells.push_back(cell);
            } else if (cell.order == minorder) {
                mincells.push_back(cell);
            }

            // Ora la logica "se coversleft, pop(0), else push_back"
            if (cell.coversleft) {
                // In Python: covering.pop(0)
                // Qui potremmo fare covering.erase(covering.begin()), se non vuota
                if (!covering.empty()) {
                    covering.erase(covering.begin());
                }
            } else {
                // covering.append(cell.halfline)
                covering.push_back(cell.halfline);
            }
        }

        std::cout << "> Expansion " << n_exp << ": Found " << mincells.size() << " mincell(s)" << std::endl;

        // 6c) Controlliamo mincells per eventuali singolari
        //     e costruiamo la lista di halflines da "espandere"
        int new_singulars = 0;
        std::vector<std::shared_ptr<HalfLine>> to_expand;
        to_expand.reserve(mincells.size());
        for (auto &mc : mincells) {
            if (mc.issingular()) {
                mincells_singular.push_back(mc);
                new_singulars++;
            } else {
                // Se l'intervallo non è singolare,
                // allora le halflines in mc.covered che sono ancora AUGMENTED
                // sono candidate per l’espansione (-> SINGULAR).
                for (auto &hl : mc.covered) {
                    if (hl->arr == Arrangement::AUGMENTED) {
                        // Evitare duplicati
                        if (std::find(to_expand.begin(), to_expand.end(), hl) == to_expand.end()) {
                            to_expand.push_back(hl);
                        }
                    }
                }
            }
        }

        if (new_singulars > 0) {
            std::cout << "> Expansion " << n_exp << ": Found "
                      << new_singulars << " singular mincell(s) with a minorder of "
                      << minorder << std::endl;
        }

        // 6d) Se non ci sono halflines da espandere, abbiamo finito:
        if (to_expand.empty()) {
            // Il MaxRank in 2D è dominators.size() + minorder + 1
            return std::make_pair(
                (int)dominators.size() + minorder + 1,
                mincells_singular
            );
        }

        // Altrimenti si continua l'espansione
        n_exp++;
        std::cout << "> Expansion " << n_exp << ": "
                  << to_expand.size() << " halfline(s) will be expanded" << std::endl;

        // 6e) Segniamo come SINGULAR le halflines in to_expand
        //     e rimuoviamo i corrispondenti punti “incomparabili”:
        for (auto &hl : to_expand) {
            hl->arr = Arrangement::SINGULAR;
            // Rimuovi il punto associato a questa halfline (hl->pnt) da incomp
            auto it = std::find_if(incomp.begin(), incomp.end(), [&](const Point &pt){
                return pt.id == hl->pnt.id;
            });
            if (it != incomp.end()) {
                incomp.erase(it);
            }
        }

        // 6f) Ricostruiamo lo skyline con i rimanenti incomparabili
        std::vector<Point> new_sky = getskyline(incomp);

        // Aggiungiamo le *nuove* halflines (rispetto allo sky precedente)
        std::vector<Point> to_insert;
        for (auto &sp : new_sky) {
            if (std::find_if(sky.begin(), sky.end(), [&](const Point &q){ return q.id == sp.id; }) == sky.end()) {
                to_insert.push_back(sp);
            }
        }
        sky = new_sky; // aggiorniamo lo sky globale

        // Per ogni punto da inserire, calcoliamo l’intersezione e creiamo un Interval
        for (auto &sp : to_insert) {
            auto sp_line = std::make_shared<HalfLine>(sp);
            Point isect = find_halflines_intersection(*p_line, *sp_line);
            if (std::isinf(isect.coord[0])) {
                continue;
            }
            double x_intersect = isect.coord[0];
            bool coversleft    = (sp_line->q < p_line->q);

            intervals.emplace_back(
                sp_line,
                std::make_pair(std::numeric_limits<double>::quiet_NaN(), x_intersect),
                coversleft
            );
        }

        if (!to_insert.empty()) {
            std::cout << "> " << to_insert.size()
                      << " halfline(s) have been inserted" << std::endl;
        }
    } // while (true)
}