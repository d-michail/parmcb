#ifndef MCB_SPTREES_HPP_
#define MCB_SPTREES_HPP_

#include <iostream>

#include <boost/throw_exception.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/property_map/function_property_map.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/serialization/vector.hpp>

#include <mcb/forestindex.hpp>
#include <mcb/dijkstra.hpp>
#include <mcb/spvecgf2.hpp>
#include <mcb/util.hpp>

#include <memory>
#include <stack>
#include <functional>

#include <tbb/tbb.h>

namespace mcb {

    template<class Graph, class WeightMap> class SPNode;
    template<class Graph, class WeightMap> class SPTree;
    template<class Graph, class WeightMap> struct SPSubtree;
    template<class Graph, class WeightMap, bool ParallelUsingTBB> class SPTrees;
    template<class Graph, class WeightMap> class CandidateCycle;
    template<class Graph> struct SerializableCandidateCycle;
    template<class Graph, class WeightMap> struct SerializableMinOddCycle;
    template<class Graph, class WeightMap> struct SerializableMinOddCycleMinOp;

    template<class Graph, class WeightMap>
    class SPNode {
    public:
        typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
        typedef typename boost::property_traits<WeightMap>::value_type WeightType;

        SPNode() :
                _parity(false), _weight(WeightType()), _pred(), _has_pred(false) {
        }

        SPNode(WeightType weight) :
                _parity(false), _weight(WeightType()), _pred(), _has_pred(false) {
        }

        SPNode(WeightType weight, const Edge &pred) :
                _parity(false), _weight(weight), _pred(pred), _has_pred(true) {
        }

        void add_child(std::shared_ptr<SPNode<Graph, WeightMap>> c) {
            _children.push_back(c);
        }

        std::vector<std::shared_ptr<SPNode<Graph, WeightMap>>>& children() {
            return _children;
        }

        bool& parity() {
            return _parity;
        }

        WeightType& weight() {
            return _weight;
        }

        const Edge& pred() {
            return _pred;
        }

        bool has_pred() {
            return _has_pred;
        }

    private:
        bool _parity;
        WeightType _weight;
        Edge _pred;
        bool _has_pred;
        std::vector<std::shared_ptr<SPNode<Graph, WeightMap>>> _children;
    };

    template<class Graph, class WeightMap>
    struct SPSubtree {
        bool parity;
        std::shared_ptr<SPNode<Graph, WeightMap>> root;

        SPSubtree(bool parity, std::shared_ptr<SPNode<Graph, WeightMap>> root) :
                parity(parity), root(root) {
        }
    };

    template<class Graph, class WeightMap>
    class SPTree {
    public:
        typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
        typedef typename boost::graph_traits<Graph>::vertex_iterator VertexIt;
        typedef typename boost::property_map<Graph, boost::vertex_index_t>::type VertexIndexMapType;
        typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
        typedef typename boost::property_traits<WeightMap>::value_type WeightType;

        SPTree(std::size_t id, const Graph &g, const WeightMap &weight_map, const Vertex &source) :
                _id(id), _g(g), _weight_map(weight_map), _index_map(boost::get(boost::vertex_index, g)), _source(
                        source), _tree_node_map(boost::num_vertices(g)) {
            initialize();
        }

        void update_parities(const std::set<Edge> &edges) {
            std::stack<SPSubtree<Graph, WeightMap>> stack;
            stack.emplace(false, _root);

            while (!stack.empty()) {
                SPSubtree<Graph, WeightMap> r = stack.top();
                stack.pop();

                r.root->parity() = r.parity;
                for (auto c : r.root->children()) {
                    bool is_signed = edges.find(c->pred()) != edges.end();
                    stack.emplace(SPSubtree<Graph, WeightMap> { static_cast<bool>(r.parity ^ is_signed), c });
                }
            }
        }

        std::shared_ptr<SPNode<Graph, WeightMap>> node(const Vertex &v) const {
            return _tree_node_map[_index_map[v]];
        }

        const Vertex& source() {
            return _source;
        }

        const Graph& graph() {
            return _g;
        }

        const std::size_t id() const {
            return _id;
        }

        template<class EdgeIterator>
        std::vector<CandidateCycle<Graph, WeightMap>> create_candidate_cycles(EdgeIterator begin, EdgeIterator end) {
            // collect tree edges
            std::set<Edge> tree_edges;
            VertexIt vi, viend;
            for (boost::tie(vi, viend) = boost::vertices(_g); vi != viend; ++vi) {
                auto v = *vi;
                auto vindex = _index_map[v];
                std::shared_ptr<SPNode<Graph, WeightMap>> n = _tree_node_map[vindex];
                if (n != nullptr && n->has_pred()) {
                    tree_edges.insert(n->pred());
                }
            }

            // loop over (non-tree) provided edges and create candidate cycles
            std::vector<CandidateCycle<Graph, WeightMap>> cycles;
            for (EdgeIterator it = begin; it != end; it++) {
                Edge e = *it;
                if (tree_edges.find(e) != tree_edges.end()) {
                    continue;
                }

                // non-tree edge
                std::shared_ptr<SPNode<Graph, WeightMap>> v = node(boost::source(e, _g));
                if (v == nullptr) {
                    continue;
                }
                std::shared_ptr<SPNode<Graph, WeightMap>> u = node(boost::target(e, _g));
                if (u == nullptr) {
                    continue;
                }

                WeightType cycle_weight = boost::get(_weight_map, e) + v->weight() + u->weight();
                cycles.emplace_back(_id, e, cycle_weight);
            }
            return cycles;
        }

        std::vector<CandidateCycle<Graph, WeightMap>> create_candidate_cycles() {
            auto itPair = boost::edges(_g);
            return create_candidate_cycles(itPair.first, itPair.second);
        }

        std::vector<SerializableCandidateCycle<Graph>> create_serializable_candidate_cycles(
                const ForestIndex<Graph> &forest_index) {
            // collect tree edges
            std::set<Edge> tree_edges;
            VertexIt vi, viend;
            for (boost::tie(vi, viend) = boost::vertices(_g); vi != viend; ++vi) {
                auto v = *vi;
                auto vindex = _index_map[v];
                std::shared_ptr<SPNode<Graph, WeightMap>> n = _tree_node_map[vindex];
                if (n != nullptr && n->has_pred()) {
                    tree_edges.insert(n->pred());
                }
            }

            // loop over all non-tree edges and create candidate cycles
            std::vector<SerializableCandidateCycle<Graph>> cycles;
            for (const auto &e : boost::make_iterator_range(boost::edges(_g))) {
                if (tree_edges.find(e) != tree_edges.end()) {
                    continue;
                }

                // non-tree edge
                std::shared_ptr<SPNode<Graph, WeightMap>> v = node(boost::source(e, _g));
                if (v == nullptr) {
                    continue;
                }
                std::shared_ptr<SPNode<Graph, WeightMap>> u = node(boost::target(e, _g));
                if (u == nullptr) {
                    continue;
                }
                cycles.emplace_back(_source, forest_index(e));
            }

            return cycles;
        }

    private:
        const std::size_t _id;
        const Graph &_g;
        const WeightMap &_weight_map;
        const VertexIndexMapType _index_map;
        const Vertex _source;

        std::vector<std::shared_ptr<SPNode<Graph, WeightMap>>> _tree_node_map;
        std::shared_ptr<SPNode<Graph, WeightMap>> _root;

        void initialize() {
            // run shortest path
            std::vector<WeightType> dist(boost::num_vertices(_g), (std::numeric_limits<WeightType>::max)());
            boost::function_property_map<mcb::detail::VertexIndexFunctor<Graph, WeightType>, Vertex, WeightType&> dist_map(
                    mcb::detail::VertexIndexFunctor<Graph, WeightType>(dist, _index_map));
            std::vector<std::tuple<bool, Edge>> pred(boost::num_vertices(_g), std::make_tuple(false, Edge()));
            boost::function_property_map<mcb::detail::VertexIndexFunctor<Graph, std::tuple<bool, Edge>>, Vertex,
                    std::tuple<bool, Edge>&> pred_map(
                    mcb::detail::VertexIndexFunctor<Graph, std::tuple<bool, Edge> >(pred, _index_map));
            dijkstra(_g, _weight_map, _source, dist_map, pred_map);

            // create tree nodes and mapping
            VertexIt vi, viend;
            for (boost::tie(vi, viend) = boost::vertices(_g); vi != viend; ++vi) {
                auto v = *vi;
                auto vindex = _index_map[v];
                auto p = boost::get(pred_map, v);
                if (v == _source) {
                    _tree_node_map[vindex] = std::shared_ptr<SPNode<Graph, WeightMap>>(
                            new SPNode<Graph, WeightMap>(dist[vindex]));
                    _root = _tree_node_map[vindex];
                } else if (std::get<0>(p)) {
                    Edge e = std::get<1>(p);
                    _tree_node_map[vindex] = std::shared_ptr<SPNode<Graph, WeightMap>>(
                            new SPNode<Graph, WeightMap>(dist[vindex], e));
                }
            }

            // link tree nodes
            for (boost::tie(vi, viend) = boost::vertices(_g); vi != viend; ++vi) {
                auto v = *vi;
                auto p = boost::get(pred_map, v);
                if (std::get<0>(p)) {
                    auto e = std::get<1>(p);
                    auto u = boost::opposite(e, v, _g);
                    auto vindex = _index_map[v];
                    auto uindex = _index_map[u];
                    _tree_node_map[uindex]->add_child(_tree_node_map[vindex]);
                }
            }
        }
    };

    template<class Graph, class WeightMap>
    class CandidateCycle {
    public:
        typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
        typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
        typedef typename boost::property_traits<WeightMap>::value_type WeightType;

        CandidateCycle(std::size_t tree, const Edge &e, WeightType weight) :
                _tree(tree), _e(e), _weight(weight) {
        }

        CandidateCycle(const CandidateCycle &c) :
                _tree(c._tree), _e(c._e), _weight(c._weight) {
        }

        CandidateCycle& operator=(const CandidateCycle &other) {
            if (this != &other) {
                _tree = other._tree;
                _e = other._e;
                _weight = other._weight;
            }
            return *this;
        }

        std::size_t tree() const {
            return _tree;
        }

        const Edge& edge() const {
            return _e;
        }

        const WeightType& weight() const {
            return _weight;
        }

    private:
        std::size_t _tree;
        Edge _e;
        WeightType _weight;
    };

    template<class Graph>
    struct SerializableCandidateCycle {
        typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
        typedef typename ForestIndex<Graph>::size_type Edge;

        SerializableCandidateCycle() {
        }

        SerializableCandidateCycle(Vertex v, Edge e) :
                v(v), e(e) {
        }

        template<typename Archive>
        void serialize(Archive &ar, const unsigned) {
            ar & v;
            ar & e;
        }

        Vertex v;
        Edge e;
    };

    template<class Graph, class WeightMap>
    struct SerializableMinOddCycle {
        typedef typename ForestIndex<Graph>::size_type Edge;
        typedef typename boost::property_traits<WeightMap>::value_type WeightType;

        SerializableMinOddCycle() :
                exists(false) {
        }

        SerializableMinOddCycle(std::vector<Edge> edges, WeightType weight, bool exists) :
                edges(edges), weight(weight), exists(exists) {
        }

        SerializableMinOddCycle(const SerializableMinOddCycle<Graph, WeightMap> &c) :
                edges(c.edges), weight(c.weight), exists(c.exists) {
        }

        SerializableMinOddCycle<Graph, WeightMap>& operator=(const SerializableMinOddCycle<Graph, WeightMap> &other) {
            if (this != &other) {
                edges = other.edges;
                weight = other.weight;
                exists = other.exists;
            }
            return *this;
        }

        template<typename Archive>
        void serialize(Archive &ar, const unsigned) {
            ar & edges;
            ar & weight;
            ar & exists;
        }

        std::vector<Edge> edges;
        WeightType weight;
        bool exists;
    };

    template<class Graph, class WeightMap>
    struct SerializableMinOddCycleMinOp {

        const SerializableMinOddCycle<Graph, WeightMap>& operator()(
                const SerializableMinOddCycle<Graph, WeightMap> &lhs,
                const SerializableMinOddCycle<Graph, WeightMap> &rhs) const {
            if (!lhs.exists || !rhs.exists) {
                if (lhs.exists) {
                    return lhs;
                } else {
                    return rhs;
                }
            }
            // both valid, compare
            if (lhs.weight < rhs.weight) {
                return lhs;
            }
            return rhs;
        }

    };

    template<class Graph, class WeightMap, bool ParallelUsingTBB>
    class SPTrees {
    public:
        typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
        typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
        typedef typename boost::property_traits<WeightMap>::value_type WeightType;

        SPTrees(const Graph &g, const WeightMap &weight_map, const std::vector<Vertex> &fvs, bool sorted_cycles) :
                g(g), weight_map(weight_map), sorted_cycles(sorted_cycles) {
            build_trees(fvs);
        }

        SPTrees(const Graph &g, const WeightMap &weight_map, const std::map<Vertex, std::vector<Edge>> &fvs,
                bool sorted_cycles) :
                g(g), weight_map(weight_map), sorted_cycles(sorted_cycles) {
            build_trees(fvs);
        }

        std::tuple<std::set<Edge>, WeightType, bool> compute_shortest_odd_cycle(const std::set<Edge> &edges) {
            return _compute_shortest_odd_cycle(edges);
        }

        std::size_t get_number_of_trees() {
            return trees.size();
        }

        mcb::SPTree<Graph, WeightMap>& get_tree(std::size_t n) {
            return trees.at(n);
        }

    private:
        const Graph &g;
        const WeightMap &weight_map;

        std::vector<mcb::SPTree<Graph, WeightMap>> trees;
        std::vector<CandidateCycle<Graph, WeightMap>> cycles;
        bool sorted_cycles;

        void build_trees(const std::map<Vertex, std::vector<Edge>> &fvs) {
            for (auto const &p : fvs) {
                SPTree<Graph, WeightMap> tree(trees.size(), g, weight_map, p.first);
                trees.push_back(tree);
                std::vector<CandidateCycle<Graph, WeightMap>> tree_cycles = tree.create_candidate_cycles(
                        p.second.begin(), p.second.end());
                cycles.insert(cycles.end(), tree_cycles.begin(), tree_cycles.end());
            }
            if (sorted_cycles) {
                // sort
                std::cout << "Sorting cycles" << std::endl;
                std::sort(cycles.begin(), cycles.end(), [](const auto &a, const auto &b) {
                    return a.weight() < b.weight();
                });
            }
        }

        void build_trees(const std::vector<Vertex> &fvs) {
            for (auto v : fvs) {
                trees.emplace_back(trees.size(), g, weight_map, v);
            }
            for (auto &tree : trees) {
                std::vector<CandidateCycle<Graph, WeightMap>> tree_cycles = tree.create_candidate_cycles();
                cycles.insert(cycles.end(), tree_cycles.begin(), tree_cycles.end());
            }
            if (sorted_cycles) {
                // sort
                std::cout << "Sorting cycles" << std::endl;
                std::sort(cycles.begin(), cycles.end(), [](const auto &a, const auto &b) {
                    return a.weight() < b.weight();
                });
            }
        }

        template<bool is_tbb_enabled = ParallelUsingTBB>
        std::tuple<std::set<Edge>, WeightType, bool> _compute_shortest_odd_cycle(const std::set<Edge> &edges,
                typename std::enable_if<!is_tbb_enabled>::type* = 0) {

            for (auto &tree : trees) {
                tree.update_parities(edges);
            }

            std::tuple<std::set<Edge>, WeightType, bool> min;

            for (CandidateCycle<Graph, WeightMap> c : cycles) {
                std::tuple<std::set<Edge>, WeightType, bool> cc = check_and_construct_candidate_cycle(c, edges,
                        std::get<2>(min), std::get<1>(min));

                if (std::get<2>(cc)) {
                    if (sorted_cycles) {
                        return cc;
                    }

                    if (!std::get<2>(min)) {
                        min = cc;
                    } else {
                        if (std::get<1>(cc) < std::get<1>(min)) {
                            min = cc;
                        }
                    }
                }
            }
            return min;
        }

        template<bool is_tbb_enabled = ParallelUsingTBB>
        std::tuple<std::set<Edge>, WeightType, bool> _compute_shortest_odd_cycle(const std::set<Edge> &edges,
                typename std::enable_if<is_tbb_enabled>::type* = 0) {

            tbb::parallel_for(tbb::blocked_range<std::size_t>(0, trees.size()),
                    [&](const tbb::blocked_range<std::size_t> &r) {
                        for (std::size_t i = r.begin(); i != r.end(); ++i) {
                            trees[i].update_parities(edges);
                        }
                    });

            std::less<WeightType> compare = std::less<WeightType>();
            typedef std::tuple<std::set<Edge>, WeightType, bool> cycle_t;
            auto cycle_min = [compare](const cycle_t &c1, const cycle_t &c2) {
                if (!std::get<2>(c1) || !std::get<2>(c2)) {
                    if (std::get<2>(c1)) {
                        return c1;
                    } else {
                        return c2;
                    }
                }
                // both valid, compare
                if (!compare(std::get<1>(c2), std::get<1>(c1))) {
                    return c1;
                }
                return c2;
            };

            return tbb::parallel_reduce(tbb::blocked_range<std::size_t>(0, cycles.size()),
                    std::make_tuple(std::set<Edge>(), (std::numeric_limits<WeightType>::max)(), false),
                    [&](tbb::blocked_range<std::size_t> r, auto running_min) {
                        for (std::size_t i = r.begin(); i < r.end(); i++) {
                            auto c = cycles[i];
                            auto cc = check_and_construct_candidate_cycle(c, edges, std::get<2>(running_min),
                                    std::get<1>(running_min));
                            if (std::get<2>(cc)) {
                                if (!std::get<2>(running_min) || compare(std::get<1>(cc), std::get<1>(running_min))) {
                                    running_min = cc;
                                }
                            }
                        }
                        return running_min;
                    },
                    cycle_min);
        }

        std::tuple<std::set<Edge>, WeightType, bool> check_and_construct_candidate_cycle(
                CandidateCycle<Graph, WeightMap> &c, const std::set<Edge> &edges, bool use_weight_limit,
                WeightType weight_limit) {
            //std::cout << "Checking candidate cycle from source: " << trees[c.tree()].source() << " with edge " << c.edge() << std::endl;

            std::shared_ptr<SPNode<Graph, WeightMap>> v = trees[c.tree()].node(boost::source(c.edge(), g));
            std::shared_ptr<SPNode<Graph, WeightMap>> u = trees[c.tree()].node(boost::target(c.edge(), g));

            Edge e = c.edge();
            if (v->parity() ^ u->parity() ^ (edges.find(e) != edges.end())) {
                // odd cycle, validate
                bool valid = true;
                WeightType cycle_weight = boost::get(weight_map, e);
                std::set<Edge> result;
                result.insert(e);

                if (use_weight_limit && cycle_weight > weight_limit) {
                    return std::make_tuple(std::set<Edge> { }, 0.0, false);
                }

                // first part
                Vertex w = boost::source(c.edge(), g);
                std::shared_ptr<SPNode<Graph, WeightMap>> ws = trees[c.tree()].node(w);
                while (ws->has_pred()) {
                    Edge a = ws->pred();
                    if (result.insert(a).second == false) {
                        valid = false;
                        break;
                    }
                    cycle_weight += boost::get(weight_map, a);
                    if (use_weight_limit && cycle_weight > weight_limit) {
                        valid = false;
                        break;
                    }
                    w = boost::opposite(a, w, g);
                    ws = trees[c.tree()].node(w);
                }

                if (!valid) {
                    return std::make_tuple(std::set<Edge> { }, 0.0, false);
                }

                // second part
                w = boost::target(c.edge(), g);
                ws = trees[c.tree()].node(w);
                while (ws->has_pred()) {
                    Edge a = ws->pred();
                    if (result.insert(a).second == false) {
                        valid = false;
                        break;
                    }
                    cycle_weight += boost::get(weight_map, a);
                    if (use_weight_limit && cycle_weight > weight_limit) {
                        valid = false;
                        break;
                    }
                    w = boost::opposite(a, w, g);
                    ws = trees[c.tree()].node(w);
                }

                if (!valid) {
                    return std::make_tuple(std::set<Edge> { }, 0.0, false);
                }

                return std::make_tuple(result, cycle_weight, true);
            }
            return std::make_tuple(std::set<Edge> { }, 0.0, false);
        }

    };

} // mcb

#endif
