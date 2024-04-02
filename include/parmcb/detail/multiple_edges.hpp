#ifndef PARMCB_DETAIL_MULTIPLE_EDGES_HPP_
#define PARMCB_DETAIL_MULTIPLE_EDGES_HPP_

//    Copyright (C) Dimitrios Michail 2019 - 2024.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#include <parmcb/config.hpp>
#include <parmcb/detail/dijkstra.hpp>
#include <parmcb/detail/bfs.hpp>
#include <parmcb/detail/approx_spanner.hpp>

#include <functional>
#include <numeric>

#ifdef PARMCB_HAVE_TBB
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/concurrent_vector.h>
#include <tbb/task_group.h>
#endif

namespace parmcb {

namespace detail {

template<typename Graph, typename WeightMap, typename CycleOutputIterator>
using MCBAlgoType = std::function<typename boost::property_traits<WeightMap>::value_type(const Graph&, const WeightMap&, CycleOutputIterator)>;


template<class Graph, class WeightMap, class CycleOutputIterator, bool ParallelUsingTBB>
class MultipleEdgesReductionAlgorithm {
public:
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename boost::graph_traits<Graph>::vertex_iterator VertexIt;
    typedef typename boost::property_map<Graph, boost::vertex_index_t>::type VertexIndexMapType;
    typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
    typedef typename boost::graph_traits<Graph>::edge_iterator EdgeIt;
    typedef typename boost::property_traits<WeightMap>::value_type WeightType;
    typedef typename boost::property_map<Graph, boost::edge_weight_t>::type EdgeWeightMapType;

    MultipleEdgesReductionAlgorithm(const Graph &g, const WeightMap &weight_map, const VertexIndexMapType& index_map, MCBAlgoType<Graph, WeightMap, CycleOutputIterator> algo) :
            _g(g), _weight_map(weight_map), _index_map(index_map), _vertex_g_to_spanner_vec(
                    boost::num_vertices(g),
                    (std::numeric_limits<std::size_t>::max)()), _vertex_g_to_spanner(
                    parmcb::detail::VertexIndexFunctor<Graph, Vertex>(
                            _vertex_g_to_spanner_vec, _index_map)), _spanner(), _spanner_index_map(
                    boost::get(boost::vertex_index, _spanner)), _weight(
                    WeightType()), 
                    _algo(algo) {
        construct_spanner();
    }

    WeightType operator()(CycleOutputIterator out) {
#ifdef PARMCB_INVARIANTS_CHECK
        // preconditions check
        check_edge_length_preconditions();
#endif

        // compute spanner MCB with exact algorithm
        EdgeWeightMapType spanner_weight_map = get(boost::edge_weight,
                _spanner);
        _weight += _algo(_spanner, spanner_weight_map, out);

        // compute remaining cycles
        parmcb::detail::NonSpannerEdgesCycleBuilder<Graph, WeightMap,
                ParallelUsingTBB> non_spanner_edges_cycle_builder(_g,
                _weight_map, _spanner, _spanner_index_map, _edge_spanner_to_g,
                spanner_weight_map, _non_spanner_edges, _vertex_g_to_spanner);
        _weight += non_spanner_edges_cycle_builder(out);

        return _weight;
    }

private:
    // graph
    const Graph &_g;
    const WeightMap &_weight_map;
    const VertexIndexMapType &_index_map;
    std::vector<Vertex> _vertex_g_to_spanner_vec;
    const boost::function_property_map<
            parmcb::detail::VertexIndexFunctor<Graph, Vertex>, Vertex, Vertex&> _vertex_g_to_spanner;
    std::vector<Edge> _non_spanner_edges;

    // spanner
    Graph _spanner;
    VertexIndexMapType _spanner_index_map;
    std::map<Edge, Edge> _edge_spanner_to_g;

    // mcb
    WeightType _weight;

    // wrapped algo
    MCBAlgoType<Graph, WeightMap, CycleOutputIterator> _algo;

    void construct_spanner() {
        if (boost::num_vertices(_spanner) > 0) {
            throw std::runtime_error("Target graph is not empty");
        }

        /*
        * Find all multiple edges and leave only the one with the lowest weight
        * in the original graph.
        */
        boost::unordered_set<Edge> heavy_multiple_edges;

        auto vRange = boost::vertices(_g);
        for (auto vit = vRange.first; vit != vRange.second; ++vit) {
            auto v = *vit;
            boost::unordered_map<Vertex, Edge> min_edge_per_neighbor;

            auto eRange = boost::out_edges(v, _g);
            for (auto eit = eRange.first; eit != eRange.second; ++eit) {
                auto e = *eit;
                auto u = boost::opposite(e, v, _g);

                auto it = min_edge_per_neighbor.find(u);
                if (it == min_edge_per_neighbor.end()) {
                    // found first edge
                    min_edge_per_neighbor[u] = e;
                } else {
                    // found another
                    if (_weight_map[e] < _weight_map[it->second]) {
                        heavy_multiple_edges.insert(it->second);
                        min_edge_per_neighbor[u] = e;
                    } else { 
                        heavy_multiple_edges.insert(e);
                    }
                }
            }
        }

#ifdef PARMCB_LOGGING
        std::cout << "Found " << heavy_multiple_edges.size() << " heavy multiple edges..." << std::endl;
#endif

        /*
         * Copy graph without the heavy edges
         */
        VertexIt vi, vi_end;
        for (std::tie(vi, vi_end) = boost::vertices(_g); vi != vi_end; ++vi) {
            Vertex v = *vi;
            Vertex spanner_v = boost::add_vertex(_spanner);
            _vertex_g_to_spanner[v] = spanner_v;
        }

        // construct edge set of spanner
        EdgeIt ei, ei_end;
        for (boost::tie(ei, ei_end) = boost::edges(_g); ei != ei_end; ++ei) {
            Edge e = *ei;
            if (heavy_multiple_edges.find(e) != heavy_multiple_edges.end()) { 
                // skip heavy edges (record missing edge)
                _non_spanner_edges.push_back(e);                
                continue;
            }
            Vertex v = boost::source(e, _g);
            Vertex u = boost::target(e, _g);
            Vertex spanner_v = _vertex_g_to_spanner[v];
            Vertex spanner_u = _vertex_g_to_spanner[u];

            if (spanner_v == spanner_u) {
                throw std::runtime_error("Self loops?");
            }

            // add edge to spanner
            Edge spanner_e = std::get<0>(
                    boost::add_edge(spanner_v, spanner_u, _spanner));
            _edge_spanner_to_g[spanner_e] = e;
        }

#ifdef PARMCB_LOGGING
        std::cout << "Graph has " << boost::num_edges(_g) << " edges"
                << std::endl;
        std::cout << "Spanner has " << boost::num_edges(_spanner) << " edges"
                << std::endl;
#endif

    }

    void check_edge_length_preconditions() {
        EdgeIt ei, ei_end;
        for (std::tie(ei, ei_end) = boost::edges(_g); ei != ei_end; ++ei) {
            Edge e = *ei;
            if (_weight_map[e] < WeightType()) {
                throw std::runtime_error("Invalid edge weight.");
            }
        }
    }

};

} // detail

} // parmcb

#endif
