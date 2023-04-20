#ifndef PARMCB_DETAIL_APPROX_SPANNER_HPP_
#define PARMCB_DETAIL_APPROX_SPANNER_HPP_

//    Copyright (C) Dimitrios Michail 2019 - 2023.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <vector>

#include <boost/throw_exception.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/property_map/function_property_map.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <parmcb/forestindex.hpp>
#include <parmcb/detail/dijkstra.hpp>
#include <parmcb/detail/bfs.hpp>

namespace parmcb {

namespace detail {

template<class Graph, class WeightMap, typename ExactAlgorithm>
class BaseApproxSpannerAlgorithm {
public:
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename boost::graph_traits<Graph>::vertex_iterator VertexIt;
    typedef typename boost::property_map<Graph, boost::vertex_index_t>::type VertexIndexMapType;
    typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
    typedef typename boost::graph_traits<Graph>::edge_iterator EdgeIt;
    typedef typename std::vector<Edge>::iterator EdgeVectorIt;
    typedef typename boost::property_traits<WeightMap>::value_type WeightType;
    typedef typename boost::property_map<Graph, boost::edge_weight_t>::type EdgeWeightMapType;

    BaseApproxSpannerAlgorithm(const Graph &g, const WeightMap &weight_map,
            std::size_t k) :
            _g(g), _weight_map(weight_map), _k(k), _index_map(
                    boost::get(boost::vertex_index, g)), _vertex_g_to_spanner_vec(
                    boost::num_vertices(g),
                    (std::numeric_limits<std::size_t>::max)()), _vertex_g_to_spanner(
                    parmcb::detail::VertexIndexFunctor<Graph, Vertex>(
                            _vertex_g_to_spanner_vec, _index_map)), _spanner(), _spanner_index_map(
                    boost::get(boost::vertex_index, _spanner)), _weight(
                    WeightType()) {
        construct_spanner();
    }

    template<class CycleOutputIterator>
    WeightType run(CycleOutputIterator out) {
        // preconditions check
        if (_k < 1) {
            throw std::runtime_error(
                    "Invalid value of k < 1 for spanner construction.");
        }
        check_edge_length_preconditions();

        // compute spanner MCB
        EdgeWeightMapType spanner_weight_map = get(boost::edge_weight,
                _spanner);
        ExactAlgorithm exact_mcb_algo;
        _weight += exact_mcb_algo(_spanner, spanner_weight_map, out);

        // compute remaining cycles
        _weight += construct_cycles_for_non_spanner_edges(out);

        return _weight;
    }

private:
    // graph
    const Graph &_g;
    const WeightMap &_weight_map;
    const std::size_t _k;
    const VertexIndexMapType &_index_map;
    std::vector<Vertex> _vertex_g_to_spanner_vec;
    const boost::function_property_map<
            parmcb::detail::VertexIndexFunctor<Graph, Vertex>, Vertex, Vertex&> _vertex_g_to_spanner;
    std::list<Edge> _non_spanner_edges;

    // spanner
    Graph _spanner;
    VertexIndexMapType _spanner_index_map;
    std::map<Edge, Edge> _edge_spanner_to_g;

    // mcb
    WeightType _weight;

    void construct_spanner() {
        if (boost::num_vertices(_spanner) > 0) {
            throw std::runtime_error("Target graph is not empty");
        }

        // sort graph edges
        std::vector<Edge> sorted_edges;
        auto eItPair = boost::edges(_g);
        std::copy(eItPair.first, eItPair.second,
                std::back_inserter(sorted_edges));
        std::sort(sorted_edges.begin(), sorted_edges.end(),
                [&](const Edge &e1, const Edge &e2) {
                    return _weight_map[e1] < _weight_map[e2];
                });

        // create vertex set of spanner
        VertexIt vi, vi_end;
        for (std::tie(vi, vi_end) = boost::vertices(_g); vi != vi_end; ++vi) {
            Vertex v = *vi;
            Vertex spanner_v = boost::add_vertex(_spanner);
            _vertex_g_to_spanner[v] = spanner_v;
        }

        // construct edge set of spanner
        EdgeVectorIt ei, ei_end;
        for (ei = sorted_edges.begin(), ei_end = sorted_edges.end();
                ei != sorted_edges.end(); ++ei) {
            Edge e = *ei;
            Vertex v = boost::source(e, _g);
            Vertex u = boost::target(e, _g);
            Vertex spanner_v = _vertex_g_to_spanner[v];
            Vertex spanner_u = _vertex_g_to_spanner[u];

            if (spanner_v == spanner_u) {
                throw std::runtime_error("Self loops?");
            }

            if (!parmcb::is_bfs_reachable(_spanner, spanner_v, spanner_u,
                    2 * _k - 1)) {
                // add edge to spanner
                Edge spanner_e = std::get<0>(
                        boost::add_edge(spanner_v, spanner_u, _spanner));
                _edge_spanner_to_g[spanner_e] = e;
            } else {
                // record missing edge from spanner
                _non_spanner_edges.push_back(e);
            }
        }

        std::cout << "Graph has " << boost::num_edges(_g) << " edges"
                << std::endl;
        std::cout << "Spanner has " << boost::num_edges(_spanner) << " edges"
                << std::endl;

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

    template<class CycleOutputIterator>
    WeightType construct_cycles_for_non_spanner_edges(CycleOutputIterator out) {
        WeightType total_weight = WeightType();

        for (auto it = _non_spanner_edges.begin();
                it != _non_spanner_edges.end(); it++) {
            auto e = *it;
            Vertex v = boost::source(e, _g);
            Vertex u = boost::target(e, _g);
            Vertex spanner_v = _vertex_g_to_spanner[v];
            Vertex spanner_u = _vertex_g_to_spanner[u];

            // compute shortest path on spanner
            std::vector<WeightType> dist(boost::num_vertices(_spanner),
                    (std::numeric_limits<WeightType>::max)());
            boost::function_property_map<
                    parmcb::detail::VertexIndexFunctor<Graph, WeightType>,
                    Vertex, WeightType&> dist_map(
                    parmcb::detail::VertexIndexFunctor<Graph, WeightType>(dist,
                            _spanner_index_map));
            std::vector<std::tuple<bool, Edge>> pred(
                    boost::num_vertices(_spanner),
                    std::make_tuple(false, Edge()));
            boost::function_property_map<
                    parmcb::detail::VertexIndexFunctor<Graph,
                            std::tuple<bool, Edge>>, Vertex,
                    std::tuple<bool, Edge>&> pred_map(
                    parmcb::detail::VertexIndexFunctor<Graph,
                            std::tuple<bool, Edge> >(pred, _spanner_index_map));
            EdgeWeightMapType spanner_weight_map = get(boost::edge_weight,
                    _spanner);

            // run dijkstra
            parmcb::dijkstra(_spanner, spanner_weight_map, spanner_v, dist_map,
                    pred_map);

            // form cycle
            std::list<Edge> cycle_edgelist;
            WeightType weight = WeightType();
            Vertex spanner_w = spanner_u;
            while (true) {
                auto pred_t = boost::get(pred_map, spanner_w);
                if (!std::get<0>(pred_t)) {
                    // no predecessor
                    break;
                }
                Edge spanner_ae = std::get<1>(pred_t);
                Edge ae = _edge_spanner_to_g[spanner_ae];
                cycle_edgelist.push_back(ae);
                weight += boost::get(_weight_map, ae);

                // go to predecessor
                auto spanner_other = boost::target(spanner_ae, _spanner);
                if (spanner_other == spanner_w) {
                    spanner_other = boost::source(spanner_ae, _spanner);
                }
                if (spanner_other == spanner_w) {
                    throw new std::runtime_error("Self loops?");
                }
                spanner_w = spanner_other;
            }
            cycle_edgelist.push_back(e);
            weight += boost::get(_weight_map, e);

            // output
            *out++ = cycle_edgelist;
            total_weight += weight;
        }

        return total_weight;
    }

};

} // detail

} // parmcb

#endif
