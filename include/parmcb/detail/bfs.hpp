#ifndef PARMCB_DETAIL_UBFS_HPP_
#define PARMCB_DETAIL_UBFS_HPP_

//    Copyright (C) Dimitrios Michail 2019 - 2023.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <queue>

#include <boost/scoped_array.hpp>
#include <boost/throw_exception.hpp>
#include <boost/functional/hash.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/property_map/function_property_map.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <parmcb/detail/util.hpp>

namespace parmcb {

template<class Graph>
bool is_bfs_reachable(const Graph &g,
        const typename boost::graph_traits<Graph>::vertex_descriptor &s,
        const typename boost::graph_traits<Graph>::vertex_descriptor &t,
        std::size_t max_hops) {
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
    typedef std::queue<Vertex> VertexQueue;

    std::vector<std::size_t> dist(boost::num_vertices(g),
            (std::numeric_limits<std::size_t>::max)());
    auto index_map = boost::get(boost::vertex_index, g);
    boost::function_property_map<
            parmcb::detail::VertexIndexFunctor<Graph, std::size_t>, Vertex,
            std::size_t&> dist_map(
            parmcb::detail::VertexIndexFunctor<Graph, std::size_t>(dist,
                    index_map));

    std::vector<std::tuple<bool, Edge>> pred(boost::num_vertices(g),
            std::make_tuple(false, Edge()));
    boost::function_property_map<
            parmcb::detail::VertexIndexFunctor<Graph, std::tuple<bool, Edge>>,
            Vertex, std::tuple<bool, Edge>&> pred_map(
            parmcb::detail::VertexIndexFunctor<Graph, std::tuple<bool, Edge> >(
                    pred, index_map));

    parmcb::detail::closed_plus<std::size_t> combine =
            parmcb::detail::closed_plus<std::size_t>();

    VertexQueue queue;
    boost::put(dist_map, s, size_t());
    boost::put(pred_map, s, std::make_tuple(false, Edge()));
    queue.push(s);

    while (!queue.empty()) {
        Vertex u = queue.front();
        queue.pop();
        size_t d_u = boost::get(dist_map, u);

        if (d_u > max_hops) {
            return false;
        }

        if (u == t) {
            return true;
        }

        auto eiRange = boost::out_edges(u, g);
        for (auto ei = eiRange.first; ei != eiRange.second; ++ei) {
            auto e = *ei;

            auto w = boost::target(e, g);
            if (w == u) {
                w = boost::source(e, g);
            }
            if (w == u) {
                // self-loop
                continue;
            }
            if (w == s) {
                continue;
            }

            const std::size_t c = combine(d_u, 1);
            bool visited_w = std::get<0>(boost::get(pred_map, w));
            if (!visited_w) {
                // first time found
                boost::put(dist_map, w, c);
                boost::put(pred_map, w, std::make_tuple(true, e));
                queue.push(w);
            }
        }
    }

    return false;
}

template<class Graph, class DistanceMap, class PredecessorMap>
void bfs(const Graph &g,
        const typename boost::graph_traits<Graph>::vertex_descriptor &s,
        DistanceMap &dist_map, PredecessorMap &pred_map) {

    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
    typedef typename boost::property_traits<DistanceMap>::value_type DistanceType;
    typedef std::queue<Vertex> VertexQueue;

    std::less<DistanceType> compare;
    parmcb::detail::closed_plus<DistanceType> combine =
            parmcb::detail::closed_plus<DistanceType>();

    VertexQueue queue;

    boost::put(dist_map, s, DistanceType());
    boost::put(pred_map, s, std::make_tuple(false, Edge()));
    queue.push(s);

    while (!queue.empty()) {
        Vertex u = queue.front();
        queue.pop();
        DistanceType d_u = boost::get(dist_map, u);

        auto eiRange = boost::out_edges(u, g);
        for (auto ei = eiRange.first; ei != eiRange.second; ++ei) {
            auto e = *ei;

            auto w = boost::target(e, g);
            if (w == u) {
                w = boost::source(e, g);
            }
            if (w == u) {
                // self-loop
                continue;
            }
            if (w == s) {
                continue;
            }

            const DistanceType c = combine(d_u, DistanceType(1));
            bool visited_w = std::get<0>(boost::get(pred_map, w));
            if (!visited_w) {
                // first time found
                boost::put(dist_map, w, c);
                boost::put(pred_map, w, std::make_tuple(true, e));
                queue.push(w);
            }
        }
    }
}

} // parmcb

#endif
