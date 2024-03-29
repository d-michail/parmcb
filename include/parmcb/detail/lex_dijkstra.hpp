#ifndef PARMCB_DETAIL_LEX_DIJKSTRA_HPP_
#define PARMCB_DETAIL_LEX_DIJKSTRA_HPP_

//    Copyright (C) Dimitrios Michail 2019 - 2021.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#include <iostream>

#include <boost/scoped_array.hpp>
#include <boost/throw_exception.hpp>
#include <boost/functional/hash.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/property_map/function_property_map.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/detail/d_ary_heap.hpp>

#include <parmcb/detail/dijkstra.hpp>
#include <parmcb/detail/util.hpp>

namespace parmcb {

    namespace detail {

        template<class Graph, class DistanceMap>
        struct LexDistance {
            typedef typename boost::property_traits<DistanceMap>::value_type DistanceType;

            LexDistance() :
                    distance(DistanceType()), edge_count(0), vertex_indices() {
            }

            LexDistance(DistanceType distance, std::size_t edge_count, std::set<std::size_t> vertex_indices) :
                    distance(distance), edge_count(edge_count), vertex_indices(vertex_indices) {
            }

            DistanceType distance;
            std::size_t edge_count;
            std::set<std::size_t> vertex_indices;
        };

        template<class Graph, class DistanceMap>
        struct LexDistanceCompare {
            typedef typename boost::property_traits<DistanceMap>::value_type DistanceType;

            bool operator()(const LexDistance<Graph, DistanceMap> &a, const LexDistance<Graph, DistanceMap> &b) {
                if (a.distance < b.distance) {
                    return true;
                }
                if (a.distance > b.distance) {
                    return false;
                }
                if (a.edge_count < b.edge_count) {
                    return true;
                }
                if (a.edge_count > b.edge_count) {
                    return false;
                }

                std::set<std::size_t> non_common_a;
                std::set_difference(a.vertex_indices.begin(), a.vertex_indices.end(), 
                            b.vertex_indices.begin(), b.vertex_indices.end(),
                            std::inserter(non_common_a, non_common_a.end()));

                std::set<std::size_t> non_common_b;
                std::set_difference(b.vertex_indices.begin(), b.vertex_indices.end(), 
                            a.vertex_indices.begin(), a.vertex_indices.end(),
                            std::inserter(non_common_b, non_common_b.end()));

                if (non_common_a.empty() && !non_common_b.empty()) {
                    return true;
                }
                if (!non_common_a.empty() && non_common_b.empty()) {
                    return false;
                }

                // If non-common elements are different, compare them
                if (!non_common_a.empty() && !non_common_b.empty()) {
                    auto min_elem_a = *non_common_a.begin();
                    auto min_elem_b = *non_common_b.begin();
                    if (min_elem_a < min_elem_b) { 
                        return true;
                    }
                }                

                return false;
            }
        };

        template<class Graph, class IndexMap, class WeightMap, class DistanceMap>
        struct LexDistanceCombine {
            typedef typename boost::property_map<Graph, boost::vertex_index_t>::type VertexIndexMapType;
            typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
            typedef typename boost::property_traits<DistanceMap>::value_type DistanceType;
            typedef typename boost::property_traits<WeightMap>::value_type WeightType;

            LexDistanceCombine(const Graph &g, const IndexMap &index_map, const WeightMap &weight_map) :
                    inf((std::numeric_limits<DistanceType>::max)()), g(g), index_map(
                            index_map), weight_map(weight_map) {
            }

            LexDistance<Graph, DistanceMap> operator()(const LexDistance<Graph, DistanceMap> &a, const Edge &e) {
                const auto index_target = index_map[boost::target(e, g)];
                const auto index_source = index_map[boost::source(e, g)];

                const WeightType e_weight = boost::get(weight_map, e);
                const WeightType sum = combine(a.distance, e_weight);

                std::set<std::size_t> vertex_indices = a.vertex_indices;
                vertex_indices.insert(index_target);
                vertex_indices.insert(index_source);

                return LexDistance<Graph, DistanceMap>(sum, a.edge_count + 1, vertex_indices);
            }

            const DistanceType inf;
            const parmcb::detail::closed_plus<DistanceType> combine;
            const Graph &g;
            const VertexIndexMapType &index_map;
            const WeightMap &weight_map;
        };

    } // detail

    template<class Graph, class WeightMap, class DistanceMap, class PredecessorMap>
    void lex_dijkstra(const Graph &g, const WeightMap &weight_map,
            const typename boost::graph_traits<Graph>::vertex_descriptor &s, DistanceMap &dist_map,
            PredecessorMap &pred_map) {

        typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
        typedef typename boost::property_map<Graph, boost::vertex_index_t>::type VertexIndexMapType;
        typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
        typedef typename boost::property_traits<DistanceMap>::value_type DistanceType;

        typedef typename parmcb::detail::LexDistance<Graph, DistanceMap> LexDistanceType;
        typedef typename parmcb::detail::LexDistanceCompare<Graph, DistanceMap> LexDistanceCompareType;

        const VertexIndexMapType &index_map = boost::get(boost::vertex_index, g);
        std::vector<std::size_t> index_in_heap(boost::num_vertices(g));
        boost::function_property_map<parmcb::detail::VertexIndexFunctor<Graph, std::size_t>, Vertex, std::size_t&> index_in_heap_map(
                parmcb::detail::VertexIndexFunctor<Graph, std::size_t>(index_in_heap, index_map));

        std::vector<LexDistanceType> lex_dist(boost::num_vertices(g));
        boost::function_property_map<parmcb::detail::VertexIndexFunctor<Graph, LexDistanceType>, Vertex,
                LexDistanceType&> lex_dist_map(
                parmcb::detail::VertexIndexFunctor<Graph, LexDistanceType>(lex_dist, index_map));

        typedef boost::d_ary_heap_indirect<Vertex, 4,
                boost::function_property_map<parmcb::detail::VertexIndexFunctor<Graph, std::size_t>, Vertex,
                        std::size_t&>,
                boost::function_property_map<parmcb::detail::VertexIndexFunctor<Graph, LexDistanceType>, Vertex,
                        LexDistanceType&>, LexDistanceCompareType> VertexQueue;

        LexDistanceCompareType compare;
        parmcb::detail::LexDistanceCombine<Graph, VertexIndexMapType, WeightMap, DistanceMap> combine(g, index_map, weight_map);

        VertexQueue queue(lex_dist_map, index_in_heap_map, compare);

        boost::put(lex_dist_map, s, LexDistanceType(DistanceType(), 0, std::set<std::size_t>({ index_map[s] })));
        boost::put(pred_map, s, std::make_tuple(false, Edge()));
        queue.push(s);

        while (!queue.empty()) {
            Vertex u = queue.top();
            queue.pop();
            LexDistanceType d_u = boost::get(lex_dist_map, u);

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

                const LexDistanceType c = combine(d_u, e);
                bool visited_w = std::get<0>(boost::get(pred_map, w));
                if (!visited_w) {
                    // first time found
                    boost::put(lex_dist_map, w, c);
                    boost::put(dist_map, w, c.distance);
                    boost::put(pred_map, w, std::make_tuple(true, e));
                    queue.push(w);
                } else if (compare(c, boost::get(lex_dist_map, w))) {
                    // already reached
                    boost::put(lex_dist_map, w, c);
                    boost::put(dist_map, w, c.distance);
                    boost::put(pred_map, w, std::make_tuple(true, e));
                    queue.update(w);
                }
            }
        }
    }

} // mcb

#endif
