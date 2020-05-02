#ifndef PARMCB_SIGNED_DIJKSTRA_HPP_
#define PARMCB_SIGNED_DIJKSTRA_HPP_

//    Copyright (C) Dimitrios Michail 2019 - 2020.
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

#include <parmcb/detail/util.hpp>

namespace std {

    std::ostream& operator <<(std::ostream &out, const std::pair<unsigned long int, bool> &v) {
        out << v.first << (v.second ? "+" : "-");
        return out;
    }

}

namespace parmcb {

    namespace detail {

        template<class Graph, class WeightMap>
        struct SignedDistanceFunctor {
            typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
            typedef typename boost::property_traits<WeightMap>::value_type DistanceType;
            typedef typename boost::property_map<Graph, boost::vertex_index_t>::type VertexIndexMapType;
            typedef std::pair<Vertex, bool> SignedVertex;

            std::size_t n;
            std::vector<DistanceType> &dist;
            const VertexIndexMapType &index_map;

            SignedDistanceFunctor(std::size_t n, std::vector<DistanceType> &dist, const VertexIndexMapType &index_map) :
                    n(n), dist(dist), index_map(index_map) {
            }

            DistanceType& operator()(const SignedVertex &v) const {
                return dist.at(index_map[v.first] + (v.second ? 0 : n));
            }
        };

        template<class Graph>
        struct SignedPredecessorFunctor {
            typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
            typedef typename boost::property_map<Graph, boost::vertex_index_t>::type VertexIndexMapType;
            typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
            typedef std::pair<Vertex, bool> SignedVertex;
            typedef std::tuple<SignedVertex, bool, Edge> Predecessor;

            std::size_t n;
            std::vector<Predecessor> &pred;
            const VertexIndexMapType &index_map;

            SignedPredecessorFunctor(std::size_t n, std::vector<Predecessor> &pred, const VertexIndexMapType &index_map) :
                    n(n), pred(pred), index_map(index_map) {
            }

            Predecessor& operator()(const SignedVertex &v) const {
                return pred.at(index_map[v.first] + (v.second ? 0 : n));
            }
        };

        template<class Graph>
        struct SignedIndexInHeapFunctor {
            typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
            typedef typename boost::property_map<Graph, boost::vertex_index_t>::type VertexIndexMapType;
            typedef std::pair<Vertex, bool> SignedVertex;

            std::size_t n;
            std::vector<std::size_t> &index_in_heap;
            const VertexIndexMapType &index_map;

            SignedIndexInHeapFunctor(std::size_t n, std::vector<std::size_t> &index_in_heap,
                    const VertexIndexMapType &index_map) :
                    n(n), index_in_heap(index_in_heap), index_map(index_map) {
            }

            std::size_t& operator()(const SignedVertex &v) const {
                return index_in_heap.at(index_map[v.first] + (v.second ? 0 : n));
            }
        };

        template<class Graph, class WeightMap>
        struct search_frontier {
            typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
            typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
            typedef typename boost::property_traits<WeightMap>::value_type WeightType;
            typedef typename boost::property_traits<WeightMap>::value_type DistanceType;
            typedef typename boost::property_map<Graph, boost::vertex_index_t>::type VertexIndexMapType;
            typedef std::pair<Vertex, bool> SignedVertex;
            typedef std::tuple<SignedVertex, bool, Edge> Predecessor;
            typedef boost::d_ary_heap_indirect<SignedVertex, 4,
                    boost::function_property_map<parmcb::detail::SignedIndexInHeapFunctor<Graph>, SignedVertex,
                            std::size_t&>,
                    boost::function_property_map<parmcb::detail::SignedDistanceFunctor<Graph, WeightMap>, SignedVertex,
                            DistanceType&>, std::less<DistanceType>> VertexQueue;

            const Graph &g;
            const VertexIndexMapType index_map;
            std::vector<std::size_t> index_in_heap;
            std::vector<DistanceType> dist;
            std::vector<Predecessor> pred;
            boost::function_property_map<parmcb::detail::SignedIndexInHeapFunctor<Graph>, SignedVertex, std::size_t&> index_in_heap_map;
            boost::function_property_map<parmcb::detail::SignedDistanceFunctor<Graph, WeightMap>, SignedVertex,
                    DistanceType&> dist_map;
            boost::function_property_map<parmcb::detail::SignedPredecessorFunctor<Graph>, SignedVertex, Predecessor&> pred_map;
            std::less<DistanceType> compare;
            VertexQueue queue;
            SignedVertex source;

            search_frontier(const Graph &g) :
                    g(g), index_map(boost::get(boost::vertex_index, g)), index_in_heap(2 * boost::num_vertices(g)), dist(
                            2 * boost::num_vertices(g), (std::numeric_limits<DistanceType>::max)()), pred(
                            2 * boost::num_vertices(g), std::make_tuple(std::make_pair(Vertex(), true), false, Edge())), index_in_heap_map(
                            parmcb::detail::SignedIndexInHeapFunctor<Graph>(boost::num_vertices(g), index_in_heap,
                                    index_map)), dist_map(
                            parmcb::detail::SignedDistanceFunctor<Graph, WeightMap>(boost::num_vertices(g), dist,
                                    index_map)), pred_map(
                            parmcb::detail::SignedPredecessorFunctor<Graph>(boost::num_vertices(g), pred, index_map)), compare(), queue(
                            dist_map, index_in_heap_map, compare) {
            }

            SignedVertex poll() {
                SignedVertex u = queue.top();
                queue.pop();
                return u;
            }

            const DistanceType& find_min() {
                SignedVertex u = queue.top();
                return boost::get(dist_map, u);
            }

            bool has_finite_dist(const SignedVertex &u) {
                return u == source || std::get<1>(boost::get(pred_map, u));
            }

            const DistanceType& get_dist(const SignedVertex &u) {
                return boost::get(dist_map, u);
            }

            const Predecessor& get_pred(const SignedVertex &u) {
                return boost::get(pred_map, u);
            }

            void push_source(SignedVertex s) {
                boost::put(dist_map, s, DistanceType());
                queue.push(s);
                source = s;
            }

            SignedVertex get_source() {
                return source;
            }

            void update(SignedVertex w, const DistanceType &c, const SignedVertex &pred, const Edge &pred_e) {
                if (w == source) {
                    return;
                }
                bool visited_w = std::get<1>(boost::get(pred_map, w));

                if (!visited_w) {
                    // first time found
                    boost::put(dist_map, w, c);
                    boost::put(pred_map, w, std::make_tuple(pred, true, pred_e));
                    queue.push(w);
                } else if (compare(c, boost::get(dist_map, w))) {
                    // already reached
                    boost::put(dist_map, w, c);
                    boost::put(pred_map, w, std::make_tuple(pred, true, pred_e));
                    queue.update(w);
                }
            }

        };

    } // detail

    template<class Graph, class WeightMap, class SignedEdges, class HiddenEdges>
    std::tuple<std::set<typename boost::graph_traits<Graph>::edge_descriptor>,
            typename boost::property_traits<WeightMap>::value_type, bool> signed_dijkstra(const Graph &g,
            const WeightMap &weight_map, const SignedEdges &signed_edges, const HiddenEdges &hidden_edges,
            bool use_hidden_edges, const typename boost::graph_traits<Graph>::vertex_descriptor &s, bool s_pos,
            const typename boost::graph_traits<Graph>::vertex_descriptor &t, bool t_pos, bool use_cycle_weight_limit,
            const typename boost::property_traits<WeightMap>::value_type &cycle_weight_limit) {

        typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
        typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
        typedef typename boost::property_traits<WeightMap>::value_type WeightType;
        typedef typename boost::property_traits<WeightMap>::value_type DistanceType;
        typedef std::pair<Vertex, bool> SignedVertex;
        typedef std::tuple<SignedVertex, bool, Edge> Predecessor;

        DistanceType distance_inf = (std::numeric_limits<DistanceType>::max)();
        std::less<DistanceType> compare;
        parmcb::detail::closed_plus<DistanceType> combine = parmcb::detail::closed_plus<DistanceType>();

        parmcb::detail::search_frontier<Graph, WeightMap> frontier(g);
        SignedVertex signed_s = std::make_pair(s, s_pos);
        SignedVertex signed_t = std::make_pair(t, t_pos);
        frontier.push_source(signed_s);

        assert(signed_s != signed_t);

        while (!frontier.queue.empty()) {
            SignedVertex signed_u = frontier.poll();
            DistanceType d_u = frontier.get_dist(signed_u);
            auto u = signed_u.first;

            if (use_cycle_weight_limit && !compare(d_u, cycle_weight_limit)) {
                // reached limit
                return std::make_tuple(std::set<Edge> { }, distance_inf, false);
            }

            if (signed_u == signed_t) { // found target
                std::set<Edge> cycle;
                WeightType cycle_weight = WeightType();
                SignedVertex signed_cur = signed_t;
                while (signed_cur != signed_s) {
                    Predecessor p = frontier.get_pred(signed_cur);
                    Edge e = std::get<2>(p);
                    if (!cycle.insert(e).second) {
                        // duplicate edge, discard cycle
                        return std::make_pair(std::set<Edge> { }, distance_inf);
                    } else {
                        cycle_weight += get(weight_map, e);
                    }
                    signed_cur = std::get<0>(p);
                }
                return std::make_tuple(cycle, cycle_weight, true);
            }

            auto eiRange = boost::out_edges(signed_u.first, g);
            for (auto ei = eiRange.first; ei != eiRange.second; ++ei) {
                auto e = *ei;
                if (use_hidden_edges && hidden_edges.find(e) != hidden_edges.end()) {
                    continue;
                }

                auto w = boost::target(e, g);
                if (w == u) {
                    w = boost::source(e, g);
                }
                if (w == u) {
                    // self-loop
                    continue;
                }
                const WeightType c = combine(d_u, get(weight_map, e));

                if (use_cycle_weight_limit && !compare(c, cycle_weight_limit)) {
                    // never insert if more than current minimum
                    continue;
                }

                bool is_signed = (signed_edges.find(e) != signed_edges.end());
                SignedVertex signed_w = std::make_pair(w, is_signed ? (!signed_u.second) : signed_u.second);

                frontier.update(signed_w, c, signed_u, e);
            }

        }

        return std::make_tuple(std::set<Edge> { }, distance_inf, false);
    }

    template<class Graph, class WeightMap, class SignedEdges, class HiddenEdges>
    std::tuple<std::set<typename boost::graph_traits<Graph>::edge_descriptor>,
            typename boost::property_traits<WeightMap>::value_type, bool> bidirectional_signed_dijkstra(const Graph &g,
            const WeightMap &weight_map, const SignedEdges &signed_edges, const HiddenEdges &hidden_edges,
            bool use_hidden_edges, const typename boost::graph_traits<Graph>::vertex_descriptor &s, bool s_pos,
            const typename boost::graph_traits<Graph>::vertex_descriptor &t, bool t_pos, bool use_cycle_weight_limit,
            const typename boost::property_traits<WeightMap>::value_type &cycle_weight_limit) {

        typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
        typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
        typedef typename boost::property_traits<WeightMap>::value_type WeightType;
        typedef typename boost::property_traits<WeightMap>::value_type DistanceType;
        typedef std::pair<Vertex, bool> SignedVertex;
        typedef std::tuple<SignedVertex, bool, Edge> Predecessor;

        DistanceType distance_inf = (std::numeric_limits<DistanceType>::max)();
        std::less<DistanceType> compare;
        parmcb::detail::closed_plus<DistanceType> combine = parmcb::detail::closed_plus<DistanceType>();

        SignedVertex signed_s = std::make_pair(s, s_pos);
        parmcb::detail::search_frontier<Graph, WeightMap> f_frontier(g);
        f_frontier.push_source(signed_s);

        SignedVertex signed_t = std::make_pair(t, t_pos);
        parmcb::detail::search_frontier<Graph, WeightMap> b_frontier(g);
        b_frontier.push_source(signed_t);

        assert(signed_s != signed_t);

        std::reference_wrapper<parmcb::detail::search_frontier<Graph, WeightMap>> frontier = std::ref(f_frontier);
        std::reference_wrapper<parmcb::detail::search_frontier<Graph, WeightMap>> other_frontier = std::ref(b_frontier);
        DistanceType best_path = distance_inf;
        bool best_path_set = false;
        SignedVertex best_path_common_vertex;

        while (true) {
            // stopping condition
            if (frontier.get().queue.empty() || other_frontier.get().queue.empty()
                    || (best_path_set
                            && !compare(combine(frontier.get().find_min(), other_frontier.get().find_min()), best_path))) {
                break;
            }

            // frontier scan
            SignedVertex signed_u = frontier.get().poll();
            DistanceType d_u = frontier.get().get_dist(signed_u);
            auto u = signed_u.first;

            if (use_cycle_weight_limit && !compare(d_u, cycle_weight_limit)) {
                // reached limit
                return std::make_tuple(std::set<Edge> { }, distance_inf, false);
            }

            auto eiRange = boost::out_edges(signed_u.first, g);
            for (auto ei = eiRange.first; ei != eiRange.second; ++ei) {
                auto e = *ei;
                if (use_hidden_edges && hidden_edges.find(e) != hidden_edges.end()) {
                    continue;
                }
                auto w = boost::target(e, g);
                if (w == u) {
                    w = boost::source(e, g);
                }
                if (w == u) {
                    // self-loop
                    continue;
                }

                const WeightType c = combine(d_u, get(weight_map, e));
                if (use_cycle_weight_limit && !frontier.get().compare(c, cycle_weight_limit)) {
                    // never insert if more than current minimum
                    continue;
                }

                bool is_signed = (signed_edges.find(e) != signed_edges.end());
                SignedVertex signed_w = std::make_pair(w, is_signed ? (!signed_u.second) : signed_u.second);

                frontier.get().update(signed_w, c, signed_u, e);

                if (other_frontier.get().has_finite_dist(signed_w)) {
                    // check path with w's distance from other frontier
                    DistanceType path_distance = combine(c, other_frontier.get().get_dist(signed_w));
                    if (compare(path_distance, best_path)) {
                        best_path_set = true;
                        best_path = path_distance;
                        best_path_common_vertex = signed_w;
                    }
                }

            }

            // swap frontiers
            std::swap(frontier, other_frontier);
        }

        if (!best_path_set || (use_cycle_weight_limit && !compare(best_path, cycle_weight_limit))) {
            return std::make_tuple(std::set<Edge> { }, distance_inf, false);
        }

        // create path if found
        std::set<Edge> cycle;
        WeightType cycle_weight = WeightType();

        SignedVertex signed_cur = best_path_common_vertex;
        SignedVertex signed_goal = frontier.get().get_source();
        while (signed_cur != signed_goal) {
            Predecessor p = frontier.get().get_pred(signed_cur);
            Edge e = std::get<2>(p);
            if (!cycle.insert(e).second) {
                // duplicate edge, discard cycle
                return std::make_tuple(std::set<Edge> { }, distance_inf, false);
            } else {
                cycle_weight += boost::get(weight_map, e);
            }
            signed_cur = std::get<0>(p);
        }

        signed_cur = best_path_common_vertex;
        signed_goal = other_frontier.get().get_source();
        while (signed_cur != signed_goal) {
            Predecessor p = other_frontier.get().get_pred(signed_cur);
            Edge e = std::get<2>(p);
            if (!cycle.insert(e).second) {
                // duplicate edge, discard cycle
                return std::make_tuple(std::set<Edge> { }, distance_inf, false);
            } else {
                cycle_weight += boost::get(weight_map, e);
            }
            signed_cur = std::get<0>(p);
        }
        return std::make_tuple(cycle, cycle_weight, true);
    }

} // parmcb

#endif
