#ifndef LIBMCB_SVA_SIGNED_PARALLEL2_HPP_
#define LIBMCB_SVA_SIGNED_PARALLEL2_HPP_

#include <cstddef>
#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <set>
#include <vector>

#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/tuple/detail/tuple_basic.hpp>
#include <boost/timer/timer.hpp>

#include <tbb/tbb.h>
#include <tbb/concurrent_vector.h>

#include <mcb/forestindex.hpp>
#include <mcb/spvecgf2.hpp>
#include <mcb/signed_dijkstra.hpp>
#include <mcb/util.hpp>

namespace mcb {

    template<class Graph, class WeightMap, class CycleOutputIterator>
    typename boost::property_traits<WeightMap>::value_type mcb_sva_parallel2(const Graph &g, WeightMap weight_map,
            CycleOutputIterator out) {

        typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
        typedef typename boost::graph_traits<Graph>::vertex_iterator VertexIt;
        typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
        typedef typename boost::property_traits<WeightMap>::value_type WeightType;

        /*
         * Index the graph
         */
        ForestIndex<Graph> forest_index(g);
        auto csd = forest_index.cycle_space_dimension();
        std::vector<Vertex> vertices;
        {
            VertexIt vi, viend;
            for (boost::tie(vi, viend) = boost::vertices(g); vi != viend; ++vi) {
                vertices.push_back(*vi);
            }
        }

        /*
         * Initialize support vectors
         */
        tbb::concurrent_vector<SpVecGF2<std::size_t>> support;
        tbb::parallel_for(tbb::blocked_range<std::size_t>(0, csd), [&](const tbb::blocked_range<std::size_t> &r) {
            for (std::size_t i = r.begin(); i != r.end(); ++i) {
                support.push_back(SpVecGF2<std::size_t> { i });
            }
        });

        boost::timer::cpu_timer cycle_timer;
        cycle_timer.stop();
        boost::timer::cpu_timer support_timer;
        support_timer.stop();

        /*
         * Main loop
         */
        WeightType mcb_weight = WeightType();
        for (std::size_t k = 0; k < csd; k++) {
            if (k % 250 == 0) {
                std::cout << k << std::endl;
            }

            /*
             * Choose the sparsest support heuristic
             */
            auto min_support = k;
            for (auto r = k + 1; r < csd; ++r) {
                if (support[r].size() < support[min_support].size())
                    min_support = r;
            }
            if (min_support != k) {  // swap
                std::swap(support[k], support[min_support]);
            }

            /**
             * TODO:
             * As long as there are a lot of remaining support vectors, find two or more that
             * are disjoint and find cycles in parallel.
             * Afterwards validate that we were lucky.
             *
             * (2 cycles case) Assume we compute shortest cycles such that <C1,S1>=1 and <C2,S2>=1.
             * If we are lucky and <C1,S2>=0, we are done. Otherwise, we need to update S1=S1+S2
             * and throw away C2.
             *
             * (3 cycles case) TODO
             */

            /*
             * Compute shortest odd cycle
             */
            cycle_timer.resume();
            std::set<Edge> signed_edges;
            convert_edges(support[k], std::inserter(signed_edges, signed_edges.end()), forest_index);
            std::less<WeightType> compare = std::less<WeightType>();
            WeightType weight_inf = (std::numeric_limits<WeightType>::max)();
            std::set<Edge> best_cycle;
            WeightType best_cycle_weight = weight_inf;

            typedef std::pair<std::set<Edge>, WeightType> cycle_t;
            auto cycle_min = [compare](const cycle_t &c1, const cycle_t &c2) {
                if (!compare(c2.second, c1.second)) {
                    return c1;
                }
                return c2;
            };

            if (signed_edges.size() >= boost::num_vertices(g)) {
                boost::tie(best_cycle, best_cycle_weight) = tbb::parallel_reduce(
                        tbb::blocked_range<std::size_t>(0, boost::num_vertices(g)),
                        std::make_pair(std::set<Edge>(), (std::numeric_limits<WeightType>::max)()),
                        [&](tbb::blocked_range<std::size_t> r, auto running_min) {
                            for (std::size_t i = r.begin(); i < r.end(); i++) {
                                auto v = vertices[i];
                                const bool use_hidden_edges = false;
                                auto res = bidirectional_signed_dijkstra(g, weight_map, signed_edges,
                                        std::set<Edge> { }, use_hidden_edges, v, true, v, false, running_min.second);
                                if (!res.first.empty() && compare(res.second, weight_inf)
                                        && compare(res.second, running_min.second)) {
                                    running_min.first = res.first;
                                    running_min.second = res.second;
                                }
                            }
                            return running_min;
                        },
                        cycle_min);
            } else if (signed_edges.size() == 1) {
                /*
                 * Heuristic in case number of signed edges is small compared to the number of vertices.
                 */
                std::set<Edge> hidden_edges;
                std::copy(signed_edges.begin(), signed_edges.end(), std::inserter(hidden_edges, hidden_edges.begin()));
                for (auto sei = signed_edges.begin(); sei != signed_edges.end(); ++sei) {
                    auto se = *sei;
                    auto se_v = boost::source(se, g);
                    auto se_u = boost::target(se, g);
                    auto res = bidirectional_signed_dijkstra(g, weight_map, signed_edges, hidden_edges, true, se_v,
                            true, se_u, true, best_cycle_weight);
                    hidden_edges.erase(hidden_edges.begin());
                    if (!res.first.empty() && compare(res.second, weight_inf)
                            && res.first.find(se) == res.first.end()) {
                        res.second += boost::get(weight_map, se);
                        if (compare(res.second, best_cycle_weight)) {
                            res.first.insert(se);
                            best_cycle = res.first;
                            best_cycle_weight = res.second;
                        }
                    }
                }
            } else {
                /*
                 * Heuristic in case number of signed edges is small compared to the number of vertices.
                 */
                std::map<Edge, std::set<Edge>> hidden_edges_per_edge;
                std::vector<Edge> signed_edges_as_vector;
                std::set<Edge> tmp_signed_edges = signed_edges;
                while (!tmp_signed_edges.empty()) {
                    auto bit = tmp_signed_edges.begin();
                    hidden_edges_per_edge.insert(std::make_pair(*bit, tmp_signed_edges));
                    signed_edges_as_vector.push_back(*bit);
                    tmp_signed_edges.erase(bit);
                }
                boost::tie(best_cycle, best_cycle_weight) = tbb::parallel_reduce(
                        tbb::blocked_range<std::size_t>(0, signed_edges_as_vector.size()),
                        std::make_pair(std::set<Edge>(), (std::numeric_limits<WeightType>::max)()),
                        [&](tbb::blocked_range<std::size_t> r, auto running_min) {
                            for (std::size_t i = r.begin(); i < r.end(); i++) {
                                auto se = signed_edges_as_vector.at(i);
                                auto se_v = boost::source(se, g);
                                auto se_u = boost::target(se, g);
                                auto hidden_edges = hidden_edges_per_edge.at(se);
                                auto res = bidirectional_signed_dijkstra(g, weight_map, signed_edges, hidden_edges,
                                        true, se_v, true, se_u, true, running_min.second);
                                if (!res.first.empty() && compare(res.second, weight_inf)
                                        && res.first.find(se) == res.first.end()) {
                                    res.second += boost::get(weight_map, se);
                                    if (compare(res.second, running_min.second)) {
                                        res.first.insert(se);
                                        running_min.first = res.first;
                                        running_min.second = res.second;
                                    }
                                }
                            }
                            return running_min;
                        },
                        cycle_min);
            }
            cycle_timer.stop();

            /*
             * Update support vectors
             */
            support_timer.resume();
            std::set<std::size_t> cyclek;
            convert_edges(best_cycle, std::inserter(cyclek, cyclek.end()), forest_index);
            tbb::parallel_for(tbb::blocked_range<std::size_t>(k + 1, csd),
                    [&](const tbb::blocked_range<std::size_t> &r) {
                        auto e = r.end();
                        for (std::size_t i = r.begin(); i != e; ++i) {
                            if (support[i] * cyclek == 1) {
                                support[i] += support[k];
                            }
                        }
                    });
            support_timer.stop();

            /*
             * Output new cycle
             */
            std::list<Edge> cyclek_edgelist;
            std::copy(best_cycle.begin(), best_cycle.end(), std::back_inserter(cyclek_edgelist));
            *out++ = cyclek_edgelist;
            mcb_weight += best_cycle_weight;
        }

        std::cout << "cycle   timer" << cycle_timer.format();
        std::cout << "support timer" << support_timer.format();

        return mcb_weight;
    }

} // namespace mcb

#endif
