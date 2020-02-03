#ifndef MCB_SVA_SIGNED_TBB_HPP_
#define MCB_SVA_SIGNED_TBB_HPP_

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
#include <tbb/task_group.h>

#include <mcb/detail/signed_dijkstra.hpp>
#include <mcb/forestindex.hpp>
#include <mcb/spvecgf2.hpp>
#include <mcb/util.hpp>

namespace mcb {

    namespace detail {

        template<class Graph, class WeightMap>
        struct OddCycleFinder {
            typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
            typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
            typedef typename boost::property_traits<WeightMap>::value_type WeightType;

            OddCycleFinder(const Graph &g, const WeightMap &weight_map, const ForestIndex<Graph> &forest_index,
                    const std::vector<Vertex> &vertices) :
                    g(g), weight_map(weight_map), forest_index(forest_index), vertices(vertices), compare(
                            std::less<WeightType>()) {
            }

            std::tuple<std::set<Edge>, WeightType, bool> find(const SpVecGF2<std::size_t> &support) {
                std::set<Edge> signed_edges;
                convert_edges(support, std::inserter(signed_edges, signed_edges.end()), forest_index);
                if (signed_edges.size() == 1) {
                    return find_single_edge(*signed_edges.begin());
                } else if (signed_edges.size() >= boost::num_vertices(g)) {
                    return find_all_vertices(signed_edges);
                } else {
                    return find_less_than_vertices(signed_edges);
                }
            }

            std::tuple<std::set<Edge>, WeightType, bool> find_single_edge(const Edge &se) {
                std::tuple<std::set<Edge>, WeightType, bool> best;

                auto se_v = boost::source(se, g);
                auto se_u = boost::target(se, g);
                auto res = bidirectional_signed_dijkstra(g, weight_map, std::set<Edge> { }, std::set<Edge> { se }, true,
                        se_v, true, se_u, true, std::get<2>(best), std::get<1>(best));
                if (std::get<2>(res) && std::get<0>(res).find(se) == std::get<0>(res).end()) {
                    std::get<1>(res) += boost::get(weight_map, se);
                    if (!std::get<2>(best) || compare(std::get<1>(res), std::get<1>(best))) {
                        std::get<0>(res).insert(se);
                        best = res;
                    }
                }
                return best;
            }

            std::tuple<std::set<Edge>, WeightType, bool> find_all_vertices(const std::set<Edge> &signed_edges) {
                typedef std::tuple<std::set<Edge>, WeightType, bool> cycle_t;
                auto cycle_min = [&](const cycle_t &c1, const cycle_t &c2) {
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

                return tbb::parallel_reduce(tbb::blocked_range<std::size_t>(0, boost::num_vertices(g)),
                        std::make_tuple(std::set<Edge>(), (std::numeric_limits<WeightType>::max)(), false),
                        [&](tbb::blocked_range<std::size_t> r, auto running_min) {
                            for (std::size_t i = r.begin(); i < r.end(); i++) {
                                auto v = vertices[i];
                                const bool use_hidden_edges = false;
                                auto res = bidirectional_signed_dijkstra(g, weight_map, signed_edges,
                                        std::set<Edge> { }, use_hidden_edges, v, true, v, false,
                                        std::get<2>(running_min), std::get<1>(running_min));
                                if (std::get<2>(res)
                                        && (!std::get<2>(running_min)
                                                || compare(std::get<1>(res), std::get<1>(running_min)))) {
                                    running_min = res;
                                }
                            }
                            return running_min;
                        },
                        cycle_min);

            }

            std::tuple<std::set<Edge>, WeightType, bool> find_less_than_vertices(const std::set<Edge> &signed_edges) {
                /*
                 * Heuristic in case number of signed edges is small compared to the number of vertices.
                 */
                typedef std::tuple<std::set<Edge>, WeightType, bool> cycle_t;
                auto cycle_min = [&](const cycle_t &c1, const cycle_t &c2) {
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

                std::map<Edge, std::set<Edge>> hidden_edges_per_edge;
                std::vector<Edge> signed_edges_as_vector;
                std::set<Edge> tmp_signed_edges = signed_edges;
                while (!tmp_signed_edges.empty()) {
                    auto bit = tmp_signed_edges.begin();
                    hidden_edges_per_edge.insert(std::make_pair(*bit, tmp_signed_edges));
                    signed_edges_as_vector.push_back(*bit);
                    tmp_signed_edges.erase(bit);
                }
                return tbb::parallel_reduce(tbb::blocked_range<std::size_t>(0, signed_edges_as_vector.size()),
                        std::make_tuple(std::set<Edge>(), (std::numeric_limits<WeightType>::max)(), false),
                        [&](tbb::blocked_range<std::size_t> r, auto running_min) {
                            for (std::size_t i = r.begin(); i < r.end(); i++) {
                                auto se = signed_edges_as_vector.at(i);
                                auto se_v = boost::source(se, g);
                                auto se_u = boost::target(se, g);
                                auto hidden_edges = hidden_edges_per_edge.at(se);
                                auto res = bidirectional_signed_dijkstra(g, weight_map, signed_edges, hidden_edges,
                                        true, se_v, true, se_u, true, std::get<2>(running_min),
                                        std::get<1>(running_min));
                                if (std::get<2>(res) && std::get<0>(res).find(se) == std::get<0>(res).end()) {
                                    std::get<1>(res) += boost::get(weight_map, se);
                                    if (!std::get<2>(running_min)
                                            || compare(std::get<1>(res), std::get<1>(running_min))) {
                                        std::get<0>(res).insert(se);
                                        running_min = res;
                                    }
                                }
                            }
                            return running_min;
                        },
                        cycle_min);
            }

            const Graph &g;
            const WeightMap &weight_map;
            const ForestIndex<Graph> &forest_index;
            const std::vector<Vertex> &vertices;
            const std::less<WeightType> compare;
        };

    }

    template<class Graph, class WeightMap, class CycleOutputIterator>
    typename boost::property_traits<WeightMap>::value_type mcb_sva_signed_tbb(const Graph &g, WeightMap weight_map,
            CycleOutputIterator out, const std::size_t hardware_concurrency_hint = 0) {

        typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
        typedef typename boost::graph_traits<Graph>::vertex_iterator VertexIt;
        typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
        typedef typename boost::property_traits<WeightMap>::value_type WeightType;

        /*
         * Index the graph
         */
        ForestIndex<Graph> forest_index(g);
        auto csd = forest_index.cycle_space_dimension();
        std::cout << "Cycle dimension " << csd << std::endl;
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
        mcb::detail::OddCycleFinder<Graph, WeightMap> odd_cycle_finder(g, weight_map, forest_index, vertices);
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

            /*
             * Compute shortest odd cycle
             */
            cycle_timer.resume();
            auto cycle = odd_cycle_finder.find(support[k]);
            cycle_timer.stop();

            /*
             * Update support vectors
             */
            support_timer.resume();
            std::set<std::size_t> cyclek;
            convert_edges(std::get<0>(cycle), std::inserter(cyclek, cyclek.end()), forest_index);
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
             * Report cycle
             */
            std::list<Edge> cyclek_edgelist;
            std::copy(std::get<0>(cycle).begin(), std::get<0>(cycle).end(), std::back_inserter(cyclek_edgelist));
            *out++ = cyclek_edgelist;
            mcb_weight += std::get<1>(cycle);
        }

        std::cout << "cycle   timer" << cycle_timer.format();
        std::cout << "support timer" << support_timer.format();

        return mcb_weight;
    }

} // namespace mcb

#endif
