#ifndef PARMCB_SVA_SIGNED_MPI_HPP_
#define PARMCB_SVA_SIGNED_MPI_HPP_

//    Copyright (C) Dimitrios Michail 2019 - 2020.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

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

#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/mpi/timer.hpp>

#include <tbb/tbb.h>
#include <tbb/concurrent_vector.h>

#include <parmcb/detail/signed_dijkstra.hpp>
#include <parmcb/mpi/sptrees.hpp>
#include <parmcb/forestindex.hpp>
#include <parmcb/spvecgf2.hpp>
#include <parmcb/util.hpp>

namespace parmcb {

    namespace detail {

        template<class Graph, class WeightMap>
        std::tuple<std::set<typename boost::graph_traits<Graph>::edge_descriptor>,
                typename boost::property_traits<WeightMap>::value_type, bool> find_shortest_odd_cycle_mpi(
                const Graph &g, const WeightMap &weight_map,
                const std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> &allVertices,
                const ForestIndex<Graph> &forest_index,
                const std::set<typename boost::graph_traits<Graph>::edge_descriptor> &signed_edges,
                boost::mpi::communicator &world) {

            typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
            typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
            typedef typename boost::property_traits<WeightMap>::value_type WeightType;

            std::less<WeightType> compare = std::less<WeightType>();
            std::tuple<std::set<Edge>, WeightType, bool> best = std::make_tuple(std::set<Edge> { },
                    (std::numeric_limits<WeightType>::max)(), false);
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

            if (signed_edges.size() == 1) {
                if (world.rank() == 0) {
                    auto se = *signed_edges.begin();
                    auto se_v = boost::source(se, g);
                    auto se_u = boost::target(se, g);
                    auto res = bidirectional_signed_dijkstra(g, weight_map, std::set<Edge> { }, signed_edges, true,
                            se_v, true, se_u, true, std::get<2>(best), std::get<1>(best));
                    if (std::get<2>(res) && std::get<0>(res).find(se) == std::get<0>(res).end()) {
                        std::get<1>(res) += boost::get(weight_map, se);
                        if (!std::get<2>(best) || compare(std::get<1>(res), std::get<1>(best))) {
                            std::get<0>(res).insert(se);
                            best = res;
                            assert(std::get<2>(best));
                        }
                    }
                }
            } else if (signed_edges.size() < boost::num_vertices(g)) {
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

                std::vector<Edge> local_signed_edges_as_vector;
                std::size_t total = signed_edges_as_vector.size();
                std::size_t stride = ceil((double) total / world.size());
                std::size_t istart = world.rank() * stride;
                std::size_t iend = istart + stride;
                for (std::size_t i = istart; i < iend && i < total; i++) {
                    local_signed_edges_as_vector.push_back(signed_edges_as_vector[i]);
                }

                std::tuple<std::set<Edge>, WeightType, bool> best_local_cycle = tbb::parallel_reduce(
                        tbb::blocked_range<std::size_t>(0, local_signed_edges_as_vector.size()),
                        std::make_tuple(std::set<Edge>(), (std::numeric_limits<WeightType>::max)(), false),
                        [&](tbb::blocked_range<std::size_t> r, auto running_min) {
                            for (std::size_t i = r.begin(); i < r.end(); i++) {
                                auto se = local_signed_edges_as_vector.at(i);
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

                std::vector<typename ForestIndex<Graph>::size_type> best_local_cycle_as_indices;
                convert_edges(std::get<0>(best_local_cycle),
                        std::inserter(best_local_cycle_as_indices, best_local_cycle_as_indices.end()), forest_index);
                SerializableMinOddCycle<Graph, WeightMap> local_min_odd_cycle(best_local_cycle_as_indices,
                        std::get<1>(best_local_cycle), std::get<2>(best_local_cycle));
                SerializableMinOddCycle<Graph, WeightMap> global_min_odd_cycle;

                boost::mpi::reduce(world, local_min_odd_cycle, global_min_odd_cycle,
                        SerializableMinOddCycleMinOp<Graph, WeightMap>(), 0);

                convert_edges(global_min_odd_cycle.edges, std::inserter(std::get<0>(best), std::get<0>(best).end()),
                        forest_index);
                std::get<1>(best) = global_min_odd_cycle.weight;
                std::get<2>(best) = global_min_odd_cycle.exists;
            } else {
                // split implicitly all vertices
                std::vector<Vertex> localVertices;
                std::size_t stride = ceil((double) allVertices.size() / world.size());
                std::size_t istart = world.rank() * stride;
                std::size_t iend = istart + stride;
                std::size_t total = allVertices.size();
                for (std::size_t i = istart; i < iend && i < total; i++) {
                    localVertices.push_back(allVertices[i]);
                }

                std::tuple<std::set<Edge>, WeightType, bool> best_local_cycle = tbb::parallel_reduce(
                        tbb::blocked_range<std::size_t>(0, localVertices.size()),
                        std::make_tuple(std::set<Edge>(), (std::numeric_limits<WeightType>::max)(), false),
                        [&](tbb::blocked_range<std::size_t> r, auto running_min) {
                            for (std::size_t i = r.begin(); i < r.end(); i++) {
                                auto v = localVertices[i];
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

                std::vector<typename ForestIndex<Graph>::size_type> best_local_cycle_as_indices;
                convert_edges(std::get<0>(best_local_cycle),
                        std::inserter(best_local_cycle_as_indices, best_local_cycle_as_indices.end()), forest_index);
                SerializableMinOddCycle<Graph, WeightMap> local_min_odd_cycle(best_local_cycle_as_indices,
                        std::get<1>(best_local_cycle), std::get<2>(best_local_cycle));
                SerializableMinOddCycle<Graph, WeightMap> global_min_odd_cycle;

                boost::mpi::reduce(world, local_min_odd_cycle, global_min_odd_cycle,
                        SerializableMinOddCycleMinOp<Graph, WeightMap>(), 0);

                convert_edges(global_min_odd_cycle.edges, std::inserter(std::get<0>(best), std::get<0>(best).end()),
                        forest_index);
                std::get<1>(best) = global_min_odd_cycle.weight;
                std::get<2>(best) = global_min_odd_cycle.exists;
            }

            return best;
        }

    }
// detail

    template<class Graph, class WeightMap, class CycleOutputIterator>
    typename boost::property_traits<WeightMap>::value_type mcb_sva_signed_mpi(const Graph &g, WeightMap weight_map,
            CycleOutputIterator out, boost::mpi::communicator &world, const std::size_t hardware_concurrency_hint = 0) {

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

        boost::mpi::timer total_timer;

        /*
         * Main loop
         */
        WeightType mcb_weight = WeightType();
        for (std::size_t k = 0; k < csd; k++) {
            if (k % 250 == 0) {
                std::cout << "Rank " << world.rank() << " at cycle " << k << std::endl;
            }

            // TODO: check if sparsest support heuristic makes sense here

            // broadcast support vector
            if (world.rank() == 0) {
                boost::mpi::broadcast(world, support[k], 0);
            } else {
                SpVecGF2<std::size_t> received;
                boost::mpi::broadcast(world, received, 0);
                support[k] = received;
            }

            std::set<Edge> signed_edges;
            convert_edges(support[k], std::inserter(signed_edges, signed_edges.end()), forest_index);
            std::tuple<std::set<Edge>, WeightType, bool> best = parmcb::detail::find_shortest_odd_cycle_mpi(g, weight_map,
                    vertices, forest_index, signed_edges, world);

            if (world.rank() == 0) {
                /*
                 * Update support vectors
                 */
                std::set<std::size_t> cyclek;
                convert_edges(std::get<0>(best), std::inserter(cyclek, cyclek.end()), forest_index);
                tbb::parallel_for(tbb::blocked_range<std::size_t>(k + 1, csd),
                        [&](const tbb::blocked_range<std::size_t> &r) {
                            auto e = r.end();
                            for (std::size_t i = r.begin(); i != e; ++i) {
                                if (support[i] * cyclek == 1) {
                                    support[i] += support[k];
                                }
                            }
                        });

                /*
                 * Output cycles
                 */
                std::list<Edge> cyclek_edgelist;
                std::copy(std::get<0>(best).begin(), std::get<0>(best).end(), std::back_inserter(cyclek_edgelist));
                *out++ = cyclek_edgelist;
                mcb_weight += std::get<1>(best);
            }

        }

        if (world.rank() == 0) {
            std::cout << "Total time: " << total_timer.elapsed() << " (sec)" << std::endl;
        }

        return mcb_weight;
    }

} // namespace parmcb

#endif
