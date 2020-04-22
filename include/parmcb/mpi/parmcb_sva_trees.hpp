#ifndef PARMCB_MPI_SVA_TREES_HPP_
#define PARMCB_MPI_SVA_TREES_HPP_

#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/tuple/detail/tuple_basic.hpp>

#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/mpi/timer.hpp>

#include <cstddef>
#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <set>
#include <vector>
#include <cmath>

#include <parmcb/forestindex.hpp>
#include <parmcb/spvecgf2.hpp>
#include <parmcb/detail/fvs.hpp>
#include <parmcb/detail/cycles.hpp>
#include <parmcb/util.hpp>
#include <parmcb/mpi/sptrees.hpp>

namespace parmcb {

    template<class Graph, class WeightMap, class CycleOutputIterator, class CyclesBuilder, bool ParallelUsingTBB>
    typename boost::property_traits<WeightMap>::value_type _mcb_sva_trees_mpi(const Graph &g, WeightMap weight_map,
            CycleOutputIterator out, boost::mpi::communicator &world) {

        typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
        typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
        typedef typename boost::property_traits<WeightMap>::value_type WeightType;

        boost::mpi::timer total_timer;

        /*
         * Index the graph
         */
        ForestIndex<Graph> forest_index(g);
        auto csd = forest_index.cycle_space_dimension();
        if (world.rank() == 0) {
            std::cout << "Cycle space dimension: " << csd << std::endl;
        }

        /*
         * Initialize support vectors
         */
        std::vector<SpVecGF2<std::size_t>> support;
        for (std::size_t k = 0; k < csd; k++) {
            support.emplace_back(k);
        }

        /**
         * Scatter candidate cycles
         */
        std::vector<SerializableCandidateCycle<Graph>> candidate_cycles;
        if (world.rank() == 0) {
            /*
             * Compute all vertex-edge candidate pairs
             */
            std::vector<parmcb::SPTree<Graph, WeightMap>> trees;
            std::vector<parmcb::CandidateCycle<Graph, WeightMap>> cycles;
            CyclesBuilder cycles_builder;
            cycles_builder(g, weight_map, trees, cycles);
            std::cout << "Total candidate cycles: " << cycles.size() << std::endl;

            std::vector<SerializableCandidateCycle<Graph>> all_candidate_cycles;
            parmcb::CandidateCycleToSerializableConverter<Graph, WeightMap> converter(trees, forest_index);
            for (const auto &c : cycles) {
                all_candidate_cycles.push_back(converter(c));
            }

            /*
             * Divide them to processes
             */
            std::size_t stride = ceil((double) all_candidate_cycles.size() / world.size());
            std::vector<std::vector<SerializableCandidateCycle<Graph>>> chunks;
            for (int p = 0; p < world.size(); p++) {
                std::size_t istart = p * stride;
                std::size_t iend = istart + stride;
                std::vector<SerializableCandidateCycle<Graph>> chunk;
                std::size_t total = all_candidate_cycles.size();
                for (std::size_t i = istart; i < iend && i < total; i++) {
                    chunk.push_back(all_candidate_cycles[i]);
                }
                chunks.push_back(chunk);
            }
            boost::mpi::scatter(world, chunks, candidate_cycles, 0);
        } else {
            boost::mpi::scatter(world, std::vector<std::vector<SerializableCandidateCycle<Graph>>> { },
                    candidate_cycles, 0);
        }

        std::cout << "Rank " << world.rank() << " received " << candidate_cycles.size() << " candidates" << std::endl;

        /*
         * Group local candidate cycles per vertex
         */
        std::map<Vertex, std::vector<Edge>> perVertexCandidates;
        for (SerializableCandidateCycle<Graph> t : candidate_cycles) {
            Vertex s = t.v;
            Edge e = forest_index(t.e);
            if (perVertexCandidates.find(s) == perVertexCandidates.end()) {
                perVertexCandidates.insert(std::make_pair(s, std::vector<Edge> { }));
            }
            perVertexCandidates[s].push_back(e);
        }

        /*
         * Build trees for local candidate cycles
         */
        std::vector<parmcb::SPTree<Graph, WeightMap>> trees;
        std::vector<parmcb::CandidateCycle<Graph, WeightMap>> cycles;
        for (auto const &p : perVertexCandidates) {
            SPTree<Graph, WeightMap> tree(trees.size(), g, weight_map, p.first);
            trees.push_back(tree);
            std::vector<CandidateCycle<Graph, WeightMap>> tree_cycles = tree.create_candidate_cycles(p.second.begin(),
                    p.second.end());
            cycles.insert(cycles.end(), tree_cycles.begin(), tree_cycles.end());
        }
        std::cout << "Total candidate cycles: " << cycles.size() << std::endl;
        const bool sorted_cycles = true;
        if (sorted_cycles) {
            // sort
            std::cout << "Sorting cycles" << std::endl;
            std::sort(cycles.begin(), cycles.end(), [](const auto &a, const auto &b) {
                return a.weight() < b.weight();
            });
        }
        ShortestOddCycleLookup<Graph, WeightMap, ParallelUsingTBB> cycle_lookup(g, weight_map, trees, cycles,
                sorted_cycles);

        /*
         * Main loop
         */
        WeightType mcb_weight = WeightType();
        for (std::size_t k = 0; k < csd; k++) {
            if (k % 250 == 0) {
                std::cout << "Rank " << world.rank() << " at cycle " << k << std::endl;
            }

            // broadcast support vector
            if (world.rank() == 0) {
                boost::mpi::broadcast(world, support[k], 0);
            } else {
                SpVecGF2<std::size_t> received;
                boost::mpi::broadcast(world, received, 0);
                support[k] = received;
            }

            //std::cout << "Rank " << world.rank() << " has support = " << support[k] << std::endl;

            std::set<Edge> signed_edges;
            convert_edges(support[k], std::inserter(signed_edges, signed_edges.end()), forest_index);

            std::tuple<std::set<Edge>, WeightType, bool> best_local_cycle = cycle_lookup(signed_edges);

            std::vector<typename ForestIndex<Graph>::size_type> best_local_cycle_as_indices;
            convert_edges(std::get<0>(best_local_cycle),
                    std::inserter(best_local_cycle_as_indices, best_local_cycle_as_indices.end()), forest_index);
            SerializableMinOddCycle<Graph, WeightMap> local_min_odd_cycle(best_local_cycle_as_indices,
                    std::get<1>(best_local_cycle), std::get<2>(best_local_cycle));
            SerializableMinOddCycle<Graph, WeightMap> global_min_odd_cycle;

            boost::mpi::reduce(world, local_min_odd_cycle, global_min_odd_cycle,
                    SerializableMinOddCycleMinOp<Graph, WeightMap>(), 0);

            if (world.rank() == 0) {
                std::set<std::size_t> cyclek;
                for (std::size_t e : global_min_odd_cycle.edges) {
                    cyclek.insert(e);
                }
                for (std::size_t l = k + 1; l < csd; l++) {
                    if (support[l] * cyclek == 1) {
                        support[l] += support[k];
                    }
                }

                std::list<Edge> cyclek_edgelist;
                convert_edges(global_min_odd_cycle.edges, std::inserter(cyclek_edgelist, cyclek_edgelist.end()),
                        forest_index);
                mcb_weight += global_min_odd_cycle.weight;
                *out++ = cyclek_edgelist;
            }
        }

        if (world.rank() == 0) {
            std::cout << "Total time: " << total_timer.elapsed() << " (sec)" << std::endl;
        }

        return mcb_weight;

    }

    template<class Graph, class WeightMap, class CycleOutputIterator>
    typename boost::property_traits<WeightMap>::value_type mcb_sva_fvs_trees_mpi(const Graph &g, WeightMap weight_map,
            CycleOutputIterator out, boost::mpi::communicator &world) {
        return _mcb_sva_trees_mpi<Graph, WeightMap, CycleOutputIterator,
                parmcb::detail::FVSCyclesBuilder<Graph, WeightMap>, false>(g, weight_map, out, world);
    }

    template<class Graph, class WeightMap, class CycleOutputIterator>
    typename boost::property_traits<WeightMap>::value_type mcb_sva_fvs_trees_tbb_mpi(const Graph &g,
            WeightMap weight_map, CycleOutputIterator out, boost::mpi::communicator &world) {
        return _mcb_sva_trees_mpi<Graph, WeightMap, CycleOutputIterator,
                parmcb::detail::FVSCyclesBuilder<Graph, WeightMap>, true>(g, weight_map, out, world);
    }

    template<class Graph, class WeightMap, class CycleOutputIterator>
    typename boost::property_traits<WeightMap>::value_type mcb_sva_iso_trees_mpi(const Graph &g, WeightMap weight_map,
            CycleOutputIterator out, boost::mpi::communicator &world) {
        return _mcb_sva_trees_mpi<Graph, WeightMap, CycleOutputIterator,
                parmcb::detail::ISOCyclesBuilder<Graph, WeightMap>, false>(g, weight_map, out, world);
    }

    template<class Graph, class WeightMap, class CycleOutputIterator>
    typename boost::property_traits<WeightMap>::value_type mcb_sva_iso_trees_tbb_mpi(const Graph &g,
            WeightMap weight_map, CycleOutputIterator out, boost::mpi::communicator &world) {
        return _mcb_sva_trees_mpi<Graph, WeightMap, CycleOutputIterator,
                parmcb::detail::ISOCyclesBuilder<Graph, WeightMap>, true>(g, weight_map, out, world);
    }

} // namespace mcb

#endif
