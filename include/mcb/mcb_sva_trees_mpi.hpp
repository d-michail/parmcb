#ifndef LIBMCB_MCB_SVA_TREES_MPI_HPP_
#define LIBMCB_MCB_SVA_TREES_MPI_HPP_

#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/tuple/detail/tuple_basic.hpp>
#include <boost/timer/timer.hpp>

#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>

#include <cstddef>
#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <set>
#include <vector>
#include <cmath>

#include <mcb/forestindex.hpp>
#include <mcb/spvecgf2.hpp>
#include <mcb/fvs.hpp>
#include <mcb/sptrees.hpp>
#include <mcb/util.hpp>

namespace mcb {

    template<class Graph, class WeightMap, class CycleOutputIterator>
    typename boost::property_traits<WeightMap>::value_type mcb_sva_trees_mpi(const Graph &g, WeightMap weight_map,
            boost::mpi::communicator &world, CycleOutputIterator out) {

        typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
        typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
        typedef typename boost::property_traits<WeightMap>::value_type WeightType;

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
            std::vector<Vertex> feedback_vertex_set;
            mcb::greedy_fvs(g, std::back_inserter(feedback_vertex_set));
            for (auto x : feedback_vertex_set) {
                std::cout << "Feedback vertex set: " << x << std::endl;
            }
            std::vector<SerializableCandidateCycle<Graph>> all_candidate_cycles;
            for (const auto &v : feedback_vertex_set) {
                const std::size_t tree_id = 0;
                mcb::SPTree<Graph, WeightMap> tree(tree_id, g, weight_map, v);
                std::vector<SerializableCandidateCycle<Graph>> tree_candidate_cycles =
                        tree.create_serializable_candidate_cycles(forest_index);
                std::ostringstream stream;
                for (auto x : tree_candidate_cycles) {
                    stream << "(" << x.v << "," << forest_index(x.e) << ")";
                }
                std::cout << "For vertex " << v << " candidate cycles are " << stream.str() << std::endl;

                all_candidate_cycles.insert(all_candidate_cycles.end(), tree_candidate_cycles.begin(),
                        tree_candidate_cycles.end());
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

        std::ostringstream stream;
        for (auto x : candidate_cycles) {
            stream << "(" << x.v << "," << forest_index(x.e) << ")";
        }
        std::cout << "Rank " << world.rank() << " received " << candidate_cycles.size() << " candidates: "
                << stream.str() << std::endl;

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

        std::cout << "Rank " << world.rank() << " grouped per source vertex" << std::endl;

        /*
         * Build trees for local candidate cycles
         */
        const bool sorted_cycles = false;
        SPTrees<Graph, WeightMap, false> sp_trees(g, weight_map, perVertexCandidates, sorted_cycles);

        /*
         * Synchronize
         */
        world.barrier();

        /*
         * Main loop
         */
        WeightType mcb_weight = WeightType();
        for (std::size_t k = 0; k < csd; k++) {
            //if (k % 250 == 0) {
            std::cout << "Rank " << world.rank() << " at cycle " << k << std::endl;
            //}

            // broadcast support vector
            if (world.rank() == 0) {
                boost::mpi::broadcast(world, support[k], 0);
            } else {
                SpVecGF2<std::size_t> received;
                boost::mpi::broadcast(world, received, 0);
                support[k] = received;
            }

            std::cout << "Rank " << world.rank() << " has support = " << support[k] << std::endl;

            std::set<Edge> signed_edges;
            convert_edges(support[k], std::inserter(signed_edges, signed_edges.end()), forest_index);

            std::tuple<std::set<Edge>, WeightType, bool> best_local_cycle = sp_trees.compute_shortest_odd_cycle(
                    signed_edges);

            if (std::get<2>(best_local_cycle)) {
                std::cout << "Rank " << world.rank() << " found best local cycle with weight = " << std::get<1>(best_local_cycle) << std::endl;
            } else {
                std::cout << "Rank " << world.rank() << " found no local" << std::endl;
            }

            world.barrier();

            std::vector<typename ForestIndex<Graph>::size_type> best_local_cycle_as_indices;
            convert_edges(std::get<0>(best_local_cycle),
                    std::inserter(best_local_cycle_as_indices, best_local_cycle_as_indices.end()), forest_index);
            SerializableMinOddCycle<Graph, WeightMap> local_min_odd_cycle(best_local_cycle_as_indices,
                    std::get<1>(best_local_cycle), std::get<2>(best_local_cycle));
            SerializableMinOddCycle<Graph, WeightMap> global_min_odd_cycle;

            // TODO: use is commutative!!
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

            world.barrier();
        }

        return mcb_weight;
    }

} // namespace mcb

#endif
