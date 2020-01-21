#ifndef LIBMCB_MCB_SVA_TREES_MPI_HPP_
#define LIBMCB_MCB_SVA_TREES_MPI_HPP_

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/lookup_edge.hpp>
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
        if (world.rank() == 0) {
            for (std::size_t k = 0; k < csd; k++) {
                support.emplace_back(k);
            }
        }

        /**
         * Scatter candidate cycles
         */
        std::vector<SerializableVertexTriplet<Graph>> candidate_triplets;
        if (world.rank() == 0) {
            /*
             * Initialize support vectors
             */
            for (std::size_t k = 0; k < csd; k++) {
                support.emplace_back(k);
            }

            /*
             * Compute all vertex-edge candidate pairs
             */
            std::vector<Vertex> feedback_vertex_set;
            mcb::greedy_fvs(g, std::back_inserter(feedback_vertex_set));
            std::vector<SerializableVertexTriplet<Graph>> all_candidate_triplets;
            for (auto v : feedback_vertex_set) {
                mcb::SPTree<Graph, WeightMap> tree(g, weight_map, v);
                std::vector<SerializableVertexTriplet<Graph>> tree_candidate_triplets =
                        tree.create_candidate_vertex_triplets();
                all_candidate_triplets.insert(all_candidate_triplets.end(), tree_candidate_triplets.begin(),
                        tree_candidate_triplets.end());
            }

            /*
             * Divide them to processes
             */
            std::size_t stride = ceil((double) all_candidate_triplets.size() / world.size());
            std::vector<std::vector<SerializableVertexTriplet<Graph>>> chunks;
            for (std::size_t p = 0; p < world.size(); p++) {
                std::size_t istart = p * stride;
                std::size_t iend = istart + stride;
                std::vector<SerializableVertexTriplet<Graph>> chunk;
                std::size_t total = all_candidate_triplets.size();
                for (std::size_t i = istart; i < iend && i < total; i++) {
                    chunk.push_back(all_candidate_triplets[i]);
                }
                chunks.push_back(chunk);
            }
            boost::mpi::scatter(world, chunks, candidate_triplets, 0);
        } else {
            boost::mpi::scatter(world, std::vector<std::vector<SerializableVertexTriplet<Graph>>> { },
                    candidate_triplets, 0);
        }

        std::cout << "Rank " << world.rank() << " received " << candidate_triplets.size() << " candidates" << std::endl;

        /*
         * Group local candidate cycles per vertex
         */
        std::map<Vertex, std::vector<Edge>> perVertexCandidates;
        for (SerializableVertexTriplet<Graph> t : candidate_triplets) {
            Vertex s = t.v;
            auto ePair = boost::lookup_edge(t.u, t.w, g);
            if (!ePair.second) {
                throw std::system_error(EIO, std::generic_category(),
                        "Input error, the graph seems to contain multiple edges");
            }
            Edge e = ePair.first;

            if (perVertexCandidates.find(s) == perVertexCandidates.end()) {
                perVertexCandidates.insert(std::make_pair(s, std::vector<Edge> { }));
            }
            perVertexCandidates[s].push_back(e);
        }

        /*
         * Build trees for local candidate cycles
         */
        const bool sorted_cycles = false;
        SPTrees<Graph, WeightMap, true> sp_trees(g, weight_map, perVertexCandidates, sorted_cycles);

        /*
         * Main loop
         */
        WeightType mcb_weight = WeightType();
        for (std::size_t k = 0; k < csd; k++) {
            if (k % 250 == 0) {
                std::cout << "Rank " << world.rank() << " at cycle " << k << std::endl;
            }

            /*
             * Compute shortest odd cycle
             */
            std::set<Edge> best_cycle;
            WeightType best_cycle_weight;
            std::set<Edge> signed_edges;
            convert_edges(support[k], std::inserter(signed_edges, signed_edges.end()), forest_index);

            // TODO

            boost::tie(best_cycle, best_cycle_weight) = sp_trees.compute_shortest_odd_cycle(signed_edges);

            /*
             * Update support vectors
             */
            std::set<std::size_t> cyclek;
            convert_edges(best_cycle, std::inserter(cyclek, cyclek.end()), forest_index);
            for (std::size_t l = k + 1; l < csd; l++) {
                if (support[l] * cyclek == 1) {
                    support[l] += support[k];
                }
            }

            /*
             * Output new cycle
             */
            std::list<Edge> cyclek_edgelist;
            std::copy(best_cycle.begin(), best_cycle.end(), std::back_inserter(cyclek_edgelist));
            *out++ = cyclek_edgelist;
            mcb_weight += best_cycle_weight;
        }

        return mcb_weight;
    }

} // namespace mcb

#endif
