#ifndef PARMCB_SVA_TREES_HPP_
#define PARMCB_SVA_TREES_HPP_

#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/tuple/detail/tuple_basic.hpp>
#include <boost/timer/timer.hpp>

#include <cstddef>
#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <set>
#include <vector>

#include <parmcb/forestindex.hpp>
#include <parmcb/spvecgf2.hpp>
#include <parmcb/util.hpp>
#include <parmcb/sptrees.hpp>
#include <parmcb/detail/cycles.hpp>

namespace parmcb {

    template<class Graph, class WeightMap, class CycleOutputIterator, class CyclesBuilder, bool ParallelUsingTBB>
    typename boost::property_traits<WeightMap>::value_type _mcb_sva_trees(const Graph &g, WeightMap weight_map,
            CycleOutputIterator out) {
        typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
        typedef typename boost::property_traits<WeightMap>::value_type WeightType;

        /*
         * Index the graph
         */
        ForestIndex<Graph> forest_index(g);
        auto csd = forest_index.cycle_space_dimension();
        std::cout << "Cycle space dimension: " << csd << std::endl;

        /*
         * Initialize support vectors
         */
        std::vector<SpVecGF2<std::size_t>> support;
        for (std::size_t k = 0; k < csd; k++) {
            support.emplace_back(k);
        }

        boost::timer::cpu_timer cycle_timer;
        cycle_timer.stop();
        boost::timer::cpu_timer support_timer;
        support_timer.stop();
        boost::timer::cpu_timer trees_timer;
        trees_timer.stop();

        /*
         * Initialize all shortest path trees
         */
        trees_timer.resume();
        std::vector<parmcb::SPTree<Graph, WeightMap>> trees;
        std::vector<parmcb::CandidateCycle<Graph, WeightMap>> cycles;
        CyclesBuilder cycles_builder;
        cycles_builder(g, weight_map, trees, cycles);
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
        trees_timer.stop();

        /*
         * Main loop
         */
        WeightType mcb_weight = WeightType();
        for (std::size_t k = 0; k < csd; k++) {
            if (k % 250 == 0) {
                std::cout << k << std::endl;
            }

            /*
             * Compute shortest odd cycle
             */
            std::set<Edge> signed_edges;
            convert_edges(support[k], std::inserter(signed_edges, signed_edges.end()), forest_index);
            cycle_timer.resume();
            std::tuple<std::set<Edge>, WeightType, bool> best = cycle_lookup(signed_edges);
            cycle_timer.stop();

            /*
             * Update support vectors
             */
            support_timer.resume();
            std::set<std::size_t> cyclek;
            convert_edges(std::get<0>(best), std::inserter(cyclek, cyclek.end()), forest_index);
            for (std::size_t l = k + 1; l < csd; l++) {
                if (support[l] * cyclek == 1) {
                    support[l] += support[k];
                }
            }
            support_timer.stop();

            /*
             * Output new cycle
             */
            std::list<Edge> cyclek_edgelist;
            std::copy(std::get<0>(best).begin(), std::get<0>(best).end(), std::back_inserter(cyclek_edgelist));
            *out++ = cyclek_edgelist;
            mcb_weight += std::get<1>(best);
        }

        std::cout << "trees   timer" << trees_timer.format();
        std::cout << "cycle   timer" << cycle_timer.format();
        std::cout << "support timer" << support_timer.format();

        return mcb_weight;
    }

    template<class Graph, class WeightMap, class CycleOutputIterator>
    typename boost::property_traits<WeightMap>::value_type mcb_sva_fvs_trees(const Graph &g, WeightMap weight_map,
            CycleOutputIterator out) {
        return _mcb_sva_trees<Graph, WeightMap, CycleOutputIterator, parmcb::detail::FVSCyclesBuilder<Graph, WeightMap>,
                false>(g, weight_map, out);
    }

    template<class Graph, class WeightMap, class CycleOutputIterator>
    typename boost::property_traits<WeightMap>::value_type mcb_sva_fvs_trees_tbb(const Graph &g, WeightMap weight_map,
            CycleOutputIterator out) {
        return _mcb_sva_trees<Graph, WeightMap, CycleOutputIterator, parmcb::detail::FVSCyclesBuilder<Graph, WeightMap>,
                true>(g, weight_map, out);
    }

    template<class Graph, class WeightMap, class CycleOutputIterator>
    typename boost::property_traits<WeightMap>::value_type mcb_sva_iso_trees(const Graph &g, WeightMap weight_map,
            CycleOutputIterator out) {
        return _mcb_sva_trees<Graph, WeightMap, CycleOutputIterator, parmcb::detail::ISOCyclesBuilder<Graph, WeightMap>,
                false>(g, weight_map, out);
    }

    template<class Graph, class WeightMap, class CycleOutputIterator>
    typename boost::property_traits<WeightMap>::value_type mcb_sva_iso_trees_tbb(const Graph &g, WeightMap weight_map,
            CycleOutputIterator out) {
        return _mcb_sva_trees<Graph, WeightMap, CycleOutputIterator, parmcb::detail::ISOCyclesBuilder<Graph, WeightMap>,
                true>(g, weight_map, out);
    }

} // namespace parmcb

#endif
