#ifndef PARMCB_SVA_SIGNED_HPP_
#define PARMCB_SVA_SIGNED_HPP_

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

#include <parmcb/detail/signed_dijkstra.hpp>

#include <parmcb/forestindex.hpp>
#include <parmcb/spvecgf2.hpp>
#include <parmcb/util.hpp>

namespace parmcb {

    template<class Graph, class WeightMap, class CycleOutputIterator>
    typename boost::property_traits<WeightMap>::value_type mcb_sva_signed(const Graph &g, WeightMap weight_map,
            CycleOutputIterator out) {

        typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
        typedef typename boost::graph_traits<Graph>::vertex_iterator VertexIt;
        typedef typename boost::property_traits<WeightMap>::value_type WeightType;

        /*
         * Index the graph
         */
        ForestIndex<Graph> forest_index(g);
        auto csd = forest_index.cycle_space_dimension();

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

        /*
         * Main loop
         */
        WeightType mcb_weight = WeightType();
        for (std::size_t k = 0; k < csd; k++) {
            /*
             * Choose the sparsest support heuristic
             */
            auto min_support = k;
            for (auto r = k + 1; r < csd; ++r) {
                if (support[r].size() < support[min_support].size())
                    min_support = r;
                if (support[min_support].size() < 5) {
                    break;
                }
            }
            if (min_support != k) {  // swap
                std::swap(support[k], support[min_support]);
            }

            /*
             * Compute shortest odd cycle
             */
            cycle_timer.resume();
            std::less<WeightType> compare = std::less<WeightType>();
            std::tuple<std::set<Edge>, WeightType, bool> best = std::make_tuple(std::set<Edge>(),
                    (std::numeric_limits<WeightType>::max)(), false);
            std::set<Edge> signed_edges;
            convert_edges(support[k], std::inserter(signed_edges, signed_edges.end()), forest_index);

            if (signed_edges.size() >= boost::num_vertices(g)) {
                VertexIt vi, viend;
                for (boost::tie(vi, viend) = boost::vertices(g); vi != viend; ++vi) {
                    auto v = *vi;
                    const bool use_hidden_edges = false;
                    auto res = bidirectional_signed_dijkstra(g, weight_map, signed_edges, std::set<Edge> { },
                            use_hidden_edges, v, true, v, false, std::get<2>(best), std::get<1>(best));
                    if (std::get<2>(res) && (!std::get<2>(best) || compare(std::get<1>(res), std::get<1>(best)))) {
                        best = res;
                        assert(std::get<2>(best));
                    }
                }
            } else {
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
                            true, se_u, true, std::get<2>(best), std::get<1>(best));
                    hidden_edges.erase(hidden_edges.begin());
                    if (std::get<2>(res) && std::get<0>(res).find(se) == std::get<0>(res).end()) {
                        std::get<1>(res) += boost::get(weight_map, se);
                        if (!std::get<2>(best) || compare(std::get<1>(res), std::get<1>(best))) {
                            std::get<0>(res).insert(se);
                            best = res;
                            assert(std::get<2>(best));
                        }
                    }
                }
            }
            assert(std::get<2>(best));
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

        std::cout << "cycle   timer" << cycle_timer.format();
        std::cout << "support timer" << support_timer.format();

        return mcb_weight;
    }

} // parmcb

#endif
