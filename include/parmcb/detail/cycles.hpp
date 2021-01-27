#ifndef PARMCB_DETAIL_CYCLES_HPP_
#define PARMCB_DETAIL_CYCLES_HPP_

//    Copyright (C) Dimitrios Michail 2019 - 2021.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <map>
#include <parmcb/util.hpp>
#include <parmcb/sptrees.hpp>
#include <parmcb/detail/fvs.hpp>

#include <boost/graph/connected_components.hpp>

namespace boost {

    struct bad_t {
        typedef boost::vertex_property_tag kind;
    };

    struct tree_t {
        typedef boost::vertex_property_tag kind;
    };

    struct edge_t {
        typedef boost::vertex_property_tag kind;
    };

}

namespace parmcb {

    namespace detail {

        template<class Graph, class WeightMap>
        struct HortonCyclesBuilder {

            void operator()(const Graph &g, const WeightMap &weight_map,
                    std::vector<parmcb::SPTree<Graph, WeightMap>> &trees,
                    std::vector<CandidateCycle<Graph, WeightMap>> &cycles) {
                typedef typename boost::graph_traits<Graph>::vertex_iterator VertexIt;
                VertexIt vi, viend;
                for (boost::tie(vi, viend) = boost::vertices(g); vi != viend; ++vi) {
                    auto v = *vi;
                    trees.emplace_back(trees.size(), g, boost::get(boost::vertex_index, g), weight_map, v);
                }
                for (auto &tree : trees) {
                    std::vector<CandidateCycle<Graph, WeightMap>> tree_cycles = tree.create_candidate_cycles();
                    cycles.insert(cycles.end(), tree_cycles.begin(), tree_cycles.end());
                }
            }

        };

        template<class Graph, class WeightMap>
        struct FVSCyclesBuilder {

            void operator()(const Graph &g, const WeightMap &weight_map,
                    std::vector<parmcb::SPTree<Graph, WeightMap>> &trees,
                    std::vector<CandidateCycle<Graph, WeightMap>> &cycles) {
                typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;

                std::vector<Vertex> feedback_vertex_set;
                parmcb::greedy_fvs(g, std::back_inserter(feedback_vertex_set));
                for (auto v : feedback_vertex_set) {
                    trees.emplace_back(trees.size(), g, boost::get(boost::vertex_index, g), weight_map, v);
                }
                for (auto &tree : trees) {
                    std::vector<CandidateCycle<Graph, WeightMap>> tree_cycles = tree.create_candidate_cycles();
                    cycles.insert(cycles.end(), tree_cycles.begin(), tree_cycles.end());
                }
            }

        };

        template<class Graph, class WeightMap>
        struct ISOCyclesBuilder {

            void operator()(const Graph &g, const WeightMap &weight_map,
                    std::vector<parmcb::SPTree<Graph, WeightMap>> &trees,
                    std::vector<CandidateCycle<Graph, WeightMap>> &cycles) {
                typedef typename boost::property_map<Graph, boost::vertex_index_t>::type VertexIndexMapType;
                typedef typename boost::graph_traits<Graph>::vertex_iterator VertexIt;
                typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
                typedef typename boost::property_traits<WeightMap>::value_type WeightType;

                const VertexIndexMapType &index_map = boost::get(boost::vertex_index, g);

                /*
                 * Build all shortest path trees
                 */
                std::vector<std::size_t> trees_index_map(boost::num_vertices(g));
                VertexIt ui, uiend;
                std::size_t next_tree = 0;
                for (boost::tie(ui, uiend) = boost::vertices(g); ui != uiend; ++ui) {
                    auto u = *ui;
                    auto uindex = index_map[u];
                    trees.emplace_back(next_tree, g, boost::get(boost::vertex_index, g), weight_map, u);
                    trees_index_map[uindex] = next_tree;
                    next_tree++;
                }

                /*
                 * Build all of Horton's candidate cycles
                 */
                std::vector<CandidateCycle<Graph, WeightMap>> allcycles;
                for (const auto &tree : trees) {
                    std::vector<CandidateCycle<Graph, WeightMap>> tree_cycles = tree.create_candidate_cycles();
                    allcycles.insert(allcycles.end(), tree_cycles.begin(), tree_cycles.end());
                }

                /*
                 * Create graph with candidate cycles
                 */
                typedef boost::property<boost::tree_t, std::size_t> TreeVertexProperty;
                typedef boost::property<boost::edge_t, Edge, TreeVertexProperty> EdgeVertexProperty;
                typedef boost::property<boost::bad_t, bool, EdgeVertexProperty> BadVertexProperty;
                typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, BadVertexProperty> graph_t;
                typedef typename boost::graph_traits<graph_t>::vertex_descriptor vertex_descriptor;
                typedef typename boost::graph_traits<graph_t>::vertex_iterator vertex_iterator;

                graph_t cycles_g;
                typename boost::property_map<graph_t, boost::tree_t>::type tree_map = boost::get(boost::tree_t(),
                        cycles_g);
                typename boost::property_map<graph_t, boost::edge_t>::type edge_map = boost::get(boost::edge_t(),
                        cycles_g);
                typename boost::property_map<graph_t, boost::bad_t>::type bad_map = boost::get(boost::bad_t(),
                        cycles_g);

                std::map<std::pair<std::size_t, Edge>, vertex_descriptor> cycle_to_vertex;
                for (const auto &cc : allcycles) {
                    vertex_descriptor newv = boost::add_vertex(cycles_g);
                    cycle_to_vertex[std::make_pair(cc.tree(), cc.edge())] = newv;
                    boost::put(tree_map, newv, cc.tree());
                    boost::put(edge_map, newv, cc.edge());
                    boost::put(bad_map, newv, false);
                }

                vertex_iterator alli, alliend;
                for (boost::tie(alli, alliend) = boost::vertices(cycles_g); alli != alliend; ++alli) {
                    auto tree = boost::get(tree_map, *alli);
                    auto e = boost::get(edge_map, *alli);

                    parmcb::SPTree<Graph, WeightMap> &tree_x = trees[tree];
                    auto x = tree_x.source();
                    auto u = boost::source(e, g);
                    auto v = boost::target(e, g);

                    auto first_x_u = tree_x.first(u);
                    auto first_x_v = tree_x.first(v);

                    if (first_x_u == first_x_v) {
                        continue;
                    }

                    if (x == u) {
                        boost::add_edge(*alli, cycle_to_vertex[std::make_pair(trees_index_map[index_map[v]], e)],
                                cycles_g);
                    } else {
                        auto xprime = tree_x.first(u);
                        parmcb::SPTree<Graph, WeightMap> &tree_xprime = trees[trees_index_map[index_map[xprime]]];

                        auto first_xprime_v = tree_xprime.first(v);

                        if (x == first_xprime_v) {
                            boost::add_edge(*alli,
                                    cycle_to_vertex[std::make_pair(trees_index_map[index_map[xprime]], e)], cycles_g);
                        } else {
                            parmcb::SPTree<Graph, WeightMap> &tree_v = trees[trees_index_map[index_map[v]]];
                            auto first_v_xprime = tree_v.first(xprime);
                            if (u == first_v_xprime) {
                                boost::add_edge(*alli,
                                        cycle_to_vertex[std::make_pair(trees_index_map[index_map[v]],
                                                tree_x.node(xprime)->pred())], cycles_g);
                            } else {
                                boost::put(bad_map, *alli, true);
                            }
                        }
                    }
                }

                /*
                 * Find components
                 */
                std::vector<std::size_t> components(boost::num_vertices(cycles_g));
                boost::function_property_map<parmcb::detail::VertexIndexFunctor<graph_t, std::size_t>,
                        vertex_descriptor, std::size_t&> components_map(
                        parmcb::detail::VertexIndexFunctor<graph_t, std::size_t>(components,
                                boost::get(boost::vertex_index, cycles_g)));
                std::size_t num_components = boost::connected_components(cycles_g, components_map);

                std::vector<bool> is_bad_component(num_components, false);
                for (boost::tie(alli, alliend) = boost::vertices(cycles_g); alli != alliend; ++alli) {
                    auto v = *alli;
                    std::size_t component = boost::get(components_map, v);
                    is_bad_component[component] = is_bad_component[component] || boost::get(bad_map, v);
                }

                std::vector<bool> is_in_output(num_components, false);
                for (boost::tie(alli, alliend) = boost::vertices(cycles_g); alli != alliend; ++alli) {
                    auto v = *alli;
                    std::size_t component = boost::get(components_map, v);
                    if (!is_bad_component[component] && !is_in_output[component]) {
                        auto tree = boost::get(tree_map, v);
                        auto e = boost::get(edge_map, v);
                        parmcb::SPTree<Graph, WeightMap> &tree_v = trees[tree];

                        std::shared_ptr<SPNode<Graph, WeightMap>> spnode_v = tree_v.node(boost::source(e, g));
                        if (spnode_v == nullptr) {
                            continue;
                        }
                        std::shared_ptr<SPNode<Graph, WeightMap>> spnode_u = tree_v.node(boost::target(e, g));
                        if (spnode_u == nullptr) {
                            continue;
                        }

                        WeightType cycle_weight = boost::get(weight_map, e) + spnode_v->weight() + spnode_u->weight();
                        cycles.emplace_back(tree, e, cycle_weight);
                        is_in_output[component] = true;
                    }
                }
            }
        };

    } // detail

} // parmcb

#endif
