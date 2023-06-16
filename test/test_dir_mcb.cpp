//    Copyright (C) Dimitrios Michail 2019 - 2023.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <iostream>

#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>

#include <parmcb/parmcb.hpp>
#include <parmcb/util.hpp>

using namespace boost;
typedef adjacency_list<vecS, vecS, bidirectionalS, no_property, property<edge_weight_t, double> > Graph;
typedef graph_traits<Graph>::edge_descriptor Edge;

void create_graph(Graph& graph) {
    property_map<Graph, edge_weight_t>::type weight = get(edge_weight, graph);

    auto v0 = add_vertex(graph);
    auto v1 = add_vertex(graph);
    auto v2 = add_vertex(graph);
    auto v3 = add_vertex(graph);
    auto v4 = add_vertex(graph);
    auto v5 = add_vertex(graph);
    auto v6 = add_vertex(graph);
    auto v7 = add_vertex(graph);
    auto v8 = add_vertex(graph);
    auto v9 = add_vertex(graph);
    auto v10 = add_vertex(graph);

    auto e01 = add_edge(v0, v1, graph).first;
    weight[e01] = 1.0;
    auto e12 = add_edge(v1, v2, graph).first;
    weight[e12] = 1.0;
    auto e23 = add_edge(v2, v3, graph).first;
    weight[e23] = 100.0;
    auto e34 = add_edge(v3, v4, graph).first;
    weight[e34] = 1.0;
    auto e45 = add_edge(v4, v5, graph).first;
    weight[e45] = 1.0;
    auto e05 = add_edge(v0, v5, graph).first;
    weight[e05] = 1.0;
    auto e56 = add_edge(v5, v6, graph).first;
    weight[e56] = 1.0;
    auto e27 = add_edge(v2, v7, graph).first;
    weight[e27] = 2.0;
    auto e79 = add_edge(v7, v9, graph).first;
    weight[e79] = 1.0;
    auto e38 = add_edge(v3, v8, graph).first;
    weight[e38] = 5.0;
    auto e108 = add_edge(v10, v8, graph).first;
    weight[e108] = 1.0;
    auto e109 = add_edge(v10, v9, graph).first;
    weight[e109] = 1.0;
}

TEST_CASE("sequential dir sva signed"){
    Graph graph;
    create_graph(graph);
    property_map<Graph, edge_weight_t>::type weight = get(edge_weight, graph);

    CHECK(num_vertices(graph) == 11);
    CHECK(num_edges(graph) == 12);

    //
    //   0 -- 1 -- 2 -- 7 -- 9
    //   |         |         |
    //   5 -- 4 -- 3 -- 8 -- 10
    //   |
    //   6
    //

    std::list<std::list<Edge>> cycles;
    boost::multiprecision::cpp_int p = 41;
    double mcb_weight = parmcb::mcb_dir_sva_signed(graph, weight, p, std::back_inserter(cycles));

//
//    for (auto it = cycles.begin(); it != cycles.end(); it++) {
//        auto cycle = *it;
//        CHECK(parmcb::is_cycle(graph, cycle));
//    }
//
//    CHECK(cycles.size() == 3);
//    CHECK(mcb_weight == 124.0);

}

