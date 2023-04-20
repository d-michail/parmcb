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
#include <parmcb/parmcb_approx_sva_signed.hpp>

using namespace boost;

typedef adjacency_list<vecS, vecS, undirectedS, no_property, property<edge_weight_t, double> > Graph;
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
    auto v11 = add_vertex(graph);
    auto v12 = add_vertex(graph);
    auto v13 = add_vertex(graph);
    auto v14 = add_vertex(graph);
    auto v15 = add_vertex(graph);
    auto v16 = add_vertex(graph);
    (void) v16;

    auto e01 = add_edge(v0, v1, graph).first;
    weight[e01] = 1.0;
    auto e12 = add_edge(v1, v2, graph).first;
    weight[e12] = 1.0;
    auto e1_4 = add_edge(v1, v4, graph).first;
    weight[e1_4] = 1.0;
    auto e23 = add_edge(v2, v3, graph).first;
    weight[e23] = 100.0;
    auto e34 = add_edge(v3, v4, graph).first;
    weight[e34] = 1.0;
    auto e45 = add_edge(v4, v5, graph).first;
    weight[e45] = 1.0;
    auto e0_5 = add_edge(v0, v5, graph).first;
    weight[e0_5] = 50.0;
    auto e56 = add_edge(v5, v6, graph).first;
    weight[e56] = 1.0;
    auto e27 = add_edge(v2, v7, graph).first;
    weight[e27] = 2.0;
    auto e7_8 = add_edge(v7, v8, graph).first;
    weight[e7_8] = 1.0;
    auto e79 = add_edge(v7, v9, graph).first;
    weight[e79] = 1.0;
    auto e38 = add_edge(v3, v8, graph).first;
    weight[e38] = 5.0;
    auto e3_12 = add_edge(v3, v12, graph).first;
    weight[e3_12] = 1.0;
    auto e108 = add_edge(v10, v8, graph).first;
    weight[e108] = 1.0;
    auto e10_9 = add_edge(v10, v9, graph).first;
    weight[e10_9] = 50.0;
    auto e1112 = add_edge(v11, v12, graph).first;
    weight[e1112] = 3.0;
    auto e12_13 = add_edge(v12, v13, graph).first;
    weight[e12_13] = 1.0;
    auto e13_14 = add_edge(v13, v14, graph).first;
    weight[e13_14] = 1.0;
    auto e14_15 = add_edge(v14, v15, graph).first;
    weight[e14_15] = 1.0;
    auto e12_15 = add_edge(v12, v15, graph).first;
    weight[e12_15] = 1.0;
}


TEST_CASE("approx_mcb_sva_signed"){
    Graph graph;
    create_graph(graph);
    property_map<Graph, edge_weight_t>::type weight = get(edge_weight, graph);

    CHECK(num_vertices(graph) == 17);
    CHECK(num_edges(graph) == 20);

    //
    //   0 -- 1 --  2 -- 7 -- 9
    //   |    |     |    |    |
    //   5 -- 4 --  3 -- 8 -- 10
    //   |          |
    //   6    11 -- 12 -- 13
    //              |      |
    //  16          15 --  14
    //

    std::list<std::list<Edge>> cycles;
    double mcb_weight = parmcb::approx_mcb_sva_signed(graph, weight, 3, std::back_inserter(cycles));

    for (auto it = cycles.begin(); it != cycles.end(); it++) {
        auto cycle = *it;
        CHECK(parmcb::is_cycle(graph, cycle));
    }

    // DEBUG
    std::for_each(cycles.begin(), cycles.end(),
            [](const std::list<Edge> &cycle) {
                std::for_each(cycle.begin(), cycle.end(),
                        [](const Edge &e) {
                            std::cout << e;
                        });
                std::cout << std::endl;
            });


    CHECK(cycles.size() == 5);
    CHECK(mcb_weight == 224.0);
}
