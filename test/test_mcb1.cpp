//    Copyright (C) Dimitrios Michail 2019 - 2024.
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
typedef adjacency_list<vecS, vecS, undirectedS, no_property, property<edge_weight_t, double> > Graph;
typedef graph_traits<Graph>::edge_descriptor Edge;

void create_graph(Graph& graph) {
    property_map<Graph, edge_weight_t>::type weight = get(edge_weight, graph);

    auto v1 = add_vertex(graph);
    auto v2 = add_vertex(graph);
    auto v3 = add_vertex(graph);
    auto v4 = add_vertex(graph);
    auto v5 = add_vertex(graph);
    auto v6 = add_vertex(graph);
    auto v7 = add_vertex(graph);
    auto v8 = add_vertex(graph);

    auto e12 = add_edge(v1, v2, graph).first;
    weight[e12] = 1.0;
    auto e15 = add_edge(v1, v5, graph).first;
    weight[e15] = 1.0;
    auto e17 = add_edge(v1, v7, graph).first;
    weight[e17] = 1.0;
    auto e23 = add_edge(v2, v3, graph).first;
    weight[e23] = 1.0;
    auto e24 = add_edge(v2, v4, graph).first;
    weight[e24] = 1.0;
    auto e38 = add_edge(v3, v8, graph).first;
    weight[e38] = 1.0;
    auto e36 = add_edge(v3, v6, graph).first;
    weight[e36] = 1.0;
    auto e57 = add_edge(v5, v7, graph).first;
    weight[e57] = 1.0;
    auto e58 = add_edge(v5, v8, graph).first;
    weight[e58] = 1.0;
    auto e84 = add_edge(v8, v4, graph).first;
    weight[e84] = 1.0;
    auto e46 = add_edge(v4, v6, graph).first;
    weight[e46] = 1.0;
    auto e76 = add_edge(v7, v6, graph).first;
    weight[e76] = 1.0;
}

TEST_CASE("sequential sva signed"){
    Graph graph;
    create_graph(graph);
    property_map<Graph, edge_weight_t>::type weight = get(edge_weight, graph);

    CHECK(num_vertices(graph) == 8);
    CHECK(num_edges(graph) == 12);

    std::list<std::list<Edge>> cycles;
    double mcb_weight = parmcb::mcb_sva_signed(graph, weight, std::back_inserter(cycles));

    for (auto it = cycles.begin(); it != cycles.end(); it++) {
        auto cycle = *it;
        CHECK(parmcb::is_cycle(graph, cycle));
    }

    CHECK(cycles.size() == 5);
    CHECK(mcb_weight == 21.0);

}

#ifdef PARMCB_HAVE_TBB
TEST_CASE("parallel sva signed"){

    Graph graph;
    create_graph(graph);
    property_map<Graph, edge_weight_t>::type weight = get(edge_weight, graph);

    CHECK(num_vertices(graph) == 8);
    CHECK(num_edges(graph) == 12);

     std::list<std::list<Edge>> cycles;
    double mcb_weight = parmcb::mcb_sva_signed_tbb(graph, weight, std::back_inserter(cycles));

    for (auto it = cycles.begin(); it != cycles.end(); it++) {
        auto cycle = *it;
        CHECK(parmcb::is_cycle(graph, cycle));
    }

    CHECK(cycles.size() == 5);
    CHECK(mcb_weight == 21.0);

}
#endif

TEST_CASE("sequential sva fvs trees"){

    Graph graph;
    create_graph(graph);
    property_map<Graph, edge_weight_t>::type weight = get(edge_weight, graph);

    CHECK(num_vertices(graph) == 8);
    CHECK(num_edges(graph) == 12);

    std::list<std::list<Edge>> cycles;
    double mcb_weight = parmcb::mcb_sva_fvs_trees(graph, weight, std::back_inserter(cycles));

    for (auto it = cycles.begin(); it != cycles.end(); it++) {
        auto cycle = *it;
        CHECK(parmcb::is_cycle(graph, cycle));
    }

    CHECK(cycles.size() == 5);
    CHECK(mcb_weight == 21.0);

}

#ifdef PARMCB_HAVE_TBB
TEST_CASE("parallel sva fvs trees"){

    Graph graph;
    create_graph(graph);
    property_map<Graph, edge_weight_t>::type weight = get(edge_weight, graph);

    CHECK(num_vertices(graph) == 8);
    CHECK(num_edges(graph) == 12);

    std::list<std::list<Edge>> cycles;
    double mcb_weight = parmcb::mcb_sva_fvs_trees_tbb(graph, weight, std::back_inserter(cycles));

    for (auto it = cycles.begin(); it != cycles.end(); it++) {
        auto cycle = *it;
        CHECK(parmcb::is_cycle(graph, cycle));
    }

    CHECK(cycles.size() == 5);
    CHECK(mcb_weight == 21.0);

}
#endif

TEST_CASE("sequential sva iso trees"){

    Graph graph;
    create_graph(graph);
    property_map<Graph, edge_weight_t>::type weight = get(edge_weight, graph);

    CHECK(num_vertices(graph) == 8);
    CHECK(num_edges(graph) == 12);

    std::list<std::list<Edge>> cycles;
    double mcb_weight = parmcb::mcb_sva_iso_trees(graph, weight, std::back_inserter(cycles));

    for (auto it = cycles.begin(); it != cycles.end(); it++) {
        auto cycle = *it;
        CHECK(parmcb::is_cycle(graph, cycle));
    }

    CHECK(cycles.size() == 5);
    CHECK(mcb_weight == 21.0);

}

#ifdef PARMCB_HAVE_TBB
TEST_CASE("parallel sva iso trees"){

    Graph graph;
    create_graph(graph);
    property_map<Graph, edge_weight_t>::type weight = get(edge_weight, graph);

    CHECK(num_vertices(graph) == 8);
    CHECK(num_edges(graph) == 12);

    std::list<std::list<Edge>> cycles;
    double mcb_weight = parmcb::mcb_sva_iso_trees_tbb(graph, weight, std::back_inserter(cycles));

    for (auto it = cycles.begin(); it != cycles.end(); it++) {
        auto cycle = *it;
        CHECK(parmcb::is_cycle(graph, cycle));
    }

    CHECK(cycles.size() == 5);
    CHECK(mcb_weight == 21.0);

}
#endif

