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

#define GRID_ROWS 10
#define GRID_COLS 10
#define NUM_VERTICES (GRID_ROWS*GRID_COLS)
#define NUM_EDGES (2*GRID_ROWS*GRID_COLS-GRID_ROWS-GRID_COLS)
#define CYCLE_SPACE_DIM (NUM_EDGES-NUM_VERTICES+1)
#define MCB_WEIGHT ((GRID_ROWS-1)*(GRID_COLS-1)*4)

void create_graph(Graph& graph) {
    int rows = GRID_ROWS;
    int cols = GRID_COLS;

    property_map<Graph, edge_weight_t>::type weight = get(edge_weight, graph);

    // Add vertices
    std::vector<Graph::vertex_descriptor> vertices;
    for (int i = 0; i < rows * cols; ++i) {
        vertices.push_back(add_vertex(graph));
    }

    // Add edges horizontally
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols - 1; ++j) {
            auto e = add_edge(vertices[i * cols + j], vertices[i * cols + j + 1], graph).first;
            weight[e] = 1.0;
        }
    }

    // Add edges vertically
    for (int i = 0; i < rows - 1; ++i) {
        for (int j = 0; j < cols; ++j) {
            auto e = add_edge(vertices[i * cols + j], vertices[(i + 1) * cols + j], graph).first;
            weight[e] = 1.0;
        }
    }
}

TEST_CASE("sequential sva signed"){
    Graph graph;
    create_graph(graph);
    property_map<Graph, edge_weight_t>::type weight = get(edge_weight, graph);

    CHECK(num_vertices(graph) == NUM_VERTICES);
    CHECK(num_edges(graph) == NUM_EDGES);

    std::list<std::list<Edge>> cycles;
    double mcb_weight = parmcb::mcb_sva_signed(graph, weight, std::back_inserter(cycles));

    for (auto it = cycles.begin(); it != cycles.end(); it++) {
        auto cycle = *it;
        CHECK(parmcb::is_cycle(graph, cycle));
    }

    CHECK(cycles.size() == CYCLE_SPACE_DIM);
    CHECK(mcb_weight == MCB_WEIGHT);

}

#ifdef PARMCB_HAVE_TBB
TEST_CASE("parallel sva signed"){

    Graph graph;
    create_graph(graph);
    property_map<Graph, edge_weight_t>::type weight = get(edge_weight, graph);

    CHECK(num_vertices(graph) == NUM_VERTICES);
    CHECK(num_edges(graph) == NUM_EDGES);

     std::list<std::list<Edge>> cycles;
    double mcb_weight = parmcb::mcb_sva_signed_tbb(graph, weight, std::back_inserter(cycles));

    for (auto it = cycles.begin(); it != cycles.end(); it++) {
        auto cycle = *it;
        CHECK(parmcb::is_cycle(graph, cycle));
    }

    CHECK(cycles.size() == CYCLE_SPACE_DIM);
    CHECK(mcb_weight == MCB_WEIGHT);

}
#endif

TEST_CASE("sequential sva fvs trees"){

    Graph graph;
    create_graph(graph);
    property_map<Graph, edge_weight_t>::type weight = get(edge_weight, graph);

    CHECK(num_vertices(graph) == NUM_VERTICES);
    CHECK(num_edges(graph) == NUM_EDGES);

    std::list<std::list<Edge>> cycles;
    double mcb_weight = parmcb::mcb_sva_fvs_trees(graph, weight, std::back_inserter(cycles));

    for (auto it = cycles.begin(); it != cycles.end(); it++) {
        auto cycle = *it;
        CHECK(parmcb::is_cycle(graph, cycle));
    }

    CHECK(cycles.size() == CYCLE_SPACE_DIM);
    CHECK(mcb_weight == MCB_WEIGHT);

}

#ifdef PARMCB_HAVE_TBB
TEST_CASE("parallel sva fvs trees"){

    Graph graph;
    create_graph(graph);
    property_map<Graph, edge_weight_t>::type weight = get(edge_weight, graph);

    CHECK(num_vertices(graph) == NUM_VERTICES);
    CHECK(num_edges(graph) == NUM_EDGES);

    std::list<std::list<Edge>> cycles;
    double mcb_weight = parmcb::mcb_sva_fvs_trees_tbb(graph, weight, std::back_inserter(cycles));

    for (auto it = cycles.begin(); it != cycles.end(); it++) {
        auto cycle = *it;
        CHECK(parmcb::is_cycle(graph, cycle));
    }

    CHECK(cycles.size() == CYCLE_SPACE_DIM);
    CHECK(mcb_weight == MCB_WEIGHT);

}
#endif

TEST_CASE("sequential sva iso trees"){

    Graph graph;
    create_graph(graph);
    property_map<Graph, edge_weight_t>::type weight = get(edge_weight, graph);

    CHECK(num_vertices(graph) == NUM_VERTICES);
    CHECK(num_edges(graph) == NUM_EDGES);

    std::list<std::list<Edge>> cycles;
    double mcb_weight = parmcb::mcb_sva_iso_trees(graph, weight, std::back_inserter(cycles));

    for (auto it = cycles.begin(); it != cycles.end(); it++) {
        auto cycle = *it;
        CHECK(parmcb::is_cycle(graph, cycle));
    }

    CHECK(cycles.size() == CYCLE_SPACE_DIM);
    CHECK(mcb_weight == MCB_WEIGHT);

}

#ifdef PARMCB_HAVE_TBB
TEST_CASE("parallel sva iso trees"){

    Graph graph;
    create_graph(graph);
    property_map<Graph, edge_weight_t>::type weight = get(edge_weight, graph);

    CHECK(num_vertices(graph) == NUM_VERTICES);
    CHECK(num_edges(graph) == NUM_EDGES);

    std::list<std::list<Edge>> cycles;
    double mcb_weight = parmcb::mcb_sva_iso_trees_tbb(graph, weight, std::back_inserter(cycles));

    for (auto it = cycles.begin(); it != cycles.end(); it++) {
        auto cycle = *it;
        CHECK(parmcb::is_cycle(graph, cycle));
    }

    CHECK(cycles.size() == CYCLE_SPACE_DIM);
    CHECK(mcb_weight == MCB_WEIGHT);

}
#endif

