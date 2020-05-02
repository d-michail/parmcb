//    Copyright (C) Dimitrios Michail 2019 - 2020.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>

#include <parmcb/detail/fvs.hpp>

TEST_CASE("Feedback Vertex Set")
{

    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property,
            boost::property<boost::edge_weight_t, double> > Graph;

    Graph graph;

    auto v0 = boost::add_vertex(graph);
    auto v1 = boost::add_vertex(graph);
    auto v2 = boost::add_vertex(graph);
    auto v3 = boost::add_vertex(graph);
    auto v4 = boost::add_vertex(graph);
    auto v5 = boost::add_vertex(graph);
    auto v6 = boost::add_vertex(graph);
    auto v7 = boost::add_vertex(graph);
    auto v8 = boost::add_vertex(graph);
    auto v9 = boost::add_vertex(graph);
    auto v10 = boost::add_vertex(graph);
    auto v11 = boost::add_vertex(graph);
    auto v12 = boost::add_vertex(graph);
    auto v13 = boost::add_vertex(graph);
    auto v14 = boost::add_vertex(graph);
    boost::add_vertex(graph);
    boost::add_vertex(graph);

    boost::add_edge(v0, v3, graph);
    boost::add_edge(v0, v2, graph);
    boost::add_edge(v1, v2, graph);
    boost::add_edge(v1, v4, graph);
    boost::add_edge(v3, v5, graph);
    boost::add_edge(v4, v5, graph);
    boost::add_edge(v4, v6, graph);
    boost::add_edge(v5, v12, graph);
    boost::add_edge(v6, v12, graph);
    boost::add_edge(v7, v8, graph);
    boost::add_edge(v8, v9, graph);
    boost::add_edge(v8, v11, graph);
    boost::add_edge(v9, v10, graph);
    boost::add_edge(v10, v11, graph);
    boost::add_edge(v9, v13, graph);
    boost::add_edge(v14, v13, graph);
    boost::add_edge(v4, v8, graph);

    CHECK(boost::num_vertices(graph) == 17);
    CHECK(boost::num_edges(graph) == 17);

    //
    //   0 -- 3 -- 5 -- 12
    //   |         |    |
    //   2 -- 1 -- 4 -- 6
    //             |
    //        7 -- 8  -- 9 -- 13 -- 14
    //             |     |
    //  15   16    11 -- 10
    //

    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    std::set<Vertex> result;
    parmcb::greedy_fvs(graph, std::inserter(result, result.end()));

    CHECK(result.size() == 2);
    CHECK(result.count(4)>0);
    CHECK(result.count(9)>0);

}
