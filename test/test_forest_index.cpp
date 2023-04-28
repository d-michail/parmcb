//    Copyright (C) Dimitrios Michail 2019 - 2021.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <iostream>
#include <set>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>

#include <parmcb/forestindex.hpp>

TEST_CASE("forest index")
{

    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property,
            boost::property<boost::edge_weight_t, double> > Graph;

    Graph graph;

    auto v0 = boost::add_vertex(graph);
    auto v3 = boost::add_vertex(graph);
    auto v5 = boost::add_vertex(graph);
    auto v7 = boost::add_vertex(graph);
    auto v8 = boost::add_vertex(graph);
    auto v9 = boost::add_vertex(graph);
    auto v10 = boost::add_vertex(graph);
    auto v11 = boost::add_vertex(graph);
    auto v12 = boost::add_vertex(graph);
    auto v13 = boost::add_vertex(graph);
    auto v14 = boost::add_vertex(graph);
    auto v15 = boost::add_vertex(graph);
    boost::add_vertex(graph);
    boost::add_vertex(graph);
    auto v18 = boost::add_vertex(graph);

    boost::add_edge(v0, v7, graph);
    boost::add_edge(v0, v5, graph);
    boost::add_edge(v3, v5, graph);
    boost::add_edge(v3, v8, graph);
    boost::add_edge(v7, v9, graph);
    boost::add_edge(v8, v9, graph);
    boost::add_edge(v8, v10, graph);
    boost::add_edge(v9, v18, graph);
    boost::add_edge(v10, v18, graph);
    boost::add_edge(v11, v12, graph);
    boost::add_edge(v12, v13, graph);
    boost::add_edge(v12, v15, graph);
    boost::add_edge(v13, v14, graph);
    boost::add_edge(v14, v15, graph);

    CHECK(boost::num_vertices(graph) == 15);
    CHECK(boost::num_edges(graph) == 14);

    //
    //   0 -- 7 -- 9 -- 18
    //   |         |     |
    //   5 -- 3 -- 8 -- 10
    //
    //       11 -- 12 -- 13
    //              |      |
    //  16   17     15 --  14
    //

    parmcb::ForestIndex<Graph> fi(graph);

    CHECK(fi.weak_connected_components() == 4);
    CHECK(fi.cycle_space_dimension() == 3);


    int not_on_forest_count = 0;
    int on_forest_count = 0;
    std::set<std::size_t> used;
    for (const auto &e : boost::make_iterator_range(boost::edges(graph))) {
        auto index = fi(e);
        CHECK (used.count(index) == 0);
        used.insert(index);
        CHECK(fi(index) == e);
        CHECK(fi(e) == index);
        if (fi.is_on_forest(e)) {
            on_forest_count++;
        } else {
            not_on_forest_count++;
        }
    }
    CHECK(used.size() == 14);

    CHECK(not_on_forest_count == 3);
    CHECK(on_forest_count == 11);
}
