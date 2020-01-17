#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>

#include <mcb/forestindex.hpp>

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

    auto e0_7 = boost::add_edge(v0, v7, graph).first;
    auto e0_5 = boost::add_edge(v0, v5, graph).first;
    auto e3_5 = boost::add_edge(v3, v5, graph).first;
    auto e3_8 = boost::add_edge(v3, v8, graph).first;
    auto e7_9 = boost::add_edge(v7, v9, graph).first;
    auto e8_9 = boost::add_edge(v8, v9, graph).first;
    auto e8_10 = boost::add_edge(v8, v10, graph).first;
    auto e9_18 = boost::add_edge(v9, v18, graph).first;
    auto e10_18 = boost::add_edge(v10, v18, graph).first;
    auto e11_12 = boost::add_edge(v11, v12, graph).first;
    auto e12_13 = boost::add_edge(v12, v13, graph).first;
    auto e12_15 = boost::add_edge(v12, v15, graph).first;
    auto e13_14 = boost::add_edge(v13, v14, graph).first;
    auto e14_15 = boost::add_edge(v14, v15, graph).first;

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

    mcb::ForestIndex<Graph> fi(graph);

    CHECK(fi.weak_connected_components() == 4);
    CHECK(fi.cycle_space_dimension() == 3);


    CHECK(fi(e3_5) == 0);
    CHECK(fi(e8_10) == 1);
    CHECK(fi(e13_14) == 2);
    CHECK(fi(e0_7) == 3);
    CHECK(fi(e0_5) == 4);
    CHECK(fi(e3_8) == 5);
    CHECK(fi(e7_9) == 6);
    CHECK(fi(e8_9) == 7);
    CHECK(fi(e9_18) == 8);
    CHECK(fi(e10_18) == 9);
    CHECK(fi(e11_12) == 10);
    CHECK(fi(e12_13) == 11);
    CHECK(fi(e12_15) == 12);
    CHECK(fi(e14_15) == 13);

    CHECK(fi(0) == e3_5);
    CHECK(fi(1) == e8_10);
    CHECK(fi(2) == e13_14);
    CHECK(fi(3) == e0_7);
    CHECK(fi(4) == e0_5);
    CHECK(fi(5) == e3_8);
    CHECK(fi(6) == e7_9);
    CHECK(fi(7) == e8_9);
    CHECK(fi(8) == e9_18);
    CHECK(fi(9) == e10_18);
    CHECK(fi(10) == e11_12);
    CHECK(fi(11) == e12_13);
    CHECK(fi(12) == e12_15);
    CHECK(fi(13) == e14_15);

    CHECK(!fi.is_on_forest(e3_5));
    CHECK(!fi.is_on_forest(e8_10));
    CHECK(!fi.is_on_forest(e13_14));
    CHECK(fi.is_on_forest(e0_7));
    CHECK(fi.is_on_forest(e0_5));
    CHECK(fi.is_on_forest(e3_8));
    CHECK(fi.is_on_forest(e7_9));
    CHECK(fi.is_on_forest(e8_9));
    CHECK(fi.is_on_forest(e9_18));
    CHECK(fi.is_on_forest(e10_18));
    CHECK(fi.is_on_forest(e11_12));
    CHECK(fi.is_on_forest(e12_13));
    CHECK(fi.is_on_forest(e12_15));
    CHECK(fi.is_on_forest(e14_15));

}
