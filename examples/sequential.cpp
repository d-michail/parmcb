//    Copyright (C) Dimitrios Michail 2019 - 2021.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <fstream>
#include <map>
#include <list>

#include <boost/config.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>

#include <parmcb/parmcb.hpp>

using namespace boost;

int main(int argc, char *argv[]) {

    typedef adjacency_list<vecS, vecS, undirectedS, no_property, property<edge_weight_t, double> > Graph;
    typedef graph_traits<Graph>::edge_descriptor Edge;


    //
    //   0 -- 1 -- 2 -- 6 -- 9
    //   |         |         |
    //   5 -- 4 -- 3 -- 7 -- 8
    //

    Graph graph;

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

    auto e01 = add_edge(v0, v1, graph).first;
    auto e12 = add_edge(v1, v2, graph).first;
    auto e23 = add_edge(v2, v3, graph).first;
    auto e34 = add_edge(v3, v4, graph).first;
    auto e45 = add_edge(v4, v5, graph).first;
    auto e05 = add_edge(v0, v5, graph).first;
    auto e26 = add_edge(v2, v6, graph).first;
    auto e37 = add_edge(v3, v7, graph).first;
    auto e69 = add_edge(v6, v9, graph).first;
    auto e78 = add_edge(v7, v8, graph).first;
    auto e89 = add_edge(v8, v9, graph).first;

    property_map<Graph, edge_weight_t>::type weight = get(edge_weight, graph);

    weight[e01] = 1.0;
    weight[e12] = 1.0;
    weight[e23] = 100.0;
    weight[e34] = 1.0;
    weight[e45] = 1.0;
    weight[e05] = 1.0;
    weight[e26] = 1.0;
    weight[e37] = 2.0;
    weight[e69] = 1.0;
    weight[e78] = 5.0;
    weight[e89] = 1.0;

    std::list<std::list<Edge>> cycles;
    double mcb_weight = parmcb::mcb_sva_signed(graph, weight, std::back_inserter(cycles));

    std::cout << "MCB weight = " << mcb_weight << std::endl;
    std::cout << "MCB cycles" << std::endl;
    for (auto it = cycles.begin(); it != cycles.end(); it++) {
        auto cycle = *it;

        for (auto eit = cycle.begin(); eit != cycle.end(); eit++) {
            std::cout << *eit << " ";
        }
        std::cout << std::endl;

    }

    return EXIT_SUCCESS;
}

