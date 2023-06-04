//    Copyright (C) Dimitrios Michail 2019 - 2023.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <iostream>
#include <set>

#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>


#include <parmcb/detail/util.hpp>
#include <parmcb/detail/dijkstra.hpp>

TEST_CASE("as undirected dijkstra for directed graphs")
{
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, boost::no_property,
            boost::property<boost::edge_weight_t, double> > Graph;
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename boost::graph_traits<Graph>::vertex_iterator VertexIt;
    typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;

    Graph graph;

    const auto &index_map = boost::get(boost::vertex_index, graph);
    boost::property_map<Graph, boost::edge_weight_t>::type weight = get(boost::edge_weight, graph);

    auto v0 = boost::add_vertex(graph);
    auto v1 = boost::add_vertex(graph);
    auto v2 = boost::add_vertex(graph);
    auto v3 = boost::add_vertex(graph);
    auto v4 = boost::add_vertex(graph);
    auto v5 = boost::add_vertex(graph);
    auto v6 = boost::add_vertex(graph);
    auto v7 = boost::add_vertex(graph);

    auto e01 = boost::add_edge(v0, v1, graph).first;
    weight[e01] = 1.0;
    auto e12 = boost::add_edge(v1, v2, graph).first;
    weight[e12] = 1.0;
    auto e23 = boost::add_edge(v2, v3, graph).first;
    weight[e23] = 1.0;
    auto e34 = boost::add_edge(v3, v4, graph).first;
    weight[e34] = 1.0;
    auto e45 = boost::add_edge(v4, v5, graph).first;
    weight[e45] = 1.0;
    auto e56 = boost::add_edge(v5, v6, graph).first;
    weight[e56] = 1.0;
    auto e67 = boost::add_edge(v6, v7, graph).first;
    weight[e67] = 1.0;
    auto e70 = boost::add_edge(v7, v0, graph).first;
    weight[e70] = 1.0;

    std::vector<double> dist(boost::num_vertices(graph), (std::numeric_limits<double>::max)());
    boost::function_property_map<
		parmcb::detail::VertexIndexFunctor<Graph, double>,
        Vertex,
		double&> dist_map(parmcb::detail::VertexIndexFunctor<Graph, double>(dist, index_map));

	std::vector<std::tuple<bool, Edge>> pred(boost::num_vertices(graph), std::make_tuple(false, Edge()));
    boost::function_property_map<
            parmcb::detail::VertexIndexFunctor<Graph,
                    std::tuple<bool, Edge>>, Vertex,
            std::tuple<bool, Edge>&> pred_map(
            parmcb::detail::VertexIndexFunctor<Graph,
                    std::tuple<bool, Edge> >(pred, index_map));

    parmcb::as_undirected_dijkstra(graph, weight, v0, dist_map, pred_map);

	VertexIt ui, uiend;
	for (boost::tie(ui, uiend) = boost::vertices(graph); ui != uiend; ++ui) {
		auto u = *ui;
		std::cout << "dist[" << u << "]=" << dist[u] << std::endl;
	}

	CHECK(dist[0] == 0.0);
	CHECK(dist[1] == 1.0);
	CHECK(dist[7] == 1.0);
	CHECK(dist[2] == 2.0);
	CHECK(dist[6] == 2.0);
	CHECK(dist[3] == 3.0);
	CHECK(dist[5] == 3.0);
	CHECK(dist[4] == 4.0);

    CHECK(boost::num_vertices(graph) == 8);
    CHECK(boost::num_edges(graph) == 8);

}
