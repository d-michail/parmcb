#ifndef PARMCB_UTIL_HPP_
#define PARMCB_UTIL_HPP_

//    Copyright (C) Dimitrios Michail 2019 - 2021.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include <set>
#include <map>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>

#include <parmcb/config.hpp>
#include <parmcb/forestindex.hpp>

#ifdef PARMCB_HAVE_TBB
#include <tbb/tbb.h>
#endif

#define BUFFER_SIZE 1024

namespace parmcb {

    template<class Graph>
    bool has_loops(const Graph &g) {
        auto eRange = boost::edges(g);
        for (auto eit = eRange.first; eit != eRange.second; ++eit) {
            auto e = *eit;
            auto v = boost::source(e, g);
            auto u = boost::target(e, g);
            if (v == u) {
                return true;
            }
        }
        return false;
    }

    template<class Graph>
    bool has_multiple_edges(const Graph &g) {
        typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;

        auto vRange = boost::vertices(g);
        for (auto vit = vRange.first; vit != vRange.second; ++vit) {
            auto v = *vit;
            std::set<Vertex> neighbors;

            auto eRange = boost::out_edges(v, g);
            for (auto eit = eRange.first; eit != eRange.second; ++eit) {
                auto e = *eit;
                auto u = boost::opposite(e, v, g);
                if (!neighbors.insert(u).second) {
                    return true;
                }
            }
        }
        return false;
    }

    template<class Graph, class WeightMap>
    bool has_non_positive_weights(const Graph &g, const WeightMap &weight_map) {
        auto eRange = boost::edges(g);
        for (auto eit = eRange.first; eit != eRange.second; ++eit) {
            auto e = *eit;
            if (boost::get(weight_map, e) <= 0.0) {
                return true;
            }
        }
        return false;
    }

    template<class Graph>
    bool is_cycle(const Graph &g, const std::list<typename boost::graph_traits<Graph>::edge_descriptor> &cycle) {
        if (cycle.size() == 0)
            return false;

        std::vector<std::size_t> degrees(boost::num_vertices(g));
        auto index_map = boost::get(boost::vertex_index, g);
        auto degrees_map = boost::make_iterator_property_map(degrees.begin(), index_map);

        for (auto eit = cycle.begin(); eit != cycle.end(); eit++) {
            auto e = *eit;
            degrees_map[boost::source(e, g)]++;
            degrees_map[boost::target(e, g)]++;
        }

        typename boost::graph_traits<Graph>::vertex_iterator vi, viend;
        for (boost::tie(vi, viend) = boost::vertices(g); vi != viend; ++vi) {
            if (degrees_map[*vi] % 2 == 1) {
                return false;
            }
        }

        return true;
    }

    template<class Graph, class InputIterator, class OutputIterator>
    void convert_edges(InputIterator inIt, OutputIterator outIt, const ForestIndex<Graph> &forest_index) {
        std::for_each(inIt.begin(), inIt.end(), [&outIt, &forest_index](auto i) {
            outIt = forest_index(i);
        });
    }

    template<class Graph>
    void read_dimacs_from_file(FILE *fp, Graph &graph) {
        typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
        typedef typename boost::graph_traits<Graph>::edge_descriptor edge_descriptor;

        char buffer[BUFFER_SIZE], problem[BUFFER_SIZE];

        std::size_t nnodes, nedges;
        std::map<std::size_t, vertex_descriptor> vertex_map;

        typename boost::property_map<Graph, boost::edge_weight_t>::type weight = get(boost::edge_weight, graph);

        while (fgets(buffer, sizeof(buffer), fp) != NULL) {
            buffer[strlen(buffer) - 1] = '\0'; // eat the newline
            if (buffer[0] == 'c' || buffer[0] == '#') {
                continue;
            } else if (buffer[0] == 'p') {
                sscanf(buffer, "p %s %lu %lu", problem, &nnodes, &nedges);
                for (std::size_t i = 1; i <= nnodes; i++) {
                    vertex_map[i] = boost::add_vertex(graph);
                }
            } else if (buffer[0] == 'a' || buffer[0] == 'e') {
                int rs, rt;
                double rw = 1;
                char fc;
                sscanf(buffer, "%c %d %d %lf", &fc, &rs, &rt, &rw);

                if (vertex_map.find(rs) == vertex_map.end()) {
                    throw std::system_error(EIO, std::generic_category(), "Vertex not found");
                }
                if (vertex_map.find(rt) == vertex_map.end()) {
                    throw std::system_error(EIO, std::generic_category(), "Vertex not found");
                }
                vertex_descriptor sDescriptor = boost::vertex(vertex_map[rs], graph);
                vertex_descriptor tDescriptor = boost::vertex(vertex_map[rt], graph);
                edge_descriptor eDescriptor = boost::add_edge(sDescriptor, tDescriptor, graph).first;
                weight[eDescriptor] = rw;
            }
        }

    }

#ifdef PARMCB_HAVE_TBB
    void set_global_tbb_concurrency(const std::size_t hardware_concurrency_hint) {
#if TBB_VERSION_MAJOR > 2020
    	oneapi::tbb::global_control global_limit(oneapi::tbb::global_control::max_allowed_parallelism, hardware_concurrency_hint);
#else
    	tbb::task_scheduler_init init(hardware_concurrency_hint);
#endif
    }
#endif

}

#endif
