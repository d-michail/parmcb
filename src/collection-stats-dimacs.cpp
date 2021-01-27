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
#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>
#include <boost/thread.hpp>

#include <parmcb/parmcb.hpp>
#include <parmcb/util.hpp>
#include <parmcb/detail/cycles.hpp>
#include <parmcb/sptrees.hpp>

using namespace boost;
namespace po = boost::program_options;

#define USAGE "Computes statistics on FVS and ISO cycle collections given a weighted undirected graph in DIMACS format."


template<class Graph, class WeightMap>
class CandidateCycleBuilder {
public:
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
    typedef typename boost::property_traits<WeightMap>::value_type WeightType;

    CandidateCycleBuilder(const Graph &g, const WeightMap &weight_map) :
            g(g), weight_map(weight_map) {
    }

    std::tuple<std::set<Edge>, WeightType> operator()(const std::vector<parmcb::SPTree<Graph, WeightMap>> &trees,
            const parmcb::CandidateCycle<Graph, WeightMap> &c) const {

        std::shared_ptr<parmcb::SPNode<Graph, WeightMap>> v = trees[c.tree()].node(boost::source(c.edge(), g));
        std::shared_ptr<parmcb::SPNode<Graph, WeightMap>> u = trees[c.tree()].node(boost::target(c.edge(), g));

        Edge e = c.edge();
        bool valid = true;
        WeightType cycle_weight = boost::get(weight_map, e);
        std::set<Edge> result;
        result.insert(e);

        // first part
        Vertex w = boost::source(c.edge(), g);
        std::shared_ptr<parmcb::SPNode<Graph, WeightMap>> ws = trees[c.tree()].node(w);
        while (ws->has_pred()) {
            Edge a = ws->pred();
            if (result.insert(a).second == false) {
                valid = false;
                break;
            }
            cycle_weight += boost::get(weight_map, a);
            w = boost::opposite(a, w, g);
            ws = trees[c.tree()].node(w);
        }

        if (!valid) {
            return std::make_tuple(std::set<Edge> { }, 0.0);
        }

        // second part
        w = boost::target(c.edge(), g);
        ws = trees[c.tree()].node(w);
        while (ws->has_pred()) {
            Edge a = ws->pred();
            if (result.insert(a).second == false) {
                valid = false;
                break;
            }
            cycle_weight += boost::get(weight_map, a);
            w = boost::opposite(a, w, g);
            ws = trees[c.tree()].node(w);
        }

        if (!valid) {
            return std::make_tuple(std::set<Edge> { }, 0.0);
        }

        return std::make_tuple(result, cycle_weight);
    }

private:
    const Graph &g;
    const WeightMap &weight_map;
};



template<class Graph, class WeightMap, class CyclesBuilder>
void print_tree_stats(const char *name, const Graph &g, WeightMap weight_map) {
    std::vector<parmcb::SPTree<Graph, WeightMap>> trees;
    std::vector<parmcb::CandidateCycle<Graph, WeightMap>> cycles;
    CyclesBuilder cycles_builder;
    cycles_builder(g, weight_map, trees, cycles);
    std::cout << name << " cycles: " << cycles.size() << std::endl;

    CandidateCycleBuilder<Graph, WeightMap> candidate_cycle_builder(g, weight_map);

    typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
    typedef typename boost::property_traits<WeightMap>::value_type WeightType;

    double total_weight = 0.0;
    long total_card = 0;
    for (parmcb::CandidateCycle<Graph, WeightMap> c : cycles) {
        std::tuple<std::set<Edge>, WeightType> cc = candidate_cycle_builder(trees, c);

        total_weight += std::get<1>(cc);
        total_card += std::get<0>(cc).size();
    }

    std::cout << name << " total cycles weight: " << total_weight << std::endl;
    std::cout << name << " total cycles cardinality: " << total_card << std::endl;
}

template <class Graph, class WeightMap>
void print_all_tree_stats(const Graph& g, WeightMap weight_map) {
    print_tree_stats<Graph, WeightMap, parmcb::detail::FVSCyclesBuilder<Graph, WeightMap>>("FVS", g, weight_map);
    print_tree_stats<Graph, WeightMap, parmcb::detail::ISOCyclesBuilder<Graph, WeightMap>>("ISO", g, weight_map);
    print_tree_stats<Graph, WeightMap, parmcb::detail::HortonCyclesBuilder<Graph, WeightMap>>("HORTON", g, weight_map);
}

int main(int argc, char *argv[]) {

    po::variables_map vm;
    try {
        po::options_description desc(USAGE);
        // @formatter:off
        desc.add_options()
                ("help,h", "Help")
                ("verbose,v", po::value<bool>()->default_value(false)->implicit_value(true), "Verbose")
                ("input-file,I",po::value<std::string>(), "Input filename");
        // @formatter:on
        po::positional_options_description pos_desc;
        pos_desc.add("input-file", 1);
        po::command_line_parser parser { argc, argv };
        parser.options(desc).positional(pos_desc).allow_unregistered();
        po::parsed_options parsed_options = parser.run();
        po::store(parsed_options, vm);

        if (vm.count("help")) {
            std::cout << desc << std::endl;
            return EXIT_SUCCESS;
        }
        po::notify(vm);

        if (!vm.count("input-file")) {
            std::cerr << "Input file missing. See usage by calling with -h ." << std::endl;
            exit(EXIT_FAILURE);
        }

    } catch (const po::error &ex) {
        std::cerr << "Invalid arguments:" << ex.what() << std::endl;
        return EXIT_FAILURE;
    }

    typedef adjacency_list<vecS, vecS, undirectedS, no_property, property<edge_weight_t, double> > graph_t;

    // create graph
    graph_t graph;
    FILE *fp = fopen(vm["input-file"].as<std::string>().c_str(), "r");
    if (fp == NULL) {
        std::cerr << "Failed to open input file." << std::endl;
        exit(EXIT_FAILURE);
    }
    parmcb::read_dimacs_from_file(fp, graph);
    fclose(fp);

    if (parmcb::has_loops(graph)) {
        std::cerr << "Graph has loops, aborting.." << std::endl;
        return EXIT_FAILURE;
    }
    if (parmcb::has_multiple_edges(graph)) {
        std::cerr << "Graph has multiple edges, aborting.." << std::endl;
        return EXIT_FAILURE;
    }
    if (parmcb::has_non_positive_weights(graph, get(boost::edge_weight, graph))) {
        std::cerr << "Graph has negative or zero weight edges, aborting.." << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "graph n: " << num_vertices(graph) << std::endl;
    std::cout << "graph m: " << num_edges(graph) << std::endl;
    std::cout << "graph n*m: " << num_vertices(graph)*num_edges(graph) << std::endl;
    std::cout << std::flush;

    print_all_tree_stats(graph, get(boost::edge_weight, graph));

    return EXIT_SUCCESS;
}

