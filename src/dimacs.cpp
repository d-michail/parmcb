#include <iostream>
#include <fstream>
#include <map>
#include <list>

#include <boost/config.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>
#include <boost/thread.hpp>

#include <mcb/mcb.hpp>
#include <mcb/util.hpp>

using namespace boost;
namespace po = boost::program_options;

#define BUFFER_SIZE 1024
#define USAGE "Computes the minimum cycle basis of a weighted undirected graph given in DIMACS format."

int main(int argc, char *argv[]) {

    po::variables_map vm;
    try {
        po::options_description desc(
                "Computes the minimum cycle basis of a weighted undirected graph given in DIMACS format");
        // @formatter:off
        desc.add_options()
                ("help,h", "Help")
                ("verbose,v", po::value<bool>()->default_value(false)->implicit_value(true), "Verbose")
                ("signed,s", po::value<bool>()->default_value(true), "Use one of either (a) signed graph or (b) shortest path trees")
                ("parallel,p", po::value<bool>()->default_value(true), "Use parallelization")
                ("printcycles", po::value<bool>()->default_value(false)->implicit_value(true), "Print cycles")
                ("cores", po::value<int>()->default_value(0), "Number of cores")
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
    typedef graph_traits<graph_t>::vertex_descriptor vertex_descriptor;
    typedef graph_traits<graph_t>::edge_descriptor edge_descriptor;

    // create graph
    graph_t graph;
    boost::property_map<graph_t, boost::edge_weight_t>::type weight = get(boost::edge_weight, graph);

    FILE *fp = fopen(vm["input-file"].as<std::string>().c_str(), "r");
    if (fp == NULL) {
        std::cerr << "Failed to open input file." << std::endl;
        exit(EXIT_FAILURE);
    }

    char buffer[BUFFER_SIZE], problem[BUFFER_SIZE];

    std::size_t nnodes, nedges;
    std::map<std::size_t, vertex_descriptor> vertex_map;

    while (fgets(buffer, sizeof(buffer), fp) != NULL) {
        buffer[strlen(buffer) - 1] = '\0'; // eat the newline
        if (buffer[0] == 'c' || buffer[0] == '#') {
            continue;
        } else if (buffer[0] == 'p') {
            sscanf(buffer, "p %s %lu %lu", problem, &nnodes, &nedges);
            for (std::size_t i = 1; i <= nnodes; i++) {
                vertex_map[i] = add_vertex(graph);
            }
        } else if (buffer[0] == 'a' || buffer[0] == 'e') {
            int rs, rt;
            double rw = 1;
            char fc;
            sscanf(buffer, "%c %d %d %lf", &fc, &rs, &rt, &rw);

            if (vertex_map.find(rs) == vertex_map.end()) {
                std::cerr << "File contains invalid vertex number " << rs << ". Bailing out.." << std::endl;
                exit(EXIT_FAILURE);
            }
            if (vertex_map.find(rt) == vertex_map.end()) {
                std::cerr << "File contains invalid vertex number " << rt << ". Bailing out.." << std::endl;
                exit(EXIT_FAILURE);
            }
            vertex_descriptor sDescriptor = vertex(vertex_map[rs], graph);
            vertex_descriptor tDescriptor = vertex(vertex_map[rt], graph);
            edge_descriptor eDescriptor = add_edge(sDescriptor, tDescriptor, graph).first;
            weight[eDescriptor] = rw;
        }
    }

    fclose(fp);

    std::cout << "Graph has " << num_vertices(graph) << " vertices" << std::endl;
    std::cout << "Graph has " << num_edges(graph) << " edges" << std::endl;
    std::cout << std::flush;

    std::size_t cores = 0;
    if (vm.count("cores")) {
        cores = vm["cores"].as<int>();
    }
    if (cores == 0) {
        cores = boost::thread::hardware_concurrency();
    }

    boost::timer::cpu_timer timer;

    std::list<std::list<edge_descriptor>> cycles;
    double mcb_weight;
    if (vm["signed"].as<bool>()) {
        if (vm["parallel"].as<bool>()) {
            std::cout << "Using PAR_MCB_SVA_SIGNED" << std::endl;
            mcb_weight = mcb::mcb_sva_signed_parallel(graph, weight, std::back_inserter(cycles), cores);
        } else {
            std::cout << "Using MCB_SVA_SIGNED" << std::endl;
            mcb_weight = mcb::mcb_sva_signed(graph, weight, std::back_inserter(cycles));
        }
    } else {
        if (vm["parallel"].as<bool>()) {
            std::cout << "Using PAR_MCB_SVA_TREES" << std::endl;
            mcb_weight = mcb::mcb_sva_trees_parallel(graph, weight, std::back_inserter(cycles));
        } else {
            std::cout << "Using MCB_SVA_TREES" << std::endl;
            mcb_weight = mcb::mcb_sva_trees(graph, weight, std::back_inserter(cycles));
        }
    }

    timer.stop();

    if (vm["printcycles"].as<bool>()) {
        for (auto it = cycles.begin(); it != cycles.end(); it++) {
            auto cycle = *it;

            for (auto eit = cycle.begin(); eit != cycle.end(); eit++) {
                std::cout << " " << *eit;
            }
            std::cout << std::endl;

            assert(mcb::is_cycle(graph, cycle));
        }
    }

    std::cout << "weight= " << mcb_weight << std::endl;

    if (vm["verbose"].as<bool>()) {
        std::cout << "time:" << timer.format();
    }

    return EXIT_SUCCESS;
}

