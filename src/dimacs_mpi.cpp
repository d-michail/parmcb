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
#include <boost/format.hpp>

#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include <mcb/mpi/mcb.hpp>

using namespace boost;
namespace po = boost::program_options;

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
                ("signed,s", po::value<bool>()->default_value(false), "Use one of either (a) signed graph or (b) shortest path trees")
                ("printcycles", po::value<bool>()->default_value(false)->implicit_value(true), "Print cycles")
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
        exit(EXIT_FAILURE);
    }

    // MPI
    boost::mpi::environment env(argc, argv, boost::mpi::threading::multiple);
    boost::mpi::communicator world;

    if (env.thread_level() < boost::mpi::threading::multiple) {
        std::cerr << "Multiple thread level unsupported, bailing out.." << std::endl;
        exit(EXIT_FAILURE);
    }

    std::cout << boost::format("processor name: %s, number of tasks: %d, rank: %d\n") % env.processor_name() % world.size() % world.rank();


    typedef adjacency_list<vecS, vecS, undirectedS, no_property, property<edge_weight_t, double> > graph_t;
    typedef graph_traits<graph_t>::edge_descriptor edge_descriptor;

    // create graph
    graph_t graph;
    FILE *fp = fopen(vm["input-file"].as<std::string>().c_str(), "r");
    if (fp == NULL) {
        std::cerr << "Failed to open input file." << std::endl;
        exit(EXIT_FAILURE);
    }
    mcb::read_dimacs_from_file(fp, graph);
    fclose(fp);

    if (world.rank() == 0) {
        std::cout << "Graph has " << num_vertices(graph) << " vertices" << std::endl;
        std::cout << "Graph has " << num_edges(graph) << " edges" << std::endl;
        std::cout << std::flush;
    }

    boost::timer::cpu_timer timer;
    std::list<std::list<edge_descriptor>> cycles;
    double mcb_weight;

    if (vm["signed"].as<bool>()) {
        mcb_weight = mcb::mcb_sva_signed_mpi(graph, get(boost::edge_weight, graph), std::back_inserter(cycles), world);
    } else {
        mcb_weight = mcb::mcb_sva_trees_mpi(graph, get(boost::edge_weight, graph), std::back_inserter(cycles), world);
    }
    timer.stop();

    if (world.rank() == 0) {
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

        std::cout << "weight = " << mcb_weight << std::endl;

        if (vm["verbose"].as<bool>()) {
            std::cout << "time:" << timer.format();
        }
    }

    return EXIT_SUCCESS;
}

