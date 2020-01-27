#ifndef MCB_VERIFY_HPP_
#define MCB_VERIFY_HPP_

#include <vector>
#include <set>
#include <map>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>

#include <mcb/forestindex.hpp>


#define BUFFER_SIZE 1024

namespace mcb {

    namespace detail {

        // The following version of the plus functor prevents
        // problems due to overflow at positive infinity.

        template<class T>
        struct closed_plus {
            const T inf;

            closed_plus() :
                    inf((std::numeric_limits<T>::max)()) {
            }
            closed_plus(T inf) :
                    inf(inf) {
            }

            T operator()(const T &a, const T &b) const {
                using namespace std;
                if (a == inf)
                    return inf;
                if (b == inf)
                    return inf;
                return a + b;
            }
        };

        template<class Graph, class V>
        struct VertexIndexFunctor {
            typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
            typedef typename boost::property_map<Graph, boost::vertex_index_t>::type VertexIndexMapType;

            std::vector<V> &values;
            const VertexIndexMapType &index_map;

            VertexIndexFunctor(std::vector<V> &values, const VertexIndexMapType &index_map) :
                    values(values), index_map(index_map) {
            }

            V& operator()(const Vertex &v) const {
                return values.at(index_map[v]);
            }
        };

    } // detail

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
    void read_dimacs_from_file(FILE *fp, Graph& graph) {
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

}

#endif
