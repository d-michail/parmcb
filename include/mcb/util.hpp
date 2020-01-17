#ifndef LIBMCB_VERIFY_HPP_
#define LIBMCB_VERIFY_HPP_

#include <vector>
#include <set>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>

#include <mcb/forestindex.hpp>

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

}

#endif
