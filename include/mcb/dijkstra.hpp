#ifndef LIBMCB_DIJKSTRA_HPP_
#define LIBMCB_DIJKSTRA_HPP_

#include <iostream>

#include <boost/scoped_array.hpp>
#include <boost/throw_exception.hpp>
#include <boost/functional/hash.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/property_map/function_property_map.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/detail/d_ary_heap.hpp>

#include <mcb/util.hpp>

namespace mcb {

    namespace detail {

        template<class Graph>
        struct IndexInHeapFunctor {
            typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
            typedef typename boost::property_map<Graph, boost::vertex_index_t>::type VertexIndexMapType;

            std::vector<std::size_t> &index_in_heap;
            const VertexIndexMapType &index_map;

            IndexInHeapFunctor(std::vector<std::size_t> &index_in_heap, const VertexIndexMapType &index_map) :
                    index_in_heap(index_in_heap), index_map(index_map) {
            }

            std::size_t& operator()(const Vertex &v) const {
                return index_in_heap.at(index_map[v]);
            }
        };

    } // detail

    template<class Graph, class WeightMap, class DistanceMap, class PredecessorMap>
    void dijkstra(const Graph &g, const WeightMap &weight_map,
            const typename boost::graph_traits<Graph>::vertex_descriptor &s, DistanceMap &dist_map,
            PredecessorMap &pred_map) {

        typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
        typedef typename boost::property_map<Graph, boost::vertex_index_t>::type VertexIndexMapType;
        typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
        typedef typename boost::property_traits<WeightMap>::value_type WeightType;
        typedef typename boost::property_traits<DistanceMap>::value_type DistanceType;
        typedef boost::d_ary_heap_indirect<Vertex, 4,
                boost::function_property_map<mcb::detail::VertexIndexFunctor<Graph, std::size_t>, Vertex, std::size_t&>,
                DistanceMap, std::less<double>> VertexQueue;

        std::less<DistanceType> compare;
        mcb::detail::closed_plus<DistanceType> combine = mcb::detail::closed_plus<DistanceType>();

        const VertexIndexMapType &index_map = boost::get(boost::vertex_index, g);
        std::vector<std::size_t> index_in_heap(boost::num_vertices(g));
        boost::function_property_map<mcb::detail::VertexIndexFunctor<Graph, std::size_t>, Vertex, std::size_t&> index_in_heap_map(
                mcb::detail::VertexIndexFunctor<Graph, std::size_t>(index_in_heap, index_map));

        VertexQueue queue(dist_map, index_in_heap_map, compare);

        boost::put(dist_map, s, DistanceType());
        boost::put(pred_map, s, std::make_tuple(false, Edge()));
        queue.push(s);

        while (!queue.empty()) {
            Vertex u = queue.top();
            queue.pop();
            DistanceType d_u = boost::get(dist_map, u);

            auto eiRange = boost::out_edges(u, g);
            for (auto ei = eiRange.first; ei != eiRange.second; ++ei) {
                auto e = *ei;

                auto w = boost::target(e, g);
                if (w == u) {
                    w = boost::source(e, g);
                }
                if (w == u) {
                    // self-loop
                    continue;
                }
                if (w == s) {
                    continue;
                }

                const WeightType c = combine(d_u, get(weight_map, e));
                bool visited_w = std::get<0>(boost::get(pred_map, w));
                if (!visited_w) {
                    // first time found
                    boost::put(dist_map, w, c);
                    boost::put(pred_map, w, std::make_tuple(true, e));
                    queue.push(w);
                } else if (compare(c, boost::get(dist_map, w))) {
                    // already reached
                    boost::put(dist_map, w, c);
                    boost::put(pred_map, w, std::make_tuple(true, e));
                    queue.update(w);
                }
            }
        }
    }

} // mcb

#endif
