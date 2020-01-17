#ifndef LIBMCB_FVS_HPP_
#define LIBMCB_FVS_HPP_

#include <iostream>

#include <boost/scoped_array.hpp>
#include <boost/throw_exception.hpp>
#include <boost/functional/hash.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/property_map/function_property_map.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/heap/pairing_heap.hpp>

namespace mcb {

    namespace detail {

        template<class Vertex>
        struct LessVertex {
            LessVertex(std::map<Vertex, double> &priority) :
                    priority(priority) {
            }

            bool operator()(const Vertex &v, const Vertex &u) const {
                return priority[v] >= priority[u];
            }

            std::map<Vertex, double> &priority;
        };

    } // detail

    template<class Graph, class VertexOutputIterator>
    void greedy_fvs(const Graph &g, VertexOutputIterator out) {

        typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
        typedef typename boost::graph_traits<Graph>::vertex_iterator VertexIt;
        typedef typename boost::property_map<Graph, boost::vertex_index_t>::type VertexIndexMapType;
        typedef typename boost::heap::pairing_heap<Vertex, boost::heap::compare<detail::LessVertex<Vertex>>>::handle_type HeapHandleType;

        std::size_t n = boost::num_vertices(g);
        std::vector<bool> exists(n);
        std::vector<std::size_t> degree(n);
        std::vector<HeapHandleType> handle(n);
        std::map<Vertex, double> priority;
        boost::heap::pairing_heap<Vertex, boost::heap::compare<detail::LessVertex<Vertex>>> heap(
                detail::LessVertex<Vertex> { priority });

        const VertexIndexMapType &index_map = boost::get(boost::vertex_index, g);
        boost::associative_property_map<std::map<Vertex, double>> priority_map(priority);

        // initialize
        std::deque<Vertex> forRemoval;
        VertexIt vi, viend;
        for (boost::tie(vi, viend) = boost::vertices(g); vi != viend; ++vi) {
            auto v = *vi;
            auto vindex = index_map[v];
            auto d = boost::out_degree(v, g);
            exists[vindex] = true;
            if ((degree[vindex] = d) <= 1) {
                forRemoval.push_front(v);
            } else {
                priority[v] = 1.0 / d;
            }
        }

        // cleanup
        // repeatedly remove degree 0 or 1
        while (!forRemoval.empty()) {
            Vertex u = forRemoval.front();
            forRemoval.pop_front();
            auto uindex = index_map[u];
            exists[uindex] = false;

            auto eiRange = boost::out_edges(u, g);
            for (auto ei = eiRange.first; ei != eiRange.second; ++ei) {
                auto w = boost::target(*ei, g);
                auto windex = index_map[w];
                if (!exists[windex]) {
                    continue;
                }
                degree[windex]--;
                if (degree[windex] <= 1) {
                    // collect for removal
                    forRemoval.push_front(w);
                } else {
                    priority[w] = 1.0 / degree[windex];
                }
            }
        }

        // add remaining vertices into the priority queue
        for (boost::tie(vi, viend) = boost::vertices(g); vi != viend; ++vi) {
            auto v = *vi;
            auto vindex = index_map[v];
            if (!exists[vindex]) {
                continue;
            }
            handle[vindex] = heap.push(v);
        }

        // main loop
        while (!heap.empty()) {
            auto v = heap.top();
            auto vindex = index_map[v];
            heap.pop();

            if (!exists[vindex]) {
                continue;
            }

            // add to feedback vertex set
            *out++ = v;

            // remove from graph
            exists[vindex] = false;

            auto eiRange = boost::out_edges(v, g);
            for (auto ei = eiRange.first; ei != eiRange.second; ++ei) {
                auto u = boost::target(*ei, g);
                auto uindex = index_map[u];
                if (!exists[uindex]) {
                    continue;
                }
                degree[uindex]--;
                if (degree[uindex] <= 1) {
                    // collect for removal
                    forRemoval.push_front(u);
                } else {
                    priority[u] = 1.0 / degree[uindex];
                    heap.decrease(handle[uindex]);
                }
            }

            // cleanup
            while (!forRemoval.empty()) {
                Vertex u = forRemoval.front();
                forRemoval.pop_front();
                auto uindex = index_map[u];
                exists[uindex] = false;

                eiRange = boost::out_edges(u, g);
                for (auto ei = eiRange.first; ei != eiRange.second; ++ei) {
                    auto w = boost::target(*ei, g);
                    auto windex = index_map[w];
                    if (!exists[windex]) {
                        continue;
                    }
                    degree[windex]--;
                    if (degree[windex] <= 1) {
                        // collect for removal
                        forRemoval.push_front(w);
                    } else {
                        priority[w] = 1.0 / degree[windex];
                        heap.decrease(handle[windex]);
                    }

                }

            }

        }

    }

} // mcb

#endif
