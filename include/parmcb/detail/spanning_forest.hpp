#ifndef PARMCB_DETAIL_SPANNING_FOREST_HPP_
#define PARMCB_DETAIL_SPANNING_FOREST_HPP_

#include <map>
#include <queue>
#include <set>

#include <boost/graph/adjacency_list.hpp>

namespace parmcb {

    namespace detail {

        template<class Graph, class OutputIterator>
        std::size_t spanning_forest(const Graph &g, OutputIterator spanning_forest_edges) {

            if (boost::num_vertices(g) == 0)
                return 0;

            typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
            typedef typename boost::graph_traits<Graph>::vertex_iterator VertexIt;
            typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;

            BOOST_CONCEPT_ASSERT(( boost::VertexListGraphConcept<Graph> ));
            BOOST_CONCEPT_ASSERT(( boost::OutputIteratorConcept<OutputIterator, Edge> ));

            std::queue<Vertex> queue;
            std::unordered_set<Vertex> unreached;
            VertexIt ui, uiend;
            for (boost::tie(ui, uiend) = boost::vertices(g); ui != uiend; ++ui) {
                unreached.insert(*ui);
            }

            std::size_t c = 0;
            while (!unreached.empty()) {
                auto vi = unreached.begin();
                auto v = *vi;
                unreached.erase(vi);
                queue.push(v);

                while (!queue.empty()) {
                    auto u = queue.front();
                    queue.pop();

                    auto eiRange = boost::out_edges(u, g);
                    for (auto ei = eiRange.first; ei != eiRange.second; ++ei) {
                        auto e = *ei;
                        auto w = boost::target(e, g);
                        if (w == u) {
                            // ignore self-loop
                            continue;
                        }

                        auto wit = unreached.find(w);
                        if (wit == unreached.end()) {
                            continue;
                        }
                        unreached.erase(wit);
                        *spanning_forest_edges++ = e;
                        queue.push(w);
                    }
                }
                c++;
            }

            return c;
        }

    } // detail

} // parmcb

#endif
