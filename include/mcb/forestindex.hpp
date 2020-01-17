#ifndef LIBMCB_FORESTINDEX_HPP_
#define LIBMCB_FORESTINDEX_HPP_

#include <vector>
#include <map>
#include <queue>
#include <set>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/concept/assert.hpp>

namespace mcb {

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

    template<class Graph>
    class ForestIndex {

    public:
        typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
        typedef typename std::size_t size_type;

        explicit ForestIndex(const Graph &g) {
            create_index(g);
        }

        ForestIndex(const ForestIndex &ei) {
            n = ei.n;
            m = ei.m;
            k = ei.k;
            index = ei.index;
            reverse_index = ei.reverse_index;
        }

        ~ForestIndex(void) {
        }

        ForestIndex& operator=(const ForestIndex &ei) {
            if (this == &ei) {
                return *this;
            }
            n = ei.n;
            m = ei.m;
            k = ei.k;
            index = ei.index;
            reverse_index = ei.reverse_index;
            return *this;
        }

        const Edge& operator()(const size_type &i) const {
            return reverse_index[i];
        }

        const size_type& operator()(const Edge &e) const {
            return index.at(e);
        }

        bool is_on_forest(const Edge &e) const {
            return index.at(e) >= cycle_space_dimension();
        }

        size_type cycle_space_dimension() const {
            return m - n + k;
        }

        size_type weak_connected_components() const {
            return k;
        }

    private:
        typedef typename boost::graph_traits<Graph>::edge_iterator EdgeIt;

        size_type n = 0;
        size_type m = 0;
        size_type k = 0;
        std::map<Edge, size_type> index;
        std::vector<Edge> reverse_index;

        void create_index(const Graph &g) {
            std::set<Edge> forest;
            n = boost::num_vertices(g);
            m = boost::num_edges(g);
            k = mcb::detail::spanning_forest(g, std::inserter(forest, forest.begin()));

            index.clear();
            reverse_index.resize(m);

            size_type csd = m - n + k; // cycle space dimension
            size_type low = 0;
            size_type high = csd;

            EdgeIt ei, eiend;
            for (boost::tie(ei, eiend) = boost::edges(g); ei != eiend; ++ei) {
                auto e = *ei;
                if (forest.find(e) == forest.end()) {
                    index[e] = low;
                    reverse_index[low] = e;
                    low++;
                } else {
                    index[e] = high;
                    reverse_index[high] = e;
                    high++;
                }
            }
        }

    };

} // namespace mcb

#endif
