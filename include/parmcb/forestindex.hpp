#ifndef PARMCB_FORESTINDEX_HPP_
#define PARMCB_FORESTINDEX_HPP_

#include <vector>
#include <map>
#include <queue>
#include <set>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/concept/assert.hpp>

#include <parmcb/detail/spanning_forest.hpp>

namespace parmcb {

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
            k = parmcb::detail::spanning_forest(g, std::inserter(forest, forest.begin()));

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

} // namespace parmcb

#endif
