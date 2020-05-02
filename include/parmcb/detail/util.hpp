#ifndef PARMCB_DETAIL_UTIL_HPP_
#define PARMCB_DETAIL_UTIL_HPP_

//    Copyright (C) Dimitrios Michail 2019 - 2020.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include <set>
#include <map>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>

namespace parmcb {

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

} // mcb

#endif
