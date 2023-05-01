#ifndef PARMCB_DETAIL_DIR_CYCLES_DIJKSTRA_HPP_
#define PARMCB_DETAIL_DIR_CYCLES_DIJKSTRA_HPP_

//    Copyright (C) Dimitrios Michail 2019 - 2023.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#include <parmcb/config.hpp>
#include <parmcb/fp.hpp>
#include <parmcb/spvecfp.hpp>

#include <functional>
#include <numeric>

#ifdef PARMCB_HAVE_TBB
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/concurrent_vector.h>
#include <tbb/task_group.h>
#endif

namespace parmcb {

namespace detail {

template<class Graph, class WeightMap, typename P, bool ParallelUsingTBB>
class ShortestNonZeroModpCycle {
public:
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
    typedef typename boost::property_map<Graph, boost::vertex_index_t>::type VertexIndexMapType;
    typedef typename boost::property_traits<WeightMap>::value_type WeightType;
    typedef typename boost::property_map<Graph, boost::edge_weight_t>::type EdgeWeightMapType;

    ShortestNonZeroModpCycle(const Graph &g, const WeightMap &weight_map,
            const VertexIndexMapType &index_map, const ForestIndex<Graph>& forest_index, P p) :
            _g(g), _weight_map(weight_map), _index_map(index_map), _forest_index(forest_index), _p(p) {
    }

    SpVecFP<P> operator()(const SpVecFP<P>& support) {
        return find_shortest_non_zero_cycle_mod_p(support);
    }

private:

    SpVecFP<P> find_shortest_non_zero_cycle_mod_p(const SpVecFP<P>& support) {
    	// TODO
    	return support;
    }

    template<class EdgeOutputIterator, bool is_tbb_enabled = ParallelUsingTBB>
    WeightType find_shortest_cycle(EdgeOutputIterator out,
            typename std::enable_if<!is_tbb_enabled>::type* = 0) {

    	// TODO

    }

    template<class EdgeOutputIterator, bool is_tbb_enabled = ParallelUsingTBB>
    WeightType find_shortest_cycle(EdgeOutputIterator out,
            typename std::enable_if<is_tbb_enabled>::type* = 0) {

    	// TODO

    }

    const Graph &_g;
    const WeightMap &_weight_map;
    const VertexIndexMapType &_index_map;
    const ForestIndex<Graph>& _forest_index;
    const P _p;

};

} // detail

template<class Graph, class WeightMap, typename P>
parmcb::SpVecFP<P> shortest_non_zero_cycle_modp(const Graph &g, WeightMap weight_map, const ForestIndex<Graph>& forest_index, P p, const parmcb::SpVecFP<P>& support) {
	const auto &index_map = boost::get(boost::vertex_index, g);
	parmcb::detail::ShortestNonZeroModpCycle<Graph, WeightMap, P, false> cycler(g, weight_map, index_map, forest_index, p);
	return cycler(support);
}

template<class Graph, class WeightMap, typename P>
parmcb::SpVecFP<P> shortest_non_zero_cycle_modp_tbb(const Graph &g, WeightMap weight_map, const ForestIndex<Graph>& forest_index, P p, const parmcb::SpVecFP<P>& support) {
	const auto &index_map = boost::get(boost::vertex_index, g);
	parmcb::detail::ShortestNonZeroModpCycle<Graph, WeightMap, P, true> cycler(g, weight_map, index_map, forest_index, p);
	return cycler(support);
}

} // parmcb

#endif
