#ifndef PARMCB_DETAIL_DIR_CYCLES_DIJKSTRA_HPP_
#define PARMCB_DETAIL_DIR_CYCLES_DIJKSTRA_HPP_

//    Copyright (C) Dimitrios Michail 2019 - 2023.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#include <parmcb/config.hpp>
#include <parmcb/fp.hpp>
#include <parmcb/spvecfp.hpp>
#include <parmcb/forestindex.hpp>
#include <parmcb/detail/dijkstra.hpp>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/adjacency_list.hpp>

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
    typedef typename boost::graph_traits<Graph>::vertex_iterator VertexIt;
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
    	std::cout << "Computing shortest cycle with non-zero mod p with support = " << support << std::endl;
    	boost::multiprecision::cpp_int minus_one(-1);

    	std::map<std::pair<Vertex,Vertex>, WeightType> first;
    	std::map<std::pair<Vertex,Vertex>, std::list<std::pair<bool,Edge>>> paths;
    	std::map<std::pair<Vertex,Vertex>, SpVecFP<P>> paths_vec;
    	std::map<std::pair<Vertex,Vertex>, P> residuals;
    	//std::map<std::pair<Vertex,Vertex>, WeightType> second;


    	VertexIt ui, uiend;
		for (boost::tie(ui, uiend) = boost::vertices(_g); ui != uiend; ++ui) {
			auto u = *ui;
			std::cout << "Running shortest paths from " << u << std::endl;

	    	// Compute one shortest path between each pair of vertices. Use the
	    	// undirected graph. Directions are only used in order to compute residue classes.
	        std::vector<WeightType> dist(boost::num_vertices(_g), (std::numeric_limits<WeightType>::max)());
	        boost::function_property_map<
	                parmcb::detail::VertexIndexFunctor<Graph, WeightType>,
	                Vertex, WeightType&> dist_map(
	                parmcb::detail::VertexIndexFunctor<Graph, WeightType>(dist,
	                        _index_map));
	    	std::vector<std::tuple<bool, Edge>> pred(boost::num_vertices(_g), std::make_tuple(false, Edge()));
	        boost::function_property_map<
	                parmcb::detail::VertexIndexFunctor<Graph,
	                        std::tuple<bool, Edge>>, Vertex,
	                std::tuple<bool, Edge>&> pred_map(
	                parmcb::detail::VertexIndexFunctor<Graph,
	                        std::tuple<bool, Edge> >(pred, _index_map));

			parmcb::as_undirected_dijkstra(_g, _weight_map, u, dist_map, pred_map);

			VertexIt vi, viend;
			for (boost::tie(vi, viend) = boost::vertices(_g); vi != viend; ++vi) {
				auto v = *vi;
				if (v == u || !std::get<0>(pred_map[v])) {
					continue;
				}

				std::pair<Vertex,Vertex> uv = std::make_pair(u, v);
				first[uv] = dist[v];

				std::cout << "Found shortest path from " << u << " to " << v << ", f=" << first[uv] << std::endl;

				// compute path and residual class
				SpVecFP<P> path_vec(_p);
				std::list<std::pair<bool,Edge>> path;

				Vertex path_cur = v;
				Vertex path_goal = u;
				std::cout << "Building inverse path from " << path_cur << " to " << path_goal << std::endl;
		        while (path_cur != path_goal) {
		        	std::cout << "cur=" << path_cur << std::endl;
		        	if (!std::get<0>(pred[path_cur])) {
						throw new std::runtime_error("No path found!");
		        	}
		        	Edge e = std::get<1>(pred[path_cur]);
		        	std::size_t e_index = _forest_index(e);
//		        	std::cout << "Found edge " << e << " with index " << e_index << std::endl;
		        	SpVecFP<P> e_vec(_p, e_index);

		        	bool direction = boost::target(e, _g) == path_cur;
//		        	std::cout << "Direction = " << direction << std::endl;
		        	if (!direction) {
		        		path_vec += e_vec * minus_one;
		        	} else {
		        		path_vec += e_vec;
		        	}
//		        	std::cout << "Current path vec = " << path_vec << std::endl;
		        	path.push_front(std::make_pair(direction, e));

		        	if (direction) {
		        		path_cur = boost::source(e, _g);
		        	} else {
		        		path_cur = boost::target(e, _g);
		        	}
		        }

				P residual = path_vec * support;
				residuals[uv] = residual;
				std::cout << "Residual r=" << residual << std::endl;

				paths_vec[uv] = path_vec;
				paths[uv] = path;
				//std::cout << "Path p=" << path[uv] << std::endl;
				std::cout << "Path p=";
				std::for_each(path.begin(), path.end(),[](const std::pair<bool,Edge> &x){
					std::cout << std::get<1>(x);
					if (!std::get<0>(x)) {
						std::cout << " (rev)";
					}
					std::cout << ", ";
				});
				std::cout << std::endl;
				std::cout << "Vec path p=" << path_vec << std::endl;


			}

		}

    	// TODO:



    	// TODO
		//
		// For each edge u-v combine it with s_{v-u} to get a cycle. If the cycle has non-zero
		// residue class use it. The minimum of all these usable cycles, is the one we are looking for.

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
