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
	typedef typename boost::graph_traits<Graph>::edge_iterator EdgeIt;
	typedef typename boost::property_map<Graph, boost::vertex_index_t>::type VertexIndexMapType;
	typedef typename boost::property_traits<WeightMap>::value_type WeightType;
	typedef typename boost::property_map<Graph, boost::edge_weight_t>::type EdgeWeightMapType;

	ShortestNonZeroModpCycle(const Graph &g, const WeightMap &weight_map,
			const VertexIndexMapType &index_map,
			const ForestIndex<Graph> &forest_index, P p) :
			_g(g), _weight_map(weight_map), _index_map(index_map), _forest_index(
					forest_index), _p(p) {
	}

	SpVecFP<P> operator()(const SpVecFP<P> &support) {
		std::cout
				<< "Computing shortest cycle with non-zero mod p with support = "
				<< support << std::endl;

		std::map<std::pair<Vertex, Vertex>, WeightType> f_weights;
		std::map<std::pair<Vertex, Vertex>, SpVecFP<P>> f_paths;
		std::map<std::pair<Vertex, Vertex>, P> f_residuals;

		compute_first_paths(support, f_weights, f_paths, f_residuals);

		std::map<std::pair<Vertex, Vertex>, WeightType> s_weights;
		std::map<std::pair<Vertex, Vertex>, SpVecFP<P>> s_paths;
		std::map<std::pair<Vertex, Vertex>, P> s_residuals;

		compute_second_paths(support, f_weights, f_paths, f_residuals,
				s_weights, s_paths, s_residuals);

		return compute_min_cycle(support, s_weights, s_paths, s_residuals);
	}

private:

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

	void compute_first_paths(const SpVecFP<P> &support,
			std::map<std::pair<Vertex, Vertex>, WeightType> &weights,
			std::map<std::pair<Vertex, Vertex>, SpVecFP<P>> &paths,
			std::map<std::pair<Vertex, Vertex>, P> &residuals) {

#ifdef PARMCB_LOGGING
        std::cout << "Computing first shortest paths" std::endl;
#endif

		// Compute one shortest path between each pair of vertices. Use the undirected graph.
		// Directions are only used in order to compute residue classes.
		VertexIt ui, uiend;
		for (boost::tie(ui, uiend) = boost::vertices(_g); ui != uiend; ++ui) {
			auto u = *ui;
			std::cout << "Running shortest paths from " << u << std::endl;

			std::vector<WeightType> dist(boost::num_vertices(_g),
					(std::numeric_limits<WeightType>::max)());
			boost::function_property_map<
					parmcb::detail::VertexIndexFunctor<Graph, WeightType>,
					Vertex, WeightType&> dist_map(
					parmcb::detail::VertexIndexFunctor<Graph, WeightType>(dist,
							_index_map));
			std::vector<std::tuple<bool, Edge>> pred(boost::num_vertices(_g),
					std::make_tuple(false, Edge()));
			boost::function_property_map<
					parmcb::detail::VertexIndexFunctor<Graph,
							std::tuple<bool, Edge>>, Vertex,
					std::tuple<bool, Edge>&> pred_map(
					parmcb::detail::VertexIndexFunctor<Graph,
							std::tuple<bool, Edge> >(pred, _index_map));

			parmcb::as_undirected_dijkstra(_g, _weight_map, u, dist_map,
					pred_map);

			VertexIt vi, viend;
			for (boost::tie(vi, viend) = boost::vertices(_g); vi != viend;
					++vi) {
				auto v = *vi;
				if (v == u || !std::get<0>(pred_map[v])) {
					// skip when self-loop or not reachable
					continue;
				}

				std::pair<Vertex, Vertex> uv = std::make_pair(u, v);
				weights[uv] = dist[v];

//				std::cout << "Found shortest path from " << u << " to " << v
//						<< ", f=" << weights[uv] << std::endl;

				// compute path and residual class
				SpVecFP<P> path(_p);

				Vertex path_cur = v;
				Vertex path_goal = u;
				while (path_cur != path_goal) {
					if (!std::get<0>(pred[path_cur])) {
						throw new std::runtime_error("Invalid path returned!");
					}
					Edge e = std::get<1>(pred[path_cur]);
					std::size_t e_idx = _forest_index(e);
					bool direction = boost::target(e, _g) == path_cur;
					SpVecFP<P> e_vec(_p, e_idx);
					if (!direction) {
						path += -e_vec;
					} else {
						path += e_vec;
					}
					if (direction) {
						path_cur = boost::source(e, _g);
					} else {
						path_cur = boost::target(e, _g);
					}
				}

				residuals[uv] = path * support;
				paths[uv] = path;
			}
		}
	}

	/**
	 * Compute a second shortest path for each pair of vertices which has a different residue class
	 * from the first path.
	 */
	void compute_second_paths(const SpVecFP<P> &support,
			const std::map<std::pair<Vertex, Vertex>, WeightType> &first_weights,
			const std::map<std::pair<Vertex, Vertex>, SpVecFP<P>> &first_paths,
			const std::map<std::pair<Vertex, Vertex>, P> &first_residuals,
			std::map<std::pair<Vertex, Vertex>, WeightType> &second_weights,
			std::map<std::pair<Vertex, Vertex>, SpVecFP<P>> &second_paths,
			std::map<std::pair<Vertex, Vertex>, P> &second_residuals) {

#ifdef PARMCB_LOGGING
        std::cout << "Computing second shortest paths" std::endl;
#endif

		VertexIt ui, uiend;
		for (boost::tie(ui, uiend) = boost::vertices(_g); ui != uiend; ++ui) {
			auto u = *ui;
			std::cout << "Running modified shortest paths from " << u
					<< std::endl;

			std::vector<WeightType> dist(boost::num_vertices(_g),
					(std::numeric_limits<WeightType>::max)());
			boost::function_property_map<
					parmcb::detail::VertexIndexFunctor<Graph, WeightType>,
					Vertex, WeightType&> dist_map(
					parmcb::detail::VertexIndexFunctor<Graph, WeightType>(dist,
							_index_map));
			std::vector<std::tuple<bool, SpVecFP<P>>> path(boost::num_vertices(_g),
					std::make_tuple(false, SpVecFP<P>(_p)));
			boost::function_property_map<
					parmcb::detail::VertexIndexFunctor<Graph,
							std::tuple<bool, SpVecFP<P>>>, Vertex,
					std::tuple<bool, SpVecFP<P>>&> path_map(
					parmcb::detail::VertexIndexFunctor<Graph,
							std::tuple<bool, SpVecFP<P>> >(path, _index_map));

			second_sp_dijkstra(support, first_weights, first_paths, first_residuals, u, dist_map, path_map);

			VertexIt vi, viend;
			for (boost::tie(vi, viend) = boost::vertices(_g); vi != viend;
					++vi) {
				auto v = *vi;
				if (v == u || !std::get<0>(path_map[v])) {
					// skip when self-loop or not reachable
					continue;
				}

				std::pair<Vertex, Vertex> uv = std::make_pair(u, v);
				second_weights[uv] = dist[v];

				std::cout << "Found second shortest path from " << u << " to "
						<< v << ", f=" << second_weights[uv] << std::endl;

				auto p = std::get<1>(path_map[v]);
				second_residuals[uv] = p * support;
				second_paths[uv] = p;
			}

		}

	}

	/**
	 * Compute final minimum cycle result.
	 */
	SpVecFP<P> compute_min_cycle(const SpVecFP<P> &support,
			std::map<std::pair<Vertex, Vertex>, WeightType> &s_weights,
			std::map<std::pair<Vertex, Vertex>, SpVecFP<P>> &s_paths,
			std::map<std::pair<Vertex, Vertex>, P> &s_residuals) {

		bool min_cycle_found;
		SpVecFP<P> min_cycle;
		WeightType min_cycle_weight;

		// For each edge u-v combine it with s_{v-u} to get a cycle. If the cycle has non-zero
		// residue class use it. The minimum of all these usable cycles, is the one we are looking for.
		EdgeIt e_it, e_it_end;
		for (boost::tie(e_it, e_it_end) = boost::edges(_g); e_it != e_it_end;
				++e_it) {
			auto e = *e_it;
			auto e_weight = _weight_map[e];

			auto u = boost::source(e, _g);
			auto v = boost::target(e, _g);
			auto vu = std::make_pair(v, u);

			if (!s_residuals.count(vu)) {
				// no path found from v->u
				continue;
			}

			auto cycle = s_paths[vu] + SpVecFP<P>(_p, _forest_index(e));
			auto cycle_weight = s_weights[vu] + e_weight;
			P cycle_residue = cycle * support;

			if (cycle_residue == 0) {
				continue;
			}

			std::cout << "Candidate for minimum is=" << cycle << " with weight="
					<< cycle_weight << " and residue=" << cycle_residue
					<< std::endl;

			if (!min_cycle_found) {
				min_cycle_found = true;
				min_cycle = cycle;
				min_cycle_weight = cycle_weight;
			} else if (cycle_weight < min_cycle_weight) {
				min_cycle = cycle;
				min_cycle_weight = cycle_weight;
			}
		}

		if (!min_cycle_found) {
			throw new std::runtime_error("Failed to find minimum cycle.");
		}

		return min_cycle;
	}

	/**
	 * Compute the second shortest path.
	 *
	 * TODO: rewrite me to work correctly
	 */
	template<class DistanceMap, class PathMap>
	void second_sp_dijkstra(
			const SpVecFP<P> &support,
			const std::map<std::pair<Vertex, Vertex>, WeightType> &f_weights,
			const std::map<std::pair<Vertex, Vertex>, SpVecFP<P>> &f_paths,
			const std::map<std::pair<Vertex, Vertex>, P> &f_residuals,
			const typename boost::graph_traits<Graph>::vertex_descriptor &u,
			DistanceMap &dist_map,
			PathMap &path_map) {

		typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
		typedef typename boost::property_traits<WeightMap>::value_type WeightType;
		typedef typename boost::property_traits<DistanceMap>::value_type DistanceType;
		typedef boost::d_ary_heap_indirect<Vertex, 4,
				boost::function_property_map<
						parmcb::detail::VertexIndexFunctor<Graph, std::size_t>,
						Vertex, std::size_t&>, DistanceMap,
				std::less<DistanceType>> VertexQueue;

		typedef typename boost::graph_traits<Graph>::directed_category DirectedCat;
		bool is_directed = boost::detail::is_directed(DirectedCat());
		if (!is_directed) {
			throw new std::runtime_error("Graph is not directed.");
		}

		std::less<DistanceType> compare;
		parmcb::detail::closed_plus<DistanceType> combine =
				parmcb::detail::closed_plus<DistanceType>();

		std::vector<std::size_t> index_in_heap(boost::num_vertices(_g));
		boost::function_property_map<
				parmcb::detail::VertexIndexFunctor<Graph, std::size_t>, Vertex,
				std::size_t&> index_in_heap_map(
				parmcb::detail::VertexIndexFunctor<Graph, std::size_t>(
						index_in_heap, _index_map));

		VertexQueue queue(dist_map, index_in_heap_map, compare);

		// initialize step
		VertexIt vi, viend;
		for (boost::tie(vi, viend) = boost::vertices(_g); vi != viend;
				++vi) {
			auto v = *vi;
			std::pair<Vertex, Vertex> uv = std::make_pair(u, v);
			if (!f_residuals.count(uv)) {
				// no path found from u->v
				continue;
			}

			WeightType f_weight = f_weights.at(uv);
			SpVecFP<P> f_path = f_paths.at(uv);
			P f_residual = f_residuals.at(uv);

			auto eiOutRange = boost::out_edges(v, _g);
			auto eiInRange = boost::in_edges(v, _g);
			auto range = boost::join(
					boost::iterator_range<decltype(eiOutRange.first)>(
							eiOutRange.first, eiOutRange.second),
					boost::iterator_range<decltype(eiInRange.first)>(
							eiInRange.first, eiInRange.second));

			bool min_found = false;
			WeightType min_weight;
			SpVecFP<P> min_path;
			for (const auto &e : range) {
				auto w = boost::target(e, _g);
				bool direction = false;
				if (w == v) {
					w = boost::source(e, _g);
					direction = true;
				}
				if (w == v) {
					throw new std::runtime_error("Self loop detected");
				}

				SpVecFP<P> e_path = SpVecFP<P>(_p, (direction?1:-1) * _forest_index(e));
				auto s_path = f_path + e_path;
				auto s_residual = s_path * support;
				if (s_residual != f_residual) {
					auto e_weight = _weight_map[e];
					auto s_weight = f_weight + e_weight;
					if (!min_found) {
						min_found = true;
						min_weight = s_weight;
						min_path = s_path;
					} else if (s_weight < min_weight) {
						min_weight = s_weight;
						min_path = s_path;
					}
				}
			}

			if (min_found) {
				boost::put(dist_map, v, min_weight);
				boost::put(path_map, v, std::make_tuple(false, min_path));
				queue.push(v);
			}
		}

		while (!queue.empty()) {
			Vertex w = queue.top();
			queue.pop();
			DistanceType d_w = boost::get(dist_map, w);

			auto eiOutRange = boost::out_edges(w, _g);
			auto eiInRange = boost::in_edges(w, _g);
			auto range = boost::join(
					boost::iterator_range<decltype(eiOutRange.first)>(
							eiOutRange.first, eiOutRange.second),
					boost::iterator_range<decltype(eiInRange.first)>(
							eiInRange.first, eiInRange.second));

			// TODO
//			for (const auto &e : range) {
//				auto o = boost::target(e, _g);
//				if (w == o) {
//					o = boost::source(e, _g);
//				}
//				if (w == o) {
//					// self-loop
//					continue;
//				}
//				if (w == u) {
//					continue;
//				}
//
//				const WeightType c = combine(d_w, get(_weight_map, e));
//				bool visited_o = std::get<0>(boost::get(pred_map, o));
//				if (!visited_o) {
//					// first time found
//					boost::put(dist_map, o, c);
//					boost::put(pred_map, o, std::make_tuple(true, e));
//					queue.push(o);
//				} else if (compare(c, boost::get(dist_map, o))) {
//					// already reached
//					boost::put(dist_map, o, c);
//					boost::put(pred_map, o, std::make_tuple(true, e));
//					queue.update(o);
//				}
//			}
		}
	}

	const Graph &_g;
	const WeightMap &_weight_map;
	const VertexIndexMapType &_index_map;
	const ForestIndex<Graph> &_forest_index;
	const P _p;

};

} // detail

template<class Graph, class WeightMap, typename P>
parmcb::SpVecFP<P> shortest_non_zero_cycle_modp(const Graph &g,
		WeightMap weight_map, const ForestIndex<Graph> &forest_index, P p,
		const parmcb::SpVecFP<P> &support) {
	const auto &index_map = boost::get(boost::vertex_index, g);
	parmcb::detail::ShortestNonZeroModpCycle<Graph, WeightMap, P, false> cycler(
			g, weight_map, index_map, forest_index, p);
	return cycler(support);
}

template<class Graph, class WeightMap, typename P>
parmcb::SpVecFP<P> shortest_non_zero_cycle_modp_tbb(const Graph &g,
		WeightMap weight_map, const ForestIndex<Graph> &forest_index, P p,
		const parmcb::SpVecFP<P> &support) {
	const auto &index_map = boost::get(boost::vertex_index, g);
	parmcb::detail::ShortestNonZeroModpCycle<Graph, WeightMap, P, true> cycler(
			g, weight_map, index_map, forest_index, p);
	return cycler(support);
}

} // parmcb

#endif
