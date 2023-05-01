#ifndef PARMCB_DIR_SVA_SIGNED_HPP_
#define PARMCB_DIR_SVA_SIGNED_HPP_

//    Copyright (C) Dimitrios Michail 2019 - 2021.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/tuple/detail/tuple_basic.hpp>
#include <boost/timer/timer.hpp>

#include <cstddef>
#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <set>
#include <vector>
#include <list>

#include <parmcb/config.hpp>
#include <parmcb/arithmetic.hpp>
#include <parmcb/fp.hpp>
#include <parmcb/forestindex.hpp>
#include <parmcb/spvecfp.hpp>
#include <parmcb/util.hpp>
#include <parmcb/detail/dir_cycles_dijkstra.hpp>

namespace parmcb {

template<class Graph, class WeightMap, class CycleOutputIterator>
typename boost::property_traits<WeightMap>::value_type mcb_dir_sva_signed(
		const Graph &g, WeightMap weight_map, parmcb::ptype p,
		CycleOutputIterator out) {

	typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
	typedef typename boost::property_traits<WeightMap>::value_type WeightType;

#ifdef PARMCB_INVARIANTS_CHECK
	if (parmcb::has_multiple_edges(g)) {
		throw new std::runtime_error("Multiple edges?");
	}
	if (parmcb::has_loops(g)) {
		throw new std::runtime_error("Self loops?");
	}
	if (parmcb::has_non_positive_weights(g, weight_map)) {
		throw new std::runtime_error("Non positive weights?");
	}

	if (!parmcb::primes<parmcb::ptype>::is_prime(p)) {
		throw new std::runtime_error("Parameter p is not a prime number.");
	}
#endif

	/*
	 * Index the graph
	 */
	ForestIndex<Graph> forest_index(g);
	auto csd = forest_index.cycle_space_dimension();
	if (csd <= 0) {
		return WeightType();
	}

	std::cout << "Cycle space dimension=" << csd << std::endl;

	/*
	 * Initialize support vectors
	 */
	std::vector<parmcb::SpVecFP<parmcb::ptype>> support;
	for (std::size_t k = 0; k < csd; k++) {
		support.emplace_back(p, k);
	}

	std::for_each(support.begin(), support.end(),[](const parmcb::SpVecFP<parmcb::ptype> &sup){
		std::cout << sup << std::endl;
	});

	boost::timer::cpu_timer cycle_timer;
	cycle_timer.stop();
	boost::timer::cpu_timer support_timer;
	support_timer.stop();

	/*
	 * Main loop
	 */
	WeightType mcb_weight = WeightType();
	parmcb::ptype tmpi;
	for (std::size_t k = 0; k < csd; k++) {
		/*
		 * Compute shortest cycle which is non-zero (mod p)
		 */
		cycle_timer.resume();
		parmcb::SpVecFP<parmcb::ptype> cycle_k = parmcb::shortest_non_zero_cycle_modp(g, weight_map, forest_index, p, support[k]);
		cycle_timer.stop();

		/*
		 * Update support vectors
		 */
		support_timer.resume();

        // precompute part
        // NOTE: we do not precompute inverses, since we don't want to have a dependency
		//       on the maximum size of an array that will store these values
        //       p is O(d^2 logd) and thus O(logd) to compute inverse at most d times,
		//       thus O(d logd) = O(m logm) in total
		parmcb::ptype tmpk = support[k] * cycle_k;
		while (tmpk < 0) tmpk += p;    // make [-i]_p = [p-i]_p
		while ( tmpk >= p ) tmpk -= p; // make [i+p]_p = [i]_p
		parmcb::SpVecFP<parmcb::ptype> tmp = support[k] * fp<parmcb::ptype>::get_mult_inverse( tmpk, p );

		// update support_j, j > k
		for(std::size_t j = k+1; j < csd; j++ ) {
			support[j] -=  tmp * (cycle_k * support[j]) ;
		}
		support_timer.stop();

		/*
		 * Output new cycle
		 */
		std::list<Edge> cyclek_edgelist;
		WeightType cyclek_weight = WeightType();
		auto ck_end = cycle_k.end();
		for (auto ck_it = cycle_k.begin(); ck_it != ck_end; ck_it++) {
			auto index = boost::get<0>(*ck_it);
			Edge e = forest_index(index);
			cyclek_weight += weight_map[e];
			cyclek_edgelist.push_back(e);
		}
		*out++ = cyclek_edgelist;
		mcb_weight += cyclek_weight;
	}

#ifdef PARMCB_LOGGING
        std::cout << "cycle   timer" << cycle_timer.format();
        std::cout << "support timer" << support_timer.format();
#endif

	return mcb_weight;
}

} // parmcb

#endif
