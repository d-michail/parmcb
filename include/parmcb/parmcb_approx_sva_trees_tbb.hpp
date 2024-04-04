#ifndef PARMCB_APPROX_SVA_TREES_TBB_HPP_
#define PARMCB_APPROX_SVA_TREES_TBB_HPP_

//    Copyright (C) Dimitrios Michail 2019 - 2023.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#include <parmcb/detail/approx_spanner.hpp>
#include <parmcb/parmcb_sva_trees.hpp>

namespace parmcb {

template<class Graph, class WeightMap, class CycleOutputIterator>
typename boost::property_traits<WeightMap>::value_type approx_mcb_sva_fvs_trees_tbb(
        const Graph &g, const WeightMap &weight, std::size_t k,
        CycleOutputIterator out) {
    auto exact_algo = parmcb::mcb_sva_fvs_trees<Graph, WeightMap, CycleOutputIterator>;
    parmcb::detail::BaseApproxSpannerAlgorithm<Graph, WeightMap, CycleOutputIterator, true> algo(g, weight, boost::get(boost::vertex_index, g), k, exact_algo);
    return algo.run(out);
}

template<class Graph, class WeightMap, class CycleOutputIterator>
typename boost::property_traits<WeightMap>::value_type approx_mcb_sva_iso_trees_tbb(
        const Graph &g, const WeightMap &weight, std::size_t k,
        CycleOutputIterator out) {
    auto exact_algo = parmcb::mcb_sva_iso_trees<Graph, WeightMap, CycleOutputIterator>;
    parmcb::detail::BaseApproxSpannerAlgorithm<Graph, WeightMap, CycleOutputIterator, true> algo(g, weight, boost::get(boost::vertex_index, g), k, exact_algo);
    return algo.run(out);
}

} // parmcb

#endif
