#ifndef PARMCB_APPROX_SVA_SIGNED_TBB_HPP_
#define PARMCB_APPROX_SVA_SIGNED_TBB_HPP_

//    Copyright (C) Dimitrios Michail 2019 - 2023.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#include <parmcb/detail/approx_spanner.hpp>
#include <parmcb/parmcb_sva_signed_tbb.hpp>

namespace parmcb {

namespace detail {

template<class Graph, class WeightMap, class CycleOutputIterator>
struct mcb_sva_signed_tbb {
    typename boost::property_traits<WeightMap>::value_type operator()(
            const Graph &g, const WeightMap &weight, CycleOutputIterator out) {
        return parmcb::mcb_sva_signed_tbb(g, weight, out);
    }
};
} // detail

template<class Graph, class WeightMap, class CycleOutputIterator>
typename boost::property_traits<WeightMap>::value_type approx_mcb_sva_signed_tbb(
        const Graph &g, const WeightMap &weight, std::size_t k,
        CycleOutputIterator out) {
    typedef typename parmcb::detail::mcb_sva_signed_tbb<Graph,WeightMap,CycleOutputIterator> ExactAlgo;

    parmcb::detail::BaseApproxSpannerAlgorithm<Graph, WeightMap, ExactAlgo, true> algo(g, weight, boost::get(boost::vertex_index, g), k);
    return algo.run(out);
}

} // parmcb

#endif
