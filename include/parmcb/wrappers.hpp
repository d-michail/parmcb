#ifndef PARMCB_WRAPPERS_HPP_
#define PARMCB_WRAPPERS_HPP_

//    Copyright (C) Dimitrios Michail 2019 - 2024.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#include <iostream>

#include <boost/throw_exception.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/property_map/function_property_map.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/serialization/vector.hpp>

#include <functional>

#include <parmcb/config.hpp>
#include <parmcb/detail/multiple_edges.hpp>
#include <parmcb/detail/self_loops.hpp>

namespace parmcb {

    /**
     * Wrap an MCB algorithm with proper preprocessing and postprocessing in order to support graphs 
     * with multiple edges. 
     */
    template<typename Graph, typename WeightMap, typename CycleOutputIterator>
    auto add_multiple_edges_support(std::function<typename boost::property_traits<WeightMap>::value_type(const Graph&, const WeightMap&, CycleOutputIterator)> algo) {
        return [=](const Graph& g, const WeightMap& weight_map, CycleOutputIterator out) -> typename boost::property_traits<WeightMap>::value_type {
            parmcb::detail::MultipleEdgesReductionAlgorithm<Graph, WeightMap, CycleOutputIterator, false> wrapped_algo(g, weight_map, boost::get(boost::vertex_index, g), algo);
            return wrapped_algo(out);
        };
    }

    /**
     * Wrap an MCB algorithm with proper preprocessing and postprocessing in order to support graphs 
     * with multiple edges. This version support parallelization with TBB.
     */ 
    template<typename Graph, typename WeightMap, typename CycleOutputIterator>
    auto add_multiple_edges_support_tbb(std::function<typename boost::property_traits<WeightMap>::value_type(const Graph&, const WeightMap&, CycleOutputIterator)> algo) {
        return [=](const Graph& g, const WeightMap& weight_map, CycleOutputIterator out) -> typename boost::property_traits<WeightMap>::value_type {
            parmcb::detail::MultipleEdgesReductionAlgorithm<Graph, WeightMap, CycleOutputIterator, true> wrapped_algo(g, weight_map, boost::get(boost::vertex_index, g), algo);
            return wrapped_algo(out);
        };
    }    

    /**
     * Wrap an MCB algorithm with proper preprocessing and postprocessing in order to support graphs 
     * with self loops. 
     */
    template<typename Graph, typename WeightMap, typename CycleOutputIterator>
    auto add_self_loops_support(std::function<typename boost::property_traits<WeightMap>::value_type(const Graph&, const WeightMap&, CycleOutputIterator)> algo) {
        return [=](const Graph& g, const WeightMap& weight_map, CycleOutputIterator out) -> typename boost::property_traits<WeightMap>::value_type {
            parmcb::detail::SelfLoopsReductionAlgorithm<Graph, WeightMap, CycleOutputIterator, false> wrapped_algo(g, weight_map, boost::get(boost::vertex_index, g), algo);
            return wrapped_algo(out);
        };
    }

    /**
     * Wrap an MCB algorithm with proper preprocessing and postprocessing in order to support graphs 
     * with self loops. 
     */
    template<typename Graph, typename WeightMap, typename CycleOutputIterator>
    auto add_self_loops_support_tbb(std::function<typename boost::property_traits<WeightMap>::value_type(const Graph&, const WeightMap&, CycleOutputIterator)> algo) {
        return [=](const Graph& g, const WeightMap& weight_map, CycleOutputIterator out) -> typename boost::property_traits<WeightMap>::value_type {
            parmcb::detail::SelfLoopsReductionAlgorithm<Graph, WeightMap, CycleOutputIterator, true> wrapped_algo(g, weight_map, boost::get(boost::vertex_index, g), algo);
            return wrapped_algo(out);
        };
    }    

} // parmcb

#endif


