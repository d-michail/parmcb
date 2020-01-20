#ifndef LIBMCB_MCB_SVA_TREES_MPI_HPP_
#define LIBMCB_MCB_SVA_TREES_MPI_HPP_

#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/tuple/detail/tuple_basic.hpp>
#include <boost/timer/timer.hpp>

#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include <cstddef>
#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <set>
#include <vector>

#include <mcb/forestindex.hpp>
#include <mcb/spvecgf2.hpp>
#include <mcb/sptrees.hpp>
#include <mcb/util.hpp>

namespace mcb {

    template<class Graph, class WeightMap, class CycleOutputIterator>
    typename boost::property_traits<WeightMap>::value_type mcb_sva_trees_mpi(const Graph &g, WeightMap weight_map,
            boost::mpi::communicator& comm,
            CycleOutputIterator out) {

        typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
        typedef typename boost::property_traits<WeightMap>::value_type WeightType;

        /*
         * Index the graph
         */
        ForestIndex<Graph> forest_index(g);
        auto csd = forest_index.cycle_space_dimension();
        if (comm.rank() == 0) {
            std::cout << "Cycle space dimension: " << csd << std::endl;
        }

        /*
         * Initialize support vectors
         */
        std::vector<SpVecGF2<std::size_t>> support;
        if (comm.rank() == 0) {
            for (std::size_t k = 0; k < csd; k++) {
                support.emplace_back(k);
            }
        }

        // TODO

        return 0;
    }

} // namespace mcb

#endif
