#ifndef LIBMCB_SPTREES_MPI_HPP_
#define LIBMCB_SPTREES_MPI_HPP_

#include <mcb/sptrees.hpp>

namespace boost {
    namespace mpi {

        template<class Graph, class WeightMap>
        struct is_commutative<mcb::SerializableMinOddCycleMinOp<Graph, WeightMap>,
                mcb::SerializableMinOddCycle<Graph, WeightMap>> : mpl::true_ {
        };

        template<class Graph>
        struct is_mpi_datatype<mcb::SerializableCandidateCycle<Graph>> : mpl::true_ {
        };

    }
}

#endif
