#ifndef PARMCB_SPTREES_MPI_HPP_
#define PARMCB_SPTREES_MPI_HPP_

#include <parmcb/sptrees.hpp>

namespace boost {
    namespace mpi {

        template<class Graph, class WeightMap>
        struct is_commutative<parmcb::SerializableMinOddCycleMinOp<Graph, WeightMap>,
                parmcb::SerializableMinOddCycle<Graph, WeightMap>> : mpl::true_ {
        };

        template<class Graph>
        struct is_mpi_datatype<parmcb::SerializableCandidateCycle<Graph>> : mpl::true_ {
        };

    }
}

#endif
