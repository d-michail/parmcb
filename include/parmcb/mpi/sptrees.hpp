#ifndef PARMCB_SPTREES_MPI_HPP_
#define PARMCB_SPTREES_MPI_HPP_

//    Copyright (C) Dimitrios Michail 2019 - 2020.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

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
