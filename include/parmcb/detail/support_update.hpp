#ifndef PARMCB_DETAIL_SUPPORT_UPDATE_HPP_
#define PARMCB_DETAIL_SUPPORT_UPDATE_HPP_

//    Copyright (C) Dimitrios Michail 2019 - 2024.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#include <parmcb/config.hpp>
#include <parmcb/spvecgf2.hpp>

// #include <functional>
// #include <numeric>

#ifdef PARMCB_HAVE_TBB
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/concurrent_vector.h>
#include <tbb/task_group.h>
#endif

namespace parmcb {

namespace detail {

template<class SupportVectorContainer, bool ParallelUsingTBB>
class SupportVectorUpdater {
public:
    SupportVectorUpdater(SupportVectorContainer &support) :
            _support(support) {
    }

    void update(std::size_t k, const std::set<std::size_t>& cyclek) {
        _update(k, cyclek);
    }

private:

    template<bool is_tbb_enabled = ParallelUsingTBB>
    void _update(std::size_t k, const std::set<std::size_t>& cyclek,
            typename std::enable_if<!is_tbb_enabled>::type* = 0) {
#ifdef PARMCB_LOGGING
        std::cout << "Updating support vectors." << std::endl;
#endif
        size_t csd = _support.size();
        for (std::size_t l = k + 1; l < csd; l++) {
                if (_support[l] * cyclek == 1) {
                        _support[l] += _support[k];
                }
        }
    }

    template<bool is_tbb_enabled = ParallelUsingTBB>
    void _update(std::size_t k, const std::set<std::size_t>& cyclek,
            typename std::enable_if<is_tbb_enabled>::type* = 0) {

#ifdef PARMCB_LOGGING
        std::cout << "Updating support vectors with TBB!" << std::endl;
#endif

        size_t csd = _support.size();
        tbb::parallel_for(tbb::blocked_range<std::size_t>(k + 1, csd),
                [&](const tbb::blocked_range<std::size_t> &r) {
                auto e = r.end();
                for (std::size_t i = r.begin(); i != e; ++i) {
                        if (_support[i] * cyclek == 1) {
                                _support[i] += _support[k];
                        }
                }
                });
    }

    SupportVectorContainer &_support;
};

} // detail

} // parmcb

#endif
