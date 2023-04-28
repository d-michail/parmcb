#ifndef PARMCB_ARITHMETIC_HPP_
#define PARMCB_ARITHMETIC_HPP_

//    Copyright (C) Dimitrios Michail 2019 - 2023.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#include <boost/multiprecision/cpp_int.hpp>

namespace parmcb {

/**
 * The prime type used for the directed cycle basis algorithm.
 */
typedef boost::multiprecision::cpp_int ptype;

// give compare for ptype inside the parmcb namespace
inline int compare(const ptype& x, const ptype& y)
{
    if (x < y) return -1;
    else if (x > y) return 1;
    else return 0;
}

}

#endif
