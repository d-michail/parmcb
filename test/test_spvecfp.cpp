//    Copyright (C) Dimitrios Michail 2019 - 2021.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <boost/multiprecision/cpp_int.hpp>
#include <parmcb/spvecfp.hpp>

TEST_CASE("spvecfp")
{
    typedef typename parmcb::SpVecFP<boost::multiprecision::cpp_int> vector_type;
    boost::multiprecision::cpp_int p = 17;

    vector_type v1(p),v2(p),v3(p),v4(p),v5(p);
    v1 = 1;
    v2 = 2;
    v3 = 3;
    v4 = 4;
    v5 = 5;

    vector_type v2_3_4 = v2 + v3 + v4;
    v2_3_4 *= boost::multiprecision::cpp_int("13");

    vector_type v1_3_5 = v1 + v3 + v5;
    v1_3_5 *= boost::multiprecision::cpp_int("15");

    vector_type res = v1_3_5 + v2_3_4;

    CHECK(res.prime() == p);
    CHECK(boost::get<1>(*(res.begin()+0)) == 15);
    CHECK(boost::get<1>(*(res.begin()+1)) == 13);
    CHECK(boost::get<1>(*(res.begin()+2)) == 11);
    CHECK(boost::get<1>(*(res.begin()+3)) == 13);
    CHECK(boost::get<1>(*(res.begin()+4)) == 15);

}

