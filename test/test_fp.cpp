//    Copyright (C) Dimitrios Michail 2019 - 2021.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <boost/multiprecision/cpp_int.hpp>
#include <parmcb/fp.hpp>

TEST_CASE("mult inverse")
{
    boost::multiprecision::cpp_int p = 17;
    CHECK(p == 17);

    boost::multiprecision::cpp_int a1 = 3;
    CHECK(parmcb::fp<boost::multiprecision::cpp_int>::get_mult_inverse(a1, p) == 6);

    boost::multiprecision::cpp_int a2 = 5;
    CHECK(parmcb::fp<boost::multiprecision::cpp_int>::get_mult_inverse(a2, p) == 7);

    boost::multiprecision::cpp_int a3 = 13;
    boost::multiprecision::cpp_int p3("911048271131448098930606524753");
    boost::multiprecision::cpp_int inv3 = parmcb::fp<boost::multiprecision::cpp_int>::get_mult_inverse(a3, p3);
    while (inv3 < 0) {
        inv3 += p3;
    }
    CHECK((a3 * inv3) % p3 == 1);

}

TEST_CASE("is_prime")
{
    boost::multiprecision::cpp_int p = 17;
    CHECK(parmcb::primes<boost::multiprecision::cpp_int>::is_prime(p));

    boost::multiprecision::cpp_int p1("7901");
    CHECK(parmcb::primes<boost::multiprecision::cpp_int>::is_prime(p1));

    boost::multiprecision::cpp_int p2("7902");
    CHECK(!parmcb::primes<boost::multiprecision::cpp_int>::is_prime(p2));
}
