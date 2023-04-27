#ifndef PARMCB_FP_HPP_
#define PARMCB_FP_HPP_

//    Copyright (C) Dimitrios Michail 2019 - 2023.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#include <boost/multiprecision/cpp_int.hpp>

namespace parmcb {

template<class T>
class fp {

public:
    // extended euclidean gcd algorithm
    static T ext_gcd(T &a, T &b, T &x, T &y);
    static T get_mult_inverse(T &a, T &p);

};

// extended euclidean gcd algorithm
template<class T>
T fp<T>::ext_gcd(T &a, T &b, T &x, T &y) {

    // initialize
    T _x[2], _y[2], _a[2], q;
    bool aneg, bneg;

    _x[0] = 1;
    _x[1] = 0;
    _y[0] = 0;
    _y[1] = 1;
    aneg = a < 0;
    bneg = b < 0;

    a = (a < 0) ? -a : a;
    b = (b < 0) ? -b : b;
    if (a == 0) {
        y = bneg ? -1 : 1;
        return b;
    }
    if (b == 0) {
        x = bneg ? -1 : 1;
        return a;
    }

    // swap arguments appropriately
    _a[0] = a;
    _a[1] = b;
    bool swap = false;
    if (b > a) {
        _a[0] = b;
        _a[1] = a;
        swap = true;
    }

    // do the work
    std::size_t i = 0;
    while (true) {
        q = _a[i] / _a[1 - i];
        if (_a[i] % _a[1 - i] == 0)
            break;
        _a[i] = _a[i] % _a[1 - i];
        _x[i] = _x[i] - q * _x[1 - i];
        _y[i] = _y[i] - q * _y[1 - i];
        i = 1 - i;
    }

    // did we swap arguments?
    if (swap) {
        x = _y[1 - i] * (aneg ? -1 : 1);
        y = _x[1 - i] * (bneg ? -1 : 1);
    } else {
        x = _x[1 - i] * (aneg ? -1 : 1);
        y = _y[1 - i] * (bneg ? -1 : 1);
    }

#ifdef PARMCB_INVARIANTS_CHECK
    assert(_a[1 - i] == ((aneg) ? (-a) : (a)) * x + ((bneg) ? (-b) : (b)) * y);
#endif

    return _a[1 - i];
}

// compute multiplication inverse of an element
template<class T>
T fp<T>::get_mult_inverse(T &a, T &p) {
#ifdef PARMCB_INVARIANTS_CHECK
    if ( p <= 0 )
        throw new std::runtime_error("p is not positive");
#endif

    T x, y;
    if (fp<T>::ext_gcd(a, p, x, y) != 1) {
        throw new std::runtime_error("mult inverse does not exist");
    }
    return x;
}

} // parmcb

#endif
