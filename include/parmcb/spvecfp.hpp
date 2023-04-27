#ifndef PARMCB_SPVECFP_HPP_
#define PARMCB_SPVECFP_HPP_

//    Copyright (C) Dimitrios Michail 2019 - 2023.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)



#include <vector>
#include <set>
#include <iostream>
#include <cassert>

#include <boost/tuple/tuple.hpp>
#include <boost/serialization/vector.hpp>

namespace parmcb {

    template<typename P>
    class SpVecFP {

    public:
        typedef typename boost::tuple<std::size_t, P> entry_type;
        typedef typename std::vector<entry_type>::size_type size_type;
        typedef typename std::vector<entry_type>::const_iterator const_iterator;

        SpVecFP(): p(3) {
        }

        SpVecFP(const P &p): p(p) {
        }

        SpVecFP(const SpVecFP<P> &v) {
            entries = v.entries;
            p = v.p;
        }

        SpVecFP(const SpVecFP<P> &&v) {
            entries = std::move(v.entries);
            p = v.p;
        }

        ~SpVecFP(void) {
        }

        SpVecFP<P>& operator=(const SpVecFP<P> &v) {
            if (this == &v) {
                return *this;
            }
            entries = v.entries;
            p = v.p;
            return *this;
        }

        SpVecFP<P>& operator=(const SpVecFP<P> &&v) {
            if (this == &v) {
                return *this;
            }
            entries = std::move(v.entries);
            p = v.p;
            return *this;
        }

        P operator*(const SpVecFP<P> &v) const {
            P res = 0;

            auto it = entries.begin(), it_e = entries.end();
            auto v_it = v.entries.begin(), v_it_e = v.entries.end();

            while (it != it_e && v_it != v_it_e) {
                entry_type entry = *it;
                std::size_t index = std::get<0>(entry);

                entry_type v_entry = *v_it;
                std::size_t v_index = std::get<0>(v_entry);

                if (index > v_index) {
                    v_it++;
                } else if (index < v_index) {
                    it++;
                } else {
                    P v = (std::get<1>(entry) * std::get<1>(v_entry)) % p;
                    res = ( res + v ) %p;
                    it++;
                    v_it++;
                }
            }
            return res;
        }

        SpVecFP<P> operator+(const SpVecFP<P> &v) const {
            SpVecFP<P> res(p);
            auto it = entries.begin(), it_e = entries.end();
            auto v_it = v.entries.begin(), v_it_e = v.entries.end();

            // now add them
            while (it != it_e && v_it != v_it_e) {
                entry_type entry = *it;
                std::size_t index = std::get<0>(entry);
                P value= std::get<1>(entry);

                entry_type v_entry = *v_it;
                std::size_t v_index = std::get<0>(v_entry);
                P v_value = std::get<1>(v_entry);

                if (index > v_index) {
                    res.entries.push_back(boost::make_tuple(v_index, v_value));
                    v_it++;
                } else if (index < v_index) {
                    res.entries.push_back(boost::make_tuple(index, value));
                    it++;
                } else {
                    P v = (value + v_value) % p;
                    while(v < 0) v += p;   // make [-i]_p = [p-i]_p
                    while (v >= p) v -= p; // make [i+p]_p = [i]_p
                    if (v != 0) {
                        res.entries.push_back(boost::make_tuple(index, v));
                    }
                    it++;
                    v_it++;
                }
            }

            // append remaining stuff
            while (it != it_e) {
                res.entries.push_back(*it);
                it++;
            }
            while (v_it != v_it_e) {
                res.entries.push_back(*v_it);
                v_it++;
            }

            return res;
        }

        SpVecFP<P>& operator+=(const SpVecFP<P> &v) {
            *this = *this + v;
            return *this;
        }

        size_type size() const {
            return entries.size();
        }

        const_iterator begin() const
        {
            return entries.begin();
        }

        const_iterator end() const
        {
            return entries.end();
        }

        void clear() {
            entries.clear();
        }

    private:
        friend class boost::serialization::access;

        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
          ar & p;
          ar & entries;
        }

        std::vector<entry_type> entries;
        P p;
    };

    template<typename P>
    std::ostream& operator<<(std::ostream &o, const SpVecFP<P> &v) {
        if (v.size() > 0)
            o << "(" << v.size() << ") ";
        std::copy(v.begin(), v.end(), std::ostream_iterator<P> { o, " " });
        o << " (mod " << v.p << ")";
        return o;
    }

} // parmcb

#endif
