#ifndef PARMCB_SPVECGF2_HPP_
#define PARMCB_SPVECGF2_HPP_

#include <vector>
#include <set>
#include <iostream>
#include <cassert>

#include <boost/serialization/vector.hpp>

namespace parmcb {

    template<typename U>
    class SpVecGF2 {

    public:

        typedef typename std::vector<U>::size_type size_type;
        typedef typename std::vector<U>::const_iterator const_iterator;

        SpVecGF2() {
        }

        SpVecGF2(const U &i) {
            ones.push_back(i);
        }

        SpVecGF2(const SpVecGF2<U> &v) {
            ones = v.ones;
        }

        SpVecGF2(const SpVecGF2<U> &&v) {
            ones = std::move(v.ones);
        }

        SpVecGF2(const std::set<U> &v) {
            std::copy(v.begin(), v.end(), std::back_inserter(ones));
        }

        ~SpVecGF2(void) {
        }

        SpVecGF2<U>& operator=(const SpVecGF2<U> &v) {
            if (this == &v) {
                return *this;
            }
            ones = v.ones;
            return *this;
        }

        SpVecGF2<U>& operator=(const SpVecGF2<U> &&v) {
            if (this == &v) {
                return *this;
            }
            ones = std::move(v.ones);
            return *this;
        }

        int operator*(const SpVecGF2<U> &v) const {
            int res = 0;

            auto it = ones.begin(), it_e = ones.end();
            auto v_it = v.ones.begin(), v_it_e = v.ones.end();

            while (it != it_e && v_it != v_it_e) {
                U index = *it;
                U v_index = *v_it;
                if (index > v_index) {
                    v_it++;
                } else if (index < v_index) {
                    it++;
                } else {
                    res = (res + 1) % 2;
                    it++;
                    v_it++;
                }
            }
            return res;
        }

        int operator*(const std::set<U> &v) const {
            int res = 0;

            auto it = ones.begin(), it_e = ones.end();
            auto v_it = v.begin(), v_it_e = v.end();

            while (it != it_e && v_it != v_it_e) {
                U index = *it;
                U v_index = *v_it;
                if (index > v_index) {
                    v_it++;
                } else if (index < v_index) {
                    it++;
                } else {
                    res = (res + 1) % 2;
                    it++;
                    v_it++;
                }
            }
            return res;
        }

        SpVecGF2<U> operator+(const SpVecGF2<U> &v) const {
            SpVecGF2<U> res;
            auto it = ones.begin(), it_e = ones.end();
            auto v_it = v.ones.begin(), v_it_e = v.ones.end();

            // now add them
            while (it != it_e && v_it != v_it_e) {
                U index = *it;
                U v_index = *v_it;

                if (index > v_index) {
                    res.ones.push_back(v_index);
                    v_it++;
                } else if (index < v_index) {
                    res.ones.push_back(index);
                    it++;
                } else {
                    // 1 + 1 = 0, don't add anything
                    it++;
                    v_it++;
                }
            }

            // append remaining stuff
            while (it != it_e) {
                res.ones.push_back(*it);
                it++;
            }
            while (v_it != v_it_e) {
                res.ones.push_back(*v_it);
                v_it++;
            }

            return res;
        }

        SpVecGF2<U>& operator+=(const SpVecGF2<U> &v) {
            *this = *this + v;
            return *this;
        }

        void add(U pos) {
            assert(ones.empty() || pos >= *(ones.end()));

            ones.push_back(pos);
        }

        size_type size() const {
            return ones.size();
        }

        const_iterator begin() const
        {
            return ones.begin();
        }

        const_iterator end() const
        {
            return ones.end();
        }

        void clear() {
            ones.clear();
        }

    private:
        friend class boost::serialization::access;

        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
          ar & ones;
        }

        std::vector<U> ones;
    };

    template<typename U>
    std::ostream& operator<<(std::ostream &o, const SpVecGF2<U> &v) {
        if (v.size() > 0)
            o << "(" << v.size() << ") ";
        std::copy(v.begin(), v.end(), std::ostream_iterator<U> { o, " " });
        return o;
    }

} // parmcb

#endif
