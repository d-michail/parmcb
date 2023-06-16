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

	SpVecFP() :
			p(3) {
	}

	SpVecFP(const P &p) :
			p(p) {
	}

	SpVecFP(const P &p, const std::size_t index) :
			p(p) {
#ifdef PARMCB_INVARIANTS_CHECK
            assert(index >= 0);
#endif
		entries.push_back(boost::make_tuple(index, 1));
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

	SpVecFP<P>& operator=(const std::size_t &index) {
		entries.clear();
#ifdef PARMCB_INVARIANTS_CHECK
            assert(index >= 0);
#endif
		entries.push_back(boost::make_tuple(index, 1));
		return *this;
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
			std::size_t index = boost::get<0>(entry);

			entry_type v_entry = *v_it;
			std::size_t v_index = boost::get<0>(v_entry);

			if (index > v_index) {
				v_it++;
			} else if (index < v_index) {
				it++;
			} else {
				P v = (boost::get<1>(entry) * boost::get<1>(v_entry)) % p;
				res = (res + v) % p;
				it++;
				v_it++;
			}
		}
		return normalize(res);
	}

	SpVecFP<P> operator+(const SpVecFP<P> &v) const {
		SpVecFP<P> res(p);
		auto it = entries.begin(), it_e = entries.end();
		auto v_it = v.entries.begin(), v_it_e = v.entries.end();

		// now add them
		while (it != it_e && v_it != v_it_e) {
			entry_type entry = *it;
			std::size_t index = boost::get<0>(entry);
			P value = boost::get<1>(entry);

			entry_type v_entry = *v_it;
			std::size_t v_index = boost::get<0>(v_entry);
			P v_value = boost::get<1>(v_entry);

			if (index > v_index) {
				res.entries.push_back(boost::make_tuple(v_index, v_value));
				v_it++;
			} else if (index < v_index) {
				res.entries.push_back(boost::make_tuple(index, value));
				it++;
			} else {
				P v = normalize((value + v_value) % p);
				if (v != 0) {
					res.entries.push_back(boost::make_tuple(index, v));
				}
				it++;
				v_it++;
			}
		}

		// append remaining stuff
		while (it != it_e) {
			entry_type entry = *it;
			std::size_t index = boost::get<0>(entry);
			P value = boost::get<1>(entry);
			res.entries.push_back(boost::make_tuple(index, value));
			it++;
		}
		while (v_it != v_it_e) {
			entry_type v_entry = *v_it;
			std::size_t v_index = boost::get<0>(v_entry);
			P v_value = boost::get<1>(v_entry);
			res.entries.push_back(boost::make_tuple(v_index, v_value));
			v_it++;
		}

		return res;
	}

	SpVecFP<P>& operator+=(const SpVecFP<P> &v) {
		*this = *this + v;
		return *this;
	}

	SpVecFP<P> operator-() const {
	    SpVecFP<P> res(p);
	    auto it = entries.begin(), it_e = entries.end();

		while (it != it_e) {
			entry_type entry = *it;
			std::size_t index = boost::get<0>(entry);
			P v = normalize(-boost::get<1>(entry));
			if (v != 0) {
				res.entries.push_back(boost::make_tuple(index, v));
			}
			it++;
		}
	    return res;
	}

	SpVecFP<P>& operator-=(const SpVecFP<P> &v) {
		*this = *this + (-v);
		return *this;
	}

	SpVecFP<P> operator*(const P &a) const {
		SpVecFP<P> res(p);
		auto it = entries.begin(), it_e = entries.end();
		while (it != it_e) {
			entry_type entry = *it;
			std::size_t index = boost::get<0>(entry);
			P value = boost::get<1>(entry);
			P v = normalize((value * a) % p);
			if (v != 0) {
				res.entries.push_back(boost::make_tuple(index, v));
			}
			it++;
		}
		return res;
	}

	SpVecFP<P>& operator*=(const P &a) {
		*this = *this * a;
		return *this;
	}

	size_type size() const {
		return entries.size();
	}

	const_iterator begin() const {
		return entries.begin();
	}

	const_iterator end() const {
		return entries.end();
	}

	P prime() const {
		return p;
	}

	void clear() {
		entries.clear();
	}

private:
	friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive &ar, const unsigned int version) {
		ar & p;
		ar & entries;
	}

	P normalize(P v) const {
		while (v < 0)
			v += p;   // make [-i]_p = [p-i]_p
		while (v >= p)
			v -= p; // make [i+p]_p = [i]_p
		return v;
	}

	std::vector<entry_type> entries;
	P p;
};

template<typename P>
std::ostream& operator<<(std::ostream &o, const SpVecFP<P> &v) {
	auto v_end = v.end();
	for (auto it = v.begin(); it != v_end; it++) {
		auto tuple = *it;
		o << "(i=" << boost::get<0>(tuple) << ",v=" << boost::get<1>(tuple) << ")"
				<< " ";
	}
	o << "(mod " << v.prime() << ")";
	return o;
}

} // parmcb

#endif
