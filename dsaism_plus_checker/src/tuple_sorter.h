#ifndef MY_TUPLE_SORTER_H
#define MY_TUPLE_SORTER_H

#include "common.h"
	
#include "stxxl/sorter"


/// \brief sort tuples by 1st component in ascending order
///
template<typename tuple_type>
struct TupleAscCmp1{

	bool operator()(const tuple_type &_a, const tuple_type &_b) const {

		return _a.first < _b.first;
	}

	tuple_type min_value() const {

		return tuple_type::min_value();
	}

	tuple_type max_value() const {

		return tuple_type::max_value();
	}
};

/// \brief sort tuples by 1st component in descending order
///
template<typename tuple_type>
struct TupleDscCmp1{

	bool operator()(const tuple_type &_a, const tuple_type &_b) const {

		return _a.first > _b.first;
	}

	tuple_type min_value() const {

		return tuple_type::max_value();
	}

	tuple_type max_value() const {

		return tuple_type::min_value();
	}
};

/// \brief sort tuples by 1st component in ascending order
///
template<typename tuple_type>
struct TupleAscCmp2{

	bool operator()(const tuple_type &_a, const tuple_type &_b) const {

		if (_a.first == _b.first) return _a.second < _b.second;

		return _a.first < _b.first;
	}

	tuple_type min_value() const {

		return tuple_type::min_value();
	}

	tuple_type max_value() const {

		return tuple_type::max_value();
	}
};

/// \brief sort tuples by 1st component in descending order
///
template<typename tuple_type>
struct TupleDscCmp2{

	bool operator()(const tuple_type &_a, const tuple_type &_b) const {

		if (_a.first == _b.first) return _a.second > _b.second;

		return _a.first > _b.first;
	}

	tuple_type min_value() const {

		return tuple_type::max_value();
	}

	tuple_type max_value() const {

		return tuple_type::min_value();
	}
};
/// \brief sort tuples by <1st, 2nd, 3rd> components in ascending order
///
template<typename tuple_type>
struct TupleAscCmp3{


	bool operator()(const tuple_type &_a, const tuple_type &_b) const {

		if (_a.first == _b.first) {

			if (_a.second == _b.second) return _a.third < _b.third;

			return _a.second < _b.second;
		}

		return _a.first < _b.first;
	}

	tuple_type & min_value() const {

		return tuple_type::min_value();
	}

	tuple_type & max_value() const {

		return tuple_type::max_value();
	}
};

/// \brief sort tuples by <1st, 2nd, 3rd> components in descending order
///
template<typename tuple_type>
struct TupleDscCmp3{

	bool operator() (const tuple_type &_a, const tuple_type &_b) const {

		if (_a.first == _b.first) {

			if (_a.second == _b.second) return _a.third > _b.third;

			return _a.second > _b.second;
		}

		return _a.first > _b.first;

	}

	tuple_type & min_value() const {

		return tuple_type::max_value();
	}

	tuple_type & max_value() const {

		return tuple_type::min_value();
	}
};


#endif // MY_TUPLE_SORTER_H
