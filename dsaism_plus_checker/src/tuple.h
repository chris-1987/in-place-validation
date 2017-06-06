#ifndef MY_TUPLE_H
#define MY_TUPLE_H

/// \brief structure pair
///
template<typename T1, typename T2>
struct Pair{

	T1 first;

	T2 second;

	/// \brief constructor, default
	Pair(): first(0), second(0) {}

	/// \brief constructor, copy
	Pair(const Pair& _item) : first(_item.first), second(_item.second) {}

	/// \brief constructor, component
	Pair(const T1& _first, const T2& _second) : first(_first), second(_second) {}

	/// \brief min value
	static Pair& min_value() {

		static Pair min_val = Pair();

		return min_val;
	}

	/// \brief max vlaue
	static Pair& max_value() {

		static Pair max_val = Pair(std::numeric_limits<T1>::max(), std::numeric_limits<T2>::max());

		return max_val;
	}

	/// \brief output
	friend std::ostream& operator << (std::ostream _os, const Pair& _item) {

		return _os << "first: " << _item.first << " second: " << _item.second << std::endl;

	}
}__attribute__((packed));


/// \brief structure triple
///
template<typename T1, typename T2, typename T3>
struct Triple{

	T1 first;

	T2 second;

	T3 third;

	/// \brief constructor, default
	Triple(): first(0), second(0), third(0) {}

	/// \brief constructor, copy
	Triple(const Triple& _item) : first(_item.first), second(_item.second), third(_item.third) {}

	/// \brief constructor, component
	Triple(const T1& _first, const T2& _second, const T3& _third) : first(_first), second(_second), third(_third) {}

	/// \brief min value
	static Triple& min_value() {

		static Triple min_val = Triple();

		return min_val;
	}

	/// \brief max value
	static Triple& max_value() {

		static Triple max_val = Triple(std::numeric_limits<T1>::max(), std::numeric_limits<T2>::max(), std::numeric_limits<T3>::max());

		return max_val;
	}

	/// \brief output
	friend std::ostream& operator << (std::ostream _os, const Triple& _item) {

		return _os << "first: " << _item.first << " second: " << _item.second << " third: " << _item.third << std::endl;
	}

}__attribute__((packed));


/// \brief structure quadruple
///
template<typename T1, typename T2, typename T3, typename T4>
struct quadruple{

	T1 first;

	T2 second;

	T3 third;

	T4 forth;

	/// \brief constructor, default
	quadruple(): first(0), second(0), third(0), forth(0) {}

	/// \brief constructor, copy
	quadruple(const quadruple& _item) : first(_item.first), second(_item.second), third(_item.third), forth(_item.forth) {}

	/// \brief constructor, component
	quadruple(const T1& _first, const T2& _second, const T3& _third, const T4& _forth) : first(_first), second(_second), third(_third), forth(_forth) {}

	/// \brief min value
	static quadruple& min_value() {

		static quadruple min_val = quadruple();

		return min_val;
	}

	/// \brief max value
	static quadruple& max_value() {

		static quadruple max_val = quadruple(std::numeric_limits<T1>::max(), std::numeric_limits<T2>::max(), std::numeric_limits<T3>::max(), std::numeric_limits<T4>::max());

		return max_val;
	}

	/// \brief output
	friend std::ostream& operator << (std::ostream _os, const quadruple& _item) {

		return _os << "first: " << _item.first << " second: " << _item.second << " third: " << _item.third << " forth: " << _item.forth << std::endl;
	}

}__attribute__((packed));



#endif // MY_TUPLE_H

