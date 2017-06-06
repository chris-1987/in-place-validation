#ifndef MY_SUBSTR_H
#define MY_SUBSTR_H

//#define TEST_DEBUG

#include "io.h"

#if _MSC_VER
#pragma pack(push, 1)
#endif
template<typename alphabet_type, uint D>
struct Substr{

private: 

	typedef typename ExVector<alphabet_type>::vector alphabet_vector_type;

	typedef typename ExVector<uint8>::vector uint8_vector_type;

public:

	std::vector<alphabet_type> data; // store at most D characters

	bool is_short; // short or long
	
public:

	/// \brief retrieve a non-size-one short substr
	///
	/// \return false is no more; otherwise, true
	/// \note short substrs in a seq are separated by 0
	bool deserializeShort(alphabet_vector_type * _substr_ch_seq, typename alphabet_vector_type::const_iterator & _substr_ch_it, 
				typename uint8_vector_type::const_iterator & _substr_len_it) {

		if (_substr_ch_it == _substr_ch_seq->end()) {

			return false;
		}

		uint8 len = static_cast<uint8>(*_substr_len_it); // len < UINT8_MAX 

		++_substr_len_it;

		data.resize(len);

		for (uint8 i = 0; i < len; ++i) { data[i] = *_substr_ch_it; ++_substr_ch_it; }

		is_short = true;

		return true;
	}

	/// \brief retrieve the first D characters from a long substr
	///
	void deserializeLong(typename alphabet_vector_type::const_iterator & _substr_ch_it) {

		data.resize(D);

		for (uint8 i = 0; i < D; ++i) { data[i] = *_substr_ch_it; ++_substr_ch_it; }	

		is_short = false;
	} 
	
	
	/// \brief compare two substrs (at most one is long)
	/// 
	int cmp(const Substr &_b) const {

		auto iter_a = data.begin(), iter_b = _b.data.begin();

		for (; iter_a != data.end() && iter_b != _b.data.end(); ++iter_a, ++iter_b) {

			if (*iter_a < *iter_b) return -1;

			if (*iter_a > *iter_b) return +1;
		}

		if (iter_a != data.end()) return -1; // b is a prefix of a, a < b
		
		if (iter_b != _b.data.end()) return +1; // a is a prefix of b, a > b

		if (is_short == false) return -1; // a is long, then b is a prefix of a, a < b

		if (_b.is_short == false) return +1; // b is long, then a is a prefix of b, a > b


		return 0; // equal
	}

	void swap(Substr<alphabet_type, D> & _b) {

		std::swap(data, _b.data);
	
		std::swap(is_short, _b.is_short);
	} 

	void output() const{

		auto it = data.begin();

		for (; it != data.end(); ++it) {

			std::cerr << (uint32)*it << " ";
		}

		std::cerr << std::endl;
	}
};



#endif // MY_SUBSTR_H
