#ifndef MY_SAIS_H
#define MY_SAIS_H

#include "common.h"
#include "utility.h"
#include "io.h"


//#define TEST_DEBUG2

class SAComputation;

/// \brief portal to SA computation on RAM
///
template<typename offset_type>
class SAIS{

private:

	typedef typename ExVector<offset_type>::vector offset_vector_type;

public:

	SAIS(offset_vector_type * _s, offset_vector_type *&_sa_reverse); 

};

/// \brief compute SA for the input string residing on EM
///
/// \note alphabet_type = offset_type > uint32
template<typename offset_type>
SAIS<offset_type>::SAIS(offset_vector_type *_s, offset_vector_type*& _sa_reverse) {

	// load _s into ram
	// note that characters in _s are named in a condensed way
	uint32 s_size = _s->size();

	uint32 *s = new uint32[s_size];

	typename offset_vector_type::bufreader_type *s_reader = new typename offset_vector_type::bufreader_type(*_s);

	uint32 max_alpha = 0;

	for (uint32 i = 0; i < s_size; ++i, ++(*s_reader)) {

		s[i] = static_cast<uint32>(*(*s_reader)); // correct

#ifdef TEST_DEBUG2
		std::cerr << "_s[i]: " << *(*s_reader) << " s[i]: " << s[i] << std::endl;
#endif
		max_alpha = std::max(max_alpha, s[i]);
	}	

	delete s_reader; s_reader = nullptr;

	// compute sa
	uint32 *sa = new uint32[s_size];

	SAComputation(s, s_size, max_alpha, sa);

	// output sa reversely
	_sa_reverse = new offset_vector_type();

#ifdef TEST_DEBUG2

	std::cerr << "sa: \n";
#endif

	for (uint32 i = s_size - 1; ; --i) {

		_sa_reverse->push_back(static_cast<offset_type>(sa[i]));

#ifdef TEST_DEBUG2

		std::cerr << sa[i] << " ";
#endif

		if (i == 0) break;
	}	

#ifdef TEST_DEBUG2
	std::cerr << std::endl;
#endif

	delete [] s;

	delete [] sa;
}



/// \brief compute SA on RAM
///
class SAComputation{

	static const uint32 uint32_MAX = std::numeric_limits<uint32>::max(); 

public:

	SAComputation(uint32 * _s, const uint32 _sLen, const uint32 _alpha, uint32 * _sa);

	void getBuckets(uint32 * _s, const uint32 _sLen, uint32 * _bkt, const uint32 _bktNum, const bool _end);

	void induceL(uint32 * _s, BitWrapper & _t, uint32 * _sa, const uint32 _sLen, uint32 * _bkt, const uint32 _alpha);

	void induceS(uint32 * _s, BitWrapper & _t, uint32 * _sa, const uint32 _sLen, uint32 * _bkt, const uint32 _alpha);  
};

/// \brief ctor
///
SAComputation::SAComputation(uint32 * _s, const uint32 _sLen, const uint32 _alpha, uint32 * _sa) {
	uint32 i, j;

	char * t_buf = new char[_sLen / 8 + 1]; BitWrapper t(t_buf);

	//compute t
	t.set(_sLen - 1, S_TYPE); t.set(_sLen - 2, L_TYPE); 
	for (i = _sLen - 3; ; --i) {
		t.set(i, (_s[i] < _s[i + 1] || (_s[i] == _s[i + 1] && t.get(i + 1) == S_TYPE)) ? S_TYPE : L_TYPE);
		if (i == 0) break;
	}

	// sort all the S-substrings
	uint32 *bkt = new uint32[_alpha + 1]; 
	getBuckets(_s, _sLen, bkt, _alpha + 1, true);


	for (i = 0; i < _sLen; ++i) _sa[i] = uint32_MAX;//init sa
	for (i = _sLen - 3; i >= 1; --i) {	// find lms-chars in _s[1, _sLen - 3]
		if (t.get(i) && !t.get(i - 1)) {
			_sa[bkt[_s[i]]--] = i;
		}
	}
	_sa[0] = _sLen - 1;	//_s[_sLen - 1] is the sentinel smallest than any other characters

	induceL(_s, t, _sa, _sLen, bkt, _alpha); //induce l-type lms-prefix
	induceS(_s, t, _sa, _sLen, bkt, _alpha); //induce s-type lms-prefix

	delete[] bkt;

	// compact all the sorted substrings into the first n1 items of s
	uint32 s1_len = 0;
	for (i = 0; i < _sLen; ++i) {
		if (_sa[i] > 0 && t.get(_sa[i]) && !t.get(_sa[i] - 1)) {
			_sa[s1_len++] = _sa[i];
		}
	}

	for (i = s1_len; i < _sLen; ++i) _sa[i] = uint32_MAX; //init

	// find the lexicographic names of all substrings
	uint32 name = 0;
	uint32 prev = uint32_MAX;
	for (i = 0; i< s1_len; ++i) {
		uint32 pos = _sa[i]; bool diff = false;
		for (uint32 d = 0; d < _sLen; ++d) {
			if (prev == uint32_MAX || pos + d == _sLen - 1 || prev + d == _sLen - 1 || _s[pos + d] != _s[prev + d] || t.get(pos + d) != t.get(prev + d)) {
				diff = true;
				break;
			}
			else {
				if (d > 0 && ((t.get(pos + d) && !t.get(pos + d - 1)) || (t.get(prev + d) && !t.get(prev + d - 1)))) {
					break;
				}
			}
		}
		if (diff) {
			++name;
			prev = pos;
		}
		pos = pos / 2;
		_sa[s1_len + pos] = name - 1;
	}
	for (i = _sLen - 1, j = _sLen - 1; i >= s1_len; --i)
		if (_sa[i] != uint32_MAX) _sa[j--] = _sa[i];

	// s1 is done now
	uint32 *sa1 = _sa;
	uint32 *s1 = _sa + _sLen - s1_len;
	// stage 2: solve the reduced problem

	// recurse if names are not yet unique
	if (name < s1_len) {
		SAComputation(s1, s1_len, name - 1, sa1);
	}
	else { // generate the suffix array of s1 directly
		for (i = 0; i < s1_len; ++i) {
			sa1[s1[i]] = i;
		}
	}

	// stage 3: induce the result for the original problem
	// put all left-most S characters into their buckets
	bkt = new uint32[_alpha + 1];
	getBuckets(_s, _sLen, bkt, _alpha + 1, true); // find ends of buckets

	j = 0;
	for (i = 1; i < _sLen; ++i) {
		if (t.get(i) && !t.get(i - 1)) s1[j++] = i;// get p1
	}
	for (i = 0; i < s1_len; ++i) sa1[i] = s1[sa1[i]]; // get index in s1
	for (i = s1_len; i < _sLen; ++i) _sa[i] = uint32_MAX; // init SA[n1..n-1]
	for (i = s1_len - 1; ; --i) {
		j = _sa[i]; _sa[i] = uint32_MAX;
		_sa[bkt[_s[j]]--] = j;
		if (i == 0) break;
	}
	induceL(_s, t, _sa, _sLen, bkt, _alpha);
	induceS(_s, t, _sa, _sLen, bkt, _alpha);

	delete[] bkt; bkt = nullptr;
	delete[] t_buf; t_buf = nullptr;
}


//@usage: compute bucket size.
void SAComputation::getBuckets(uint32 * _s, const uint32 _sLen, uint32 * _bkt, const uint32 _bktNum, const bool _end) {
	uint32 i;
	uint32 sum = 0;
	for (i = 0; i < _bktNum; ++i) _bkt[i] = 0; // clear all buckets
	for (i = 0; i < _sLen; ++i) ++_bkt[_s[i]]; // compute the size of each bucket
	for (i = 0; i < _bktNum; ++i) { sum += _bkt[i]; _bkt[i] = _end ? sum - 1 : sum - _bkt[i]; }
}

//@usage: induce L.	
void SAComputation::induceL(uint32 * _s, BitWrapper & _t, uint32 * _sa, const uint32 _sLen, uint32 * _bkt, const uint32 _alpha) {
	uint32 i, j;
	getBuckets(_s, _sLen, _bkt, _alpha + 1, false); // find heads of buckets
	for (i = 0; i < _sLen; ++i) {
		if (_sa[i] != uint32_MAX && _sa[i] != 0) {
			j = _sa[i] - 1;
			if (!_t.get(j)) { // _sa[i] is unsigned, if _sa[i] == 0, then j does not exist. 
				_sa[_bkt[_s[j]]++] = j;
			}
		}
	}
}

//@usage: induce S.
void SAComputation::induceS(uint32 * _s, BitWrapper & _t, uint32 * _sa, const uint32 _sLen, uint32 * _bkt, const uint32 _alpha) {
	uint32 i, j;
	getBuckets(_s, _sLen, _bkt, _alpha + 1, true); // find ends of buckets
	for (i = _sLen - 1; ; --i) {
		if (_sa[i] != uint32_MAX && _sa[i] != 0) {
			j = _sa[i] - 1;
			if (_t.get(j)) { // _sa[i] is unsigned, if _sa[i] == 0, then j == EMPTY, which does not exist. 
				_sa[_bkt[_s[j]]--] = j;
			}
		}
		if (i == 0) break;
	}
}

#endif // MY_SAIS_H
