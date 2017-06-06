#ifndef BYTESTREAM_H
#define BYTESTREAM_H


#include "namespace.h"
#include "mycommon.h"

WTL_BEG_NAMESPACE

class ByteStream {
public:
	//alias
	//type of allocation strategy

	//type of block identifier
	typedef stxxl::BID<BLOCK_SIZE> bid_type;

	//type of block identifer vector
	typedef std::vector<bid_type> bid_vector_type;

	//type of iterator for block identifier vector
	typedef typename bid_vector_type::iterator bid_iterator_type;

	//type of uint8 block
	typedef stxxl::typed_block<BLOCK_SIZE, uint8> block_type;

	//type of uint8 buf_ostream
	typedef stxxl::buf_ostream<block_type, bid_iterator_type> buf_ostream_type;

	//type of uint8 buf_istream
	typedef stxxl::buf_istream<block_type, bid_iterator_type> buf_istream_type;

	//type of vector for uint8 buf_ostream
	typedef std::vector<buf_ostream_type*> buf_ostream_vector_type;

	//type of vector for uint8 buf_istream 
	typedef std::vector<buf_istream_type*> buf_istream_vector_type;

	//data members
	//put data into the buf_ostream
	template<typename dataT> static buf_ostream_type& putData(buf_ostream_type & _os, const dataT & _data) {
		for (size_t i = 0; i < sizeof(dataT) / sizeof(typename buf_ostream_type::block_type::type); ++i) {
			_os << (((typename buf_ostream_type::block_type::type*)&_data)[i]);
		}
		return _os;
	}

	//get data from the buf_istream
	template<typename dataT> static buf_istream_type& getData(buf_istream_type & _is, const dataT & _data) {
		for (size_t i = 0; i < sizeof(dataT) / sizeof(typename buf_istream_type::block_type::type); ++i) {
			_is >> (((typename buf_istream_type::block_type::type*)&_data)[i]);
		}
		return _is;
	}
};


WTL_END_NAMESPACE


#endif
