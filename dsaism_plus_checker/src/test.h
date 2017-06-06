#ifndef MY_TEST_H
#define MY_TEST_H

#include "common.h"

template<typename alphabet_vector_type>
void test_output_vector(alphabet_vector_type *_vector) {

	auto it = _vector->begin();

	for (; it != _vector->end(); ++it) {

		if (sizeof(*it) == sizeof(uint8)) {
			
			std::cerr << (uint32)*it << " ";
		}
		else {

			std::cerr << *it << " ";
		}
	}

	std::cerr << std::endl;
}

template<typename BlockInfo>
void test_output_blocks_info(std::vector<BlockInfo> & _blocks_info) {

	for (uint8 i = 0; i < _blocks_info.size(); ++i) {

		const BlockInfo & bi = _blocks_info[i];

		std::cerr << "id: " << (uint32)bi.m_id;
	
		std::cerr << " bpos: " << bi.m_beg_pos;
	
		std::cerr << " epos: " << bi.m_end_pos;

		std::cerr << " lms_num: " << bi.m_lms_num;

		std::cerr << " size: " << bi.m_size;

		std::cerr << " cap: " << bi.m_capacity << std::endl;

	}
}

template<typename T>
void test_output_array_ram(T * _arr, uint32 _size) {

	for (uint32 i = 0; i < _size; ++i) {

		if (sizeof (T) == sizeof(uint8)) {

			std::cerr << (uint32)_arr[i] <<" ";
		}
		else {
			
			std::cerr << _arr[i] << " ";
		}
	}

	std::cerr << std::endl;
}

template<typename T>
void test_output_element(T & _elem) {

	std::cerr << _elem << " ";
}





#endif // MY_TEST_H
