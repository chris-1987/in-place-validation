#ifndef MY_IO_H
#define MY_IO_H

#include "common.h"

/// \brief stxxl vector generator
///
template<typename element_type>
struct ExVector{

	typedef typename stxxl::VECTOR_GENERATOR<element_type, 8 / sizeof(element_type) + 1, 2, K_512 * sizeof(element_type)>::result vector;
};


/// \brief stxxl sorter generator
///
template<typename tuple_type, typename comparator_type>
struct ExSorter{

	typedef typename stxxl::sorter<tuple_type, comparator_type> sorter;
};


/// \brief stxxl priority queue
///
template<typename tuple_type, typename comparator_type, uint32 mem_size, uint64 max_items, uint64 tune = 6>
struct ExHeap{

	typedef typename stxxl::PRIORITY_QUEUE_GENERATOR<tuple_type, comparator_type, mem_size, max_items, tune>::result heap;
};

/// \brief stxxl queue
template<typename tuple_type>
struct ExQueue{

	typedef typename stxxl::queue<tuple_type> queue;
};

/// \brief stxxl deque
template<typename tuple_type>
struct ExDeque{

	typedef typename stxxl::deque<tuple_type> deque;
};

#endif
