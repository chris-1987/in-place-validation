#ifndef MY_CHECKER_H
#define MY_CHECHER_H

#include "common.h"
#include "tuple.h"
#include "tuple_sorter.h"
#include "io.h"

/// \brief perform checking after construction
///
template<typename alphabet_type, typename offset_type>
class Checker{

private:

	typedef typename ExVector<alphabet_type>::vector alphabet_vector_type;

	typedef typename ExVector<offset_type>::vector offset_vector_type; 

private:

	std::string m_s_fname;

	std::string m_sa_fname;

public:

	/// \brief ctor
	///
	Checker(const std::string & _s_fname, const std::string & _sa_fname) :
		m_s_fname(_s_fname), m_sa_fname(_sa_fname) {}


	/// \brief check
	///
	bool run() {

		// sort <sa[i], i> by sa[i]
		typedef Pair<offset_type, offset_type> pair_type;

		typedef TupleAscCmp1<pair_type> pair_comparator_type;

		typedef typename ExSorter<pair_type, pair_comparator_type>::sorter sorter_type;

		sorter_type *sorter = new sorter_type(pair_comparator_type(), MAX_MEM / 2);

		stxxl::syscall_file *sa_file = new stxxl::syscall_file(m_sa_fname, stxxl::syscall_file::RDWR | stxxl::syscall_file::DIRECT);

		offset_vector_type *sa = new offset_vector_type(sa_file);

		typename offset_vector_type::const_iterator it = sa->begin();

		for (offset_type i = 1; it != sa->end(); ++it, ++i) {

			sorter->push(pair_type(*it, i)); // start ranking from 1 (0 is reserved for the sentinel)
		}

		delete sa; sa = nullptr;

		delete sa_file; sa_file = nullptr;

		sorter->sort();

		// sort <isa[i], s[i], isa[i + 1]> by isa[i]
		typedef Triple<offset_type, alphabet_type, offset_type> triple_type;

		typedef TupleAscCmp1<triple_type> triple_comparator_type;

		typedef typename ExSorter<triple_type, triple_comparator_type>::sorter sorter_type2;

		sorter_type2 *sorter2 = new sorter_type2(triple_comparator_type(), MAX_MEM / 2);

		stxxl::syscall_file *s_file = new stxxl::syscall_file(m_s_fname, stxxl::syscall_file::RDWR | stxxl::syscall_file::DIRECT);

		alphabet_vector_type *s = new alphabet_vector_type(s_file);

		typename alphabet_vector_type::const_iterator it2 = s->begin();

		alphabet_type cur_ch;

		offset_type cur_rank, next_rank;

		offset_type cur_pos, next_pos;

		cur_rank = (*sorter)->second, cur_pos = (*sorter)->first, ++(*sorter);

		for (offset_type i = 0; !sorter->empty(); ++it2, ++(*sorter), ++i) {

			if (cur_pos != i) {
			
				std::cerr << "not a permutation\n";

				return false;
			}

			cur_ch = *it2, next_rank = (*sorter)->second, next_pos = (*sorter)->first;

			sorter2->push(triple_type(cur_rank, cur_ch, next_rank));

			cur_rank = next_rank, cur_pos = next_pos;
		}

		cur_ch = *it2, next_rank = offset_type(0); // the rank of the sentinel is set to 0

		sorter2->push(triple_type(cur_rank, cur_ch, next_rank));

		delete sorter; sorter = nullptr;

		delete s; s = nullptr;

		delete s_file; s_file = nullptr;
		
		sorter2->sort();

		// scan to check the result
		{

			stxxl::syscall_file *sa_file = new stxxl::syscall_file(m_sa_fname, stxxl::syscall_file::RDWR | stxxl::syscall_file::DIRECT);

			offset_vector_type *sa = new offset_vector_type(sa_file);

			typename offset_vector_type::const_iterator it = sa->begin();
				
			alphabet_type last_ch, cur_ch;

			offset_type last_rank, cur_rank;

			offset_type last_suc_rank, cur_suc_rank;

			offset_type last_pos, cur_pos;

			last_rank = (*sorter2)->first, last_ch = (*sorter2)->second, last_suc_rank = (*sorter2)->third, ++(*sorter2);

			last_pos = *it, ++it;

			for (; !sorter2->empty(); ++(*sorter2), ++it) {

				cur_rank = (*sorter2)->first, cur_ch = (*sorter2)->second, cur_suc_rank = (*sorter2)->third;
			
				cur_pos = *it;

				if (last_ch > cur_ch || (last_ch == cur_ch && last_suc_rank > cur_suc_rank)) {

					std::cerr << "last_rank: " << last_rank << " ";

					std::cerr << "last_ch: " << (uint32)last_ch << " ";

					std::cerr << "last_suc_rank: " << last_suc_rank << " ";

					std::cerr << "last_pos: " << last_pos << std::endl;

					std::cerr << "cur_rank: " << cur_rank << " ";
				
					std::cerr << "cur_ch: " << (uint32)cur_ch << " ";

					std::cerr << "cur_suc_rank: " << cur_suc_rank << " ";

					std::cerr << "cur_pos: " << cur_pos << std::endl;

					return false;	
				}

				last_rank = cur_rank, last_ch = cur_ch, last_suc_rank = cur_suc_rank, last_pos = cur_pos;
			}	

			delete sorter2; sorter2 = nullptr;

			delete sa; sa = nullptr;

			delete sa_file; sa_file = nullptr;

			return true;
		}

	}
};



#endif // my checker
