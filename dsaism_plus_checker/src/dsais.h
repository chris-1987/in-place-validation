#ifndef MY_DSAIS_H
#define MY_DSAIS_H

#include "common.h"
#include "io.h"
#include "sais.h"
#include "tuple.h"
#include "tuple_sorter.h"
#include "substr_sorter.h"
#include "test.h"

#include <string>
#include <fstream>
#include <cassert>

//#define TEST_DEBUG1 // output
//#define TEST_DEBUG2 // assertion

template<typename alphabet_type, typename offset_type, uint8 D>
class DSAComputation;

/// \brief portal to DSAComputation
///
template<typename alphabet_type, typename offset_type>
class DSAIS{

private:

	typedef typename ExVector<alphabet_type>::vector alphabet_vector_type;

	typedef typename ExVector<offset_type>::vector offset_vector_type;

private:

	const std::string m_s_fname; ///< input string file name
	
	const std::string m_sa_fname; ///< output SA file name

public:

	/// \brief ctor 
	///
	DSAIS(const std::string & _s_fname, const std::string & _sa_fname) : m_s_fname(_s_fname), m_sa_fname(_sa_fname) {}

	/// \brief create DSAComputation
	///
	void run() {

		// append a sentinel to the input string
		stxxl::syscall_file *s_file = new stxxl::syscall_file(m_s_fname, stxxl::syscall_file::RDWR | stxxl::syscall_file::DIRECT);

		alphabet_vector_type *s_origin = new alphabet_vector_type(s_file);

		typename alphabet_vector_type::bufreader_type s_origin_reader(*s_origin);

		alphabet_vector_type *s_target = new alphabet_vector_type(); s_target->resize(s_origin->size() + 1); // append a sentinel

		typename alphabet_vector_type::bufwriter_type s_target_writer(*s_target);

		uint64 s_origin_len = s_origin->size();

		for (; !s_origin_reader.empty(); ++s_origin_reader) {

			s_target_writer << *s_origin_reader;
		}
	
		s_target_writer << alphabet_type(0);

		s_target_writer.finish();

		// clear
		delete s_origin; s_origin = nullptr;

		delete s_file; s_file = nullptr;

		// compute sa_reverse
		offset_vector_type *sa_reverse = nullptr;
	
		DSAComputation<alphabet_type, offset_type, D_LOW> dsac(s_target, 0, sa_reverse);
	
		dsac.run();
	
		// clear
		delete s_target; s_target = nullptr;
	
		// produce sa from sa_reverse
		typename offset_vector_type::bufreader_reverse_type sa_reader(*sa_reverse);
		
		stxxl::syscall_file *sa_file = new stxxl::syscall_file(m_sa_fname, 
						stxxl::syscall_file::CREAT | stxxl::syscall_file::RDWR | stxxl::syscall_file::DIRECT);
	
		offset_vector_type *sa = new offset_vector_type(sa_file); sa->resize(s_origin_len);
	
		typename offset_vector_type::bufwriter_type sa_writer(*sa);

		++sa_reader; // skip the sentinel
	
		for (; !sa_reader.empty(); ++sa_reader) {

			sa_writer << *sa_reader;
		}
		
		sa_writer.finish();
	
		// clear	
		delete sa; sa = nullptr;
	
		delete sa_file; sa_file = nullptr;
		
		delete sa_reverse; sa_reverse = nullptr;
	
		return;
	}	
};

/// \brief compute SA by DSA-IS
///
template<typename alphabet_type, typename offset_type, uint8 D>
class DSAComputation{

private:
	
	typedef typename ExVector<alphabet_type>::vector alphabet_vector_type;

	typedef typename ExVector<offset_type>::vector offset_vector_type;

	typedef typename ExVector<uint8>::vector uint8_vector_type;

	typedef typename ExVector<uint32>::vector uint32_vector_type;

private:
	
	const alphabet_type ALPHA_MAX;

	const alphabet_type ALPHA_MIN;

	const offset_type OFFSET_MAX;

	const offset_type OFFSET_MIN;

private:

	/// \brief block infomation collector
	///
	struct BlockInfo{

	public:

		uint32 m_capacity; ///< capacity for multi-block

		uint64 m_beg_pos; ///< starting position, global

		uint64 m_end_pos; //< ending position, global

		uint64 m_size; ///< number of characters in the block, non-multi-block may exceed the capacity

		uint32 m_lms_num; ///< number of S*-substrs in the block

		uint8 m_id; ///< id

	public:

		/// \brief ctor
		///
		BlockInfo(const uint64 & _end_pos, const uint8 _id) : m_end_pos(_end_pos), m_id(_id) {

			if (D == D_LOW) {

				double div = sizeof(alphabet_type) + sizeof(uint32) + sizeof(uint32) + sizeof(uint32) + (double)1 / 4;

				m_capacity = MAX_MEM2 / div;
			}
			else {

				double div = std::max((double)sizeof(alphabet_type) + sizeof(uint32) + sizeof(alphabet_type) + sizeof(uint32), sizeof(alphabet_type) + sizeof(uint32) + sizeof(uint32) + sizeof(uint32) + (double)1 / 4);

				m_capacity = MAX_MEM2 / div;
			}

			m_size = 0;

			m_lms_num = 0;
		}

		/// \brief assign copy
		///
		BlockInfo& operator=(const BlockInfo & _bi) {

			m_capacity = _bi.m_capacity;

			m_beg_pos = _bi.m_beg_pos;

			m_end_pos = _bi.m_end_pos;

			m_size = _bi.m_size;

			m_lms_num = _bi.m_lms_num;

			m_id = _bi.m_id;

			return *this;
		}


		/// \brief check if the block is empty
		///
		bool is_empty() const {

			return m_size == 0;
		}

		/// \brief check if the block contains more than one S*-substr
		///
		bool is_multi() const {

			return m_lms_num > 1;
		}

		/// \brief check if the block contains only single S*-substr
		///
		bool is_single() const {

			return m_lms_num == 1;
		}

		/// \brief check if the block can afford the S*-substr of a certain size
		///
		/// \note if the block is empty, it can afford an S*-substr of any size
		bool try_fill(const uint64 _lms_size) {

			if (is_empty() || m_size + _lms_size - 1 <= m_capacity) {

				return true;
			}

			return false;
		}

		/// \brief fill the S*-substr into the block
		///
		/// \note If empty, insert all the characters; otherwise, insert all but the last character
		void fill(const uint64 _lms_size) {

			if (is_empty()) {

				m_size = _lms_size;
			}
			else {

				m_size = m_size + _lms_size - 1;
			}

			++m_lms_num;

			return;
		}
		
		/// \brief close the block
		///
		void close() {

			m_beg_pos = m_end_pos - m_size + 1;

			return;
		}
	};

private:

	std::vector<BlockInfo> m_blocks_info;

	std::vector<uint8> m_block_id_of_samplings; ///< block_id of samplings
	
	alphabet_vector_type *& m_s;

	const uint64 m_s_len; 

	const uint32 m_level; ///< recursion level

	offset_vector_type *& m_sa_reverse; ///< output SA
	
	std::vector<alphabet_vector_type*> m_short_ch_seqs; ///< characters for short S*-substrs in sorted order (blockwise)

	std::vector<uint8_vector_type*> m_short_len_seqs; ///< lengths for short S*-substrs in sorted order (blockwise)

	std::vector<uint32_vector_type*> m_short_pos_seqs; ///< starting positions for short S*-substrs in sorted order (blockwise) 

	std::vector<alphabet_vector_type*> m_long_ch_seqs; ///< heading D characters for long S*-substrs in sorted order (blockwise)

	std::vector<alphabet_vector_type*> m_long_l_bwt_seqs; ///< prececding characters required for inducing L-type long substrs (blockwise)

	std::vector<alphabet_vector_type*> m_long_s_bwt_seqs; ///< prececding characters required for inducing S-type long substrs (blockwise)

	offset_vector_type *m_long_pos_seq; ///< starting positions for long S*-substrs in sorted order (global)

	uint8_vector_type *m_long_aux_seq; ///< aux for long S*-substrs in sorted order, 7 bits for block_id, 1 bit for diff (global)

	offset_vector_type *m_s1; ///< reduced string

	std::vector<offset_vector_type*> m_lms_rank_seqs; ///< S*-suffixes' ranks (blockwise)

	std::vector<alphabet_vector_type*> m_l_bwt_seqs; ///< preceding characters required for inducing L-type suffixes (blockwise)

	std::vector<alphabet_vector_type*> m_s_bwt_seqs; ///< preceding characters required for inducing S-type suffixes (blockwise)

#ifdef TEST_DEBUG2
	std::vector<offset_vector_type*> m_l_bpos_seqs;

	std::vector<offset_vector_type*> m_s_bpos_seqs;

#endif

public:

	DSAComputation(alphabet_vector_type *& _s, const uint32 _level, offset_vector_type *& _sa_reverse);

	void run();

	bool sortSStarGlobal();

	uint64 partitionS();

	void compute_block_id_of_samplings();

	uint8 getBlockId(const offset_type & _pos);

	void sortSStarBlock(const BlockInfo & _block_info);

	void sortSStarSingleBlock(const BlockInfo & _block_info);

	template<bool FORMAT>
	void sortSStarMultiBlock(const BlockInfo & _block_info);

	template<typename U>
	void getBuckets(const U * _s, const uint32 _s_size, uint32 * _bkt, const uint32 _bkt_num, const bool _end);

	void formatMultiBlock(const BlockInfo & _block_info, alphabet_type *_block, uint32 *_fblock, uint32 & _max_alpha);

	void sortLongSStarGlobal();

	bool mergeSortedSStarGlobal();
	
	void sortSuffixGlobal(offset_vector_type *&_sa1_reverse);

	void sortSuffixSingleBlock();

	void sortSuffixLeftmostBlock(const BlockInfo & _block_info);

	template<bool FORMAT>
	void sortSuffixMultiBlock(const BlockInfo & _block_info);

	void mergeSortedSuffixGlobal();
};


/// \brief ctor 
///
template<typename alphabet_type, typename offset_type, uint8 D>
DSAComputation<alphabet_type, offset_type, D>::DSAComputation(alphabet_vector_type *& _s, const uint32 _level, offset_vector_type *& _sa_reverse) : ALPHA_MAX(std::numeric_limits<alphabet_type>::max()), ALPHA_MIN(std::numeric_limits<alphabet_type>::min()),
	OFFSET_MAX(std::numeric_limits<offset_type>::max()), OFFSET_MIN(std::numeric_limits<offset_type>::min()),
	m_s(_s), m_s_len(m_s->size()), m_level(_level), m_sa_reverse(_sa_reverse){}


/// \brief run
///
template<typename alphabet_type, typename offset_type, uint8 D>
void DSAComputation<alphabet_type, offset_type, D>::run() {


#ifdef TEST_DEBUG1

	test_output_vector<alphabet_vector_type>(m_s);
#endif

	// sort S*-substrs
	bool is_unique = sortSStarGlobal();

	
	offset_vector_type *sa1_reverse = nullptr;

	// check recursion condition
	if (is_unique == false) {

		if (MAX_MEM >= m_s1->size() * (sizeof(uint32) + sizeof(uint32) + sizeof(uint32) + (double)1 / 8)) {

#ifdef TEST_DEBUG1
			std::cerr << "recurse in RAM\n";
#endif

			SAIS<offset_type>(m_s1, sa1_reverse);
		}
		else {

#ifdef TEST_DEBUG1
			std::cerr << "recurse in EM\n";
#endif

			DSAComputation<offset_type, offset_type, D_HIGH> dsac_recurse(m_s1, m_level + 1, sa1_reverse);

			dsac_recurse.run();
		}
	}

	// sort suffixes
	sortSuffixGlobal(sa1_reverse);

	return;
}

/// \brief sort S*-substr
///
template<typename alphabet_type, typename offset_type, uint8 D>
bool DSAComputation<alphabet_type, offset_type, D>::sortSStarGlobal() {

	uint64 lms_num = partitionS(); // partition s into blocks

#ifdef TEST_DEBUG2

	std::cerr << "lms_num: " << lms_num << std::endl;

	test_output_blocks_info(m_blocks_info);

#endif

	if (lms_num == 0) { // no S*-substrs except for the sentinel, the reduced string only contain the sentinel, return true

		m_s1 = new offset_vector_type();

		m_s1->push_back(OFFSET_MIN);		

		return true;
	}
	
	for (uint8 i = 0; i < m_blocks_info.size(); ++i) { // sort S*-substrs in each block

		sortSStarBlock(m_blocks_info[i]);
	}

	sortLongSStarGlobal(); // sort all the long S*-substrs using external memory

	bool is_unique = mergeSortedSStarGlobal(); // merge long and short S*-substrs

	return is_unique;
}

/// \brief partition s into blocks
///
/// \return number of S*-substrs 
template<typename alphabet_type, typename offset_type, uint8 D>
uint64 DSAComputation<alphabet_type, offset_type, D>::partitionS() {

	// scan s leftward to find all the S*-substrs
	alphabet_type cur_ch, last_ch;

	uint8 cur_t, last_t;

	uint64 lms_end_pos, lms_size;

	typename alphabet_vector_type::bufreader_reverse_type s_rev_reader(*m_s);
	
	uint64 lms_num = 0;

#ifdef TEST_DEBUG2

	l_num = 0, s_num = 0;
#endif


	lms_end_pos = m_s->size() - 1;

	BlockInfo block_info = BlockInfo(lms_end_pos, m_blocks_info.size());

	lms_size = 1, ++s_rev_reader; // rightmost is the sentinel

#ifdef TEST_DEBUG2
	
	s_num++;
#endif

	cur_ch = *s_rev_reader, cur_t = L_TYPE, ++lms_size, ++s_rev_reader; // next on the left is L-type

#ifdef TEST_DEBUG2

	l_num++;
#endif 

	last_ch = cur_ch, last_t = cur_t;

	for (; !s_rev_reader.empty(); ++s_rev_reader) {

		cur_ch = *s_rev_reader;

		cur_t = (cur_ch < last_ch || (cur_ch == last_ch && last_t == S_TYPE)) ? S_TYPE : L_TYPE;

#ifdef TEST_DEBUG2
		if (cur_t == S_TYPE) s_num++;
		else l_num++;
#endif

		if (cur_t == L_TYPE && last_t == S_TYPE) { // last_ch is S*-type, find an S*-substr

			++lms_num;

			if (block_info.try_fill(lms_size) == false) {

				block_info.close();

				m_blocks_info.push_back(block_info);

				block_info = BlockInfo(lms_end_pos, m_blocks_info.size());
			}

			block_info.fill(lms_size);

			lms_end_pos = lms_end_pos - lms_size + 1;

			lms_size = 1; // overlap an S*-character
		}

		++lms_size; // include current L-type

		last_ch = cur_ch, last_t = cur_t;
	}

	if (block_info.is_empty() == true) { // current block is empty, and there's no remaining S*-substrs, then it is leftmost block

		block_info.m_size = lms_end_pos + 1; // lms_end_pos - 0 + 1

		block_info.close();

		m_blocks_info.push_back(block_info);
	}
	else { // current block is not empty, close it and create a new block to carry the remaining characters

		block_info.close();

		m_blocks_info.push_back(block_info);

		block_info = BlockInfo(lms_end_pos, m_blocks_info.size());

		block_info.m_size = lms_end_pos + 1;

		block_info.close();

		m_blocks_info.push_back(block_info);
	}

	compute_block_id_of_samplings(); // compute block_id for samplings

#ifdef TEST_DEBUG2
	std::cerr << "lms_num: " << lms_num << " l_num: " << l_num << " s_num: " << s_num << std::endl;
#endif

	return lms_num;
}

/// \brief compute samplings' block_id
///
/// \note samplings are positioned at 0, 1 * m_capacity, 2 * m_capacity, ... m_s_len - 1
template<typename alphabet_type, typename offset_type, uint8 D>
void DSAComputation<alphabet_type, offset_type, D>::compute_block_id_of_samplings() {

	uint64 sampling_pos = 0;

	uint8 block_id = m_blocks_info.size() - 1;

	while (sampling_pos < m_s_len - 1) {

		while (true) {

			if (sampling_pos <= m_blocks_info[block_id].m_end_pos) {

				m_block_id_of_samplings.push_back(block_id);

				break;
			}

			--block_id;
		}

		sampling_pos += m_blocks_info[0].m_capacity;
	}

	sampling_pos = m_s_len - 1; // the sentinel is also a sampling

	while (true) {

		if (sampling_pos <= m_blocks_info[block_id].m_end_pos) {

			m_block_id_of_samplings.push_back(block_id);

			break;
		}

		--block_id;
	}

	return;
}

/// \brief return the block idx of the character at the given position
///
template<typename alphabet_type, typename offset_type, uint8 D>
uint8 DSAComputation<alphabet_type, offset_type, D>::getBlockId(const offset_type & _pos) {

	uint8 sampling_id = static_cast<uint8>(_pos / m_blocks_info[0].m_capacity); // _pos >= sampling_id * m_capacity

	uint8 block_id = m_block_id_of_samplings[sampling_id];

	while (true) {

		if (m_blocks_info[block_id].m_end_pos >= _pos) {

			return block_id;
		}

		--block_id;
	}

	std::cerr << "getBlockId is wrong\n";

	return 0;
}

/// \brief sort S*-substrs in the block
///
template<typename alphabet_type, typename offset_type, uint8 D>
void DSAComputation<alphabet_type, offset_type, D>::sortSStarBlock(const BlockInfo & _block_info) {

	if (_block_info.is_single() == true) {

		sortSStarSingleBlock(_block_info);
	}

	if (_block_info.is_multi() == true) {

		if (sizeof(alphabet_type) <= sizeof(uint32)) { // no need to format input

			sortSStarMultiBlock<false>(_block_info);
		}
		else {
		
			sortSStarMultiBlock<true>(_block_info);
		}
	}

	return;
}

/// \brief contain single S*-substr. If long, output the heading D characters; otherwise, output the whole substr
///
template<typename alphabet_type, typename offset_type, uint8 D>
void DSAComputation<alphabet_type, offset_type, D>::sortSStarSingleBlock(const BlockInfo & _block_info) {

	typename alphabet_vector_type::const_iterator s_it = m_s->begin() + _block_info.m_beg_pos;

	alphabet_vector_type *short_ch_seq = new alphabet_vector_type(), *long_ch_seq = new alphabet_vector_type();

	uint8_vector_type *short_len_seq = new uint8_vector_type();

	uint32_vector_type *short_pos_seq = new uint32_vector_type();

	if (_block_info.m_size > D) { // is long
		
		// store the heading D characters
		for (uint8 i = 0; i < D; ++i, ++s_it) {

			long_ch_seq->push_back(*s_it);
		}

		// directly read bwt from m_s, no need to compute bwt
	}
	else {

		// store the whole substr
		short_len_seq->push_back(static_cast<uint8>(_block_info.m_size)); // len

		for (uint8 i = 0; i < _block_info.m_size; ++i, ++s_it) {

			short_ch_seq->push_back(*s_it); // characters
		}

		short_pos_seq->push_back(uint32(0)); // starting position (relative)
	}

	m_short_ch_seqs.push_back(short_ch_seq), m_short_len_seqs.push_back(short_len_seq), m_short_pos_seqs.push_back(short_pos_seq);

	m_long_ch_seqs.push_back(long_ch_seq), m_long_l_bwt_seqs.push_back(nullptr), m_long_s_bwt_seqs.push_back(nullptr);

	return;
}

/// \brief contain multiple S*-substr. 
///
/// \note If alphabet_type > uint32, then format the block before inducing in RAM
template<typename alphabet_type, typename offset_type, uint8 D>
template<bool FORMAT>
void DSAComputation<alphabet_type, offset_type, D>::sortSStarMultiBlock(const BlockInfo & _block_info) {
	
	// load s into RAM, format s
	uint32 block_size = _block_info.m_size;

	alphabet_type *s = new alphabet_type[block_size];

	uint32 *fs = nullptr;

	uint32 max_alpha;

	if (FORMAT == true) {

		fs = new uint32[block_size];

		formatMultiBlock(_block_info, s, fs, max_alpha);
	}
	else {

		typename alphabet_vector_type::const_iterator it = m_s->begin() + _block_info.m_beg_pos;

		for (uint32 i = 0; i < block_size; ++i, ++it) {
		
			s[i] = *it;
		}

		max_alpha = ALPHA_MAX;
	}

	char *t_buf = new char[block_size / 8 + 1]; BitWrapper t(t_buf); // L & S type array
	
	char *h_buf = new char[block_size / 8 + 1]; BitWrapper h(h_buf); // long & short type array

	uint32 *sa = new uint32[block_size]; 

	for (uint32 i = 0; i < block_size; ++i) {

		sa[i] = UINT32_MAX;
	}

	uint32 *bkt = new uint32[max_alpha + 1];

	if (FORMAT == true) {
	
		getBuckets<uint32>(fs, block_size, bkt, max_alpha + 1, true);
	}
	else {
		
		getBuckets<alphabet_type>(s, block_size, bkt, max_alpha + 1, true);
	}

	// compute t and h
	uint32 lms_end_pos, lms_beg_pos, lms_size;

	lms_end_pos = block_size - 1, t.set(block_size - 1, S_TYPE), lms_size = 1; // rightmost is the sentinel

	if (FORMAT == true) {

		sa[bkt[fs[lms_end_pos]]--] = lms_end_pos;
	}
	else {

		sa[bkt[s[lms_end_pos]]--] = lms_end_pos;
	}

	t.set(block_size - 2, L_TYPE), ++lms_size; // next on the left is L-type

	for (uint32 i = block_size - 3; ; --i) {

		t.set(i, (s[i] < s[i + 1] || (s[i] == s[i + 1] && t.get(i + 1) == S_TYPE)) ? S_TYPE : L_TYPE);

		if (t.get(i) == L_TYPE && t.get(i + 1) == S_TYPE) { // s[i + 1] is S*-type

			lms_beg_pos = lms_end_pos - lms_size + 1;

			if (lms_size > D) { // long

				for (uint32 j = lms_beg_pos; j < lms_end_pos; ++j) { // set long except for the rightmost
				
					h.set(j, LONG_TYPE);
				}
			}
			else { // short
				
				for (uint32 j = lms_beg_pos; j < lms_end_pos; ++j) { // set short except for the rightmost
			
					h.set(j, SHORT_TYPE);
				}
			}

			lms_end_pos = lms_beg_pos;

			lms_size = 1; // overlap an S*-character

			if (FORMAT == true) {

				sa[bkt[fs[lms_end_pos]]--] = lms_end_pos;
			}
			else {
			
				sa[bkt[s[lms_end_pos]]--] = lms_end_pos;
			}
		}

		++lms_size; // include current character

		if (i == 0) break;
	}

	lms_beg_pos = lms_end_pos - lms_size + 1; // leftmost is S*-type

	if (lms_size > D) {

		for (uint32 j = lms_beg_pos; j < lms_end_pos; ++j) {

			h.set(j, LONG_TYPE);
		}
	}
	else {

		for (uint32 j = lms_beg_pos; j < lms_end_pos; ++j) {

			h.set(j, SHORT_TYPE);
		}
	}	

	// no need to insert the leftmost character (also S*-type) into sa, because it induces no preceding substr in current block
	
	// induce L-type substrs
	alphabet_vector_type *short_ch_seq = new alphabet_vector_type(), *long_ch_seq = new alphabet_vector_type();

	uint8_vector_type *short_len_seq = new uint8_vector_type();

	uint32_vector_type *short_pos_seq = new uint32_vector_type();

	alphabet_vector_type *long_l_bwt_seq = new alphabet_vector_type(), *long_s_bwt_seq = new alphabet_vector_type();

	if (FORMAT == true) {

		getBuckets<uint32>(fs, block_size, bkt, max_alpha + 1, false);
	}
	else {

		getBuckets<alphabet_type>(s, block_size, bkt, max_alpha + 1, false);
	}

	for (uint32 i = 0; i < block_size; ++i) {

		if (sa[i] != UINT32_MAX) { // sa[i] != 0

			uint32 j = sa[i] - 1;

			if (h.get(j) == LONG_TYPE) {

				long_l_bwt_seq->push_back(s[j]); // h.get(j) not h.get(sa[i])
			}

			if (t.get(j) == L_TYPE) {

				if (FORMAT == true) {

					sa[bkt[fs[j]]++] = j;
				}
				else {
		
					sa[bkt[s[j]]++] = j;
				}

				sa[i] = UINT32_MAX; // only reserve L*-type
			}
		}
	} 

	if (FORMAT == true) {
		
		getBuckets<uint32>(fs, block_size, bkt, max_alpha + 1, true);
	}
	else {
		
		getBuckets<alphabet_type>(s, block_size, bkt, max_alpha + 1, true);
	}

	for (uint32 i = block_size - 1; ; --i) {

		if (sa[i] != UINT32_MAX && sa[i] != 0) {

			uint32 j = sa[i] - 1;

			if (h.get(sa[i]) == LONG_TYPE) {

				long_s_bwt_seq->push_back(s[j]); // h.get(sa[i]) not h.get(j)
			}

			if (t.get(j) == S_TYPE) {

				if (FORMAT == true) {

					sa[bkt[fs[j]]--] = j;
				}
				else {

					sa[bkt[s[j]]--] = j;
				}

				sa[i] = UINT32_MAX; // only reserve S*-type
			}
		}

		if (i == 0) break;
	}

	delete bkt; bkt = nullptr;

	// output the sorted S*-substrs
	for (uint32 i = 0; i < block_size; ++i) {

		if (sa[i] != UINT32_MAX) {

			if (h.get(sa[i]) == LONG_TYPE) {
	
				for (uint8 j = 0; j < D; ++j) {

					long_ch_seq->push_back(s[sa[i] + j]);
				}
			}
			else {

				// find next S*-ch
				uint32 j = sa[i] + 1;

				while (t.get(j) == L_TYPE || t.get(j - 1) == S_TYPE) {

					++j;
				}

				uint8 len = static_cast<uint8>(j - sa[i] + 1);

				short_len_seq->push_back(len);

				for (uint32 k = sa[i]; k <= j; ++k) {

					short_ch_seq->push_back(s[k]);
				}

				short_pos_seq->push_back(sa[i]);	
			}
		}
	}

	m_short_ch_seqs.push_back(short_ch_seq), m_short_len_seqs.push_back(short_len_seq), m_short_pos_seqs.push_back(short_pos_seq);

	m_long_ch_seqs.push_back(long_ch_seq), m_long_l_bwt_seqs.push_back(long_l_bwt_seq), m_long_s_bwt_seqs.push_back(long_s_bwt_seq);

	delete [] s; s = nullptr;

	delete [] fs; fs = nullptr;

	delete [] t_buf; t_buf = nullptr;

	delete [] h_buf; h_buf = nullptr;

	delete [] sa;

	return;
}

/// \brief get bucket end/start
///
template<typename alphabet_type, typename offset_type, uint8 D>
template<typename U>
void DSAComputation<alphabet_type, offset_type, D>::getBuckets(const U *_s, const uint32 _s_size, uint32 *_bkt, const uint32 _bkt_num, const bool _end) {

	uint32 sum = 0;

	for (uint32 i = 0; i < _bkt_num; ++i) {

		_bkt[i] = 0;
	}

	for (uint32 i = 0; i < _s_size; ++i) {

		++_bkt[_s[i]];
	}

	for (uint32 i = 0; i < _bkt_num; ++i) {

		sum += _bkt[i]; 
	
		_bkt[i] = _end ? (sum - 1) : (sum - _bkt[i]);
	}

	return;
}

/// \brief format block
///
template<typename alphabet_type, typename offset_type, uint8 D>
void DSAComputation<alphabet_type, offset_type, D>::formatMultiBlock(const BlockInfo & _block_info, alphabet_type *_block, uint32 *_fblock, uint32 & _max_alpha){

	uint32 block_size = _block_info.m_size;

	typename alphabet_vector_type::const_iterator it = m_s->begin() + _block_info.m_beg_pos;

	typedef Pair<alphabet_type, uint32> pair_type;

	typedef TupleAscCmp1<pair_type> pair_comparator_type;

	std::vector<pair_type> container(block_size);

	for (uint32 i = 0; i < block_size; ++i, ++it) {

		container[i] = pair_type(*it, i);
	}

	std::sort(container.begin(), container.end(), pair_comparator_type());

	alphabet_type pre_ch = container[0].first;

	uint32 pre_name = 0;

	_fblock[container[0].second] = pre_name;

	_block[container[0].second] = pre_ch;

	for (uint32 i = 1; i < block_size; ++i) {

		if (container[i].first != pre_ch) {

			pre_ch = container[i].first;

			pre_name = pre_name + 1;
		}

		_fblock[container[i].second] = pre_name;

		_block[container[i].second] = pre_ch;
	}	

	_max_alpha = pre_name;

	return;
}

/// \brief sort long S*-substrs using external memory
///
template<typename alphabet_type, typename offset_type, uint8 D>
void DSAComputation<alphabet_type, offset_type, D>::sortLongSStarGlobal() {

	// sort the ending characters of S*-substrs
	typedef Pair<alphabet_type, offset_type> pair_type;

	typedef TupleAscCmp2<pair_type> pair_comparator_type; // sort by <ch, pos> in ascending order

	typedef typename ExSorter<pair_type, pair_comparator_type>::sorter sorter_type;

	sorter_type *sorter_lms = new sorter_type(pair_comparator_type(), MAX_MEM / 2);

	typename alphabet_vector_type::bufreader_reverse_type s_rev_reader(*m_s);

	{
	
		uint8 cur_t, last_t;

		alphabet_type cur_ch, last_ch, lms_end_ch;

		uint64 lms_size, lms_end_pos;

#ifdef TEST_DEBUG2

		uint64 lms_num = 0;
#endif

		lms_end_pos = m_s->size() - 1, lms_end_ch = *s_rev_reader, lms_size = 1, ++s_rev_reader; // rightmost is sentinel

		cur_ch = *s_rev_reader, cur_t = L_TYPE, ++lms_size, ++s_rev_reader; // next on the left is L-type

		last_ch = cur_ch, last_t = cur_t;

		for (; !s_rev_reader.empty(); ++s_rev_reader) {

			cur_ch = *s_rev_reader;

			cur_t = (cur_ch < last_ch || (cur_ch == last_ch && last_t == S_TYPE)) ? S_TYPE : L_TYPE;

			if (cur_t == L_TYPE && last_t == S_TYPE) { // find an S*-substr

				if (lms_size > D) {

					sorter_lms->push(pair_type(lms_end_ch, lms_end_pos));
				}

#ifdef TEST_DEBUG2

				++lms_num;
#endif
				lms_end_pos = lms_end_pos - lms_size + 1, lms_end_ch = last_ch,	lms_size = 1;
			}	
	
			++lms_size;
	
			last_ch = cur_ch, last_t = cur_t;
		}

		sorter_lms->sort();

#ifdef TEST_DEBUG1

	std::cerr << "lms num: " << lms_num << " long lms num: " << sorter_lms->size() << std::endl;
#endif

	}

	// induce L-type substrs
	typedef Triple<alphabet_type, offset_type, offset_type> triple_type;

	typedef TupleDscCmp3<triple_type> triple_comparator_type; // sort <ch, rank, pos> by the first three components

	typedef typename ExHeap<triple_type, triple_comparator_type, MAX_MEM / 4, MAX_ITEM>::heap heap_type;

	typedef typename heap_type::block_type heap_block_type;

	typedef typename stxxl::read_write_pool<heap_block_type> heap_pool_type;

	heap_pool_type * pq_l_pool = new heap_pool_type(MAX_MEM / 8 / heap_block_type::raw_size, MAX_MEM / 8 / heap_block_type::raw_size);

	heap_type *pq_l = new heap_type(*pq_l_pool);

	alphabet_vector_type *sorted_l_ch = new alphabet_vector_type();

	offset_vector_type *sorted_l_pos = new offset_vector_type();

	uint8_vector_type *sorted_l_diff = new uint8_vector_type();

	std::vector<typename alphabet_vector_type::const_iterator> long_l_bwt_its(m_blocks_info.size());

	std::vector<typename alphabet_vector_type::const_reverse_iterator> s_rits(m_blocks_info.size());

	for (uint8 i = 0; i < m_blocks_info.size(); ++i) {

		if (m_blocks_info[i].is_single()) {

			s_rits[i] = m_s->rbegin() + m_s->size() - m_blocks_info[i].m_end_pos;
		}

		if (m_blocks_info[i].is_multi()) {
	
			long_l_bwt_its[i] = m_long_l_bwt_seqs[i]->begin();
		}
	}

	{

		alphabet_type cur_bkt = ALPHA_MAX;

		offset_type rank_cnt = 0;

		alphabet_type last_l_ch, last_s_ch, last_lml_ch;

		offset_type last_l_rank0, last_l_rank1, last_s_rank0, last_s_rank1, last_lml_rank0;

		last_l_rank0 = last_l_rank1 = last_s_rank0 = last_s_rank1 = last_lml_rank0 = OFFSET_MAX;

		offset_type cur_l_rank0, cur_l_rank1, cur_s_rank0, cur_s_rank1;

		cur_l_rank0 = cur_l_rank1 = cur_s_rank0 = cur_s_rank1 = OFFSET_MAX;

		last_l_ch = last_lml_ch = ALPHA_MIN;

		last_s_ch = ALPHA_MAX;

		if (!sorter_lms->empty()) { // may contain no long S*-substrs

			cur_bkt = (*sorter_lms)->first; 
		}

		while (true) {

			while (!pq_l->empty() && pq_l->top().first == cur_bkt) {

				const triple_type cur_str = pq_l->top();

				cur_l_rank1 = (cur_bkt == last_l_ch && cur_str.second == last_l_rank0) ? last_l_rank1 : rank_cnt;

				rank_cnt = rank_cnt + 1;

				uint8 block_id = getBlockId(cur_str.third);

				alphabet_type pre_ch;

				if (m_blocks_info[block_id].is_multi()) {

					pre_ch = *long_l_bwt_its[block_id], ++long_l_bwt_its[block_id];	
				}
				else {

					pre_ch = *s_rits[block_id], ++s_rits[block_id];
				}

				if (pre_ch >= cur_bkt) {

					pq_l->push(triple_type(pre_ch, cur_l_rank1, cur_str.third - 1));
				}
				else { // L*-type

					(*sorted_l_ch).push_back(cur_bkt);

					(*sorted_l_pos).push_back(cur_str.third);

					(*sorted_l_diff).push_back((cur_bkt == last_lml_ch && cur_str.second == last_lml_rank0) ? 0 : 1);

					last_lml_ch = cur_bkt, last_lml_rank0 = cur_str.second;
				}

				last_l_ch = cur_bkt, last_l_rank0 = cur_str.second, last_l_rank1 = cur_l_rank1;

				pq_l->pop();
			}

			while (!sorter_lms->empty() && (*sorter_lms)->first == cur_bkt) {

				const pair_type cur_str = *(*sorter_lms);

				cur_s_rank1 = (cur_bkt == last_s_ch) ? last_s_rank1 : rank_cnt;

				rank_cnt = rank_cnt + 1;

				uint8 block_id = getBlockId(cur_str.second);

				alphabet_type pre_ch;

				if (m_blocks_info[block_id].is_multi()) {

					pre_ch = *long_l_bwt_its[block_id]; ++long_l_bwt_its[block_id];
				}
				else {

					pre_ch = *s_rits[block_id], ++s_rits[block_id];
				}

				pq_l->push(triple_type(pre_ch, cur_s_rank1, cur_str.second - 1));

				last_s_ch = cur_bkt, last_s_rank1 = cur_s_rank1;

				++(*sorter_lms);
			}

			if (!pq_l->empty()) {

				cur_bkt = pq_l->top().first;

				if (!sorter_lms->empty() && (*sorter_lms)->first < cur_bkt) {

					cur_bkt = (*sorter_lms)->first;
				}
			}
			else if (!sorter_lms->empty()) {

				cur_bkt = (*sorter_lms)->first;
			}
			else {

				break;
			}
		}
	}

	delete sorter_lms; sorter_lms = nullptr;

	delete pq_l_pool; pq_l_pool = nullptr;

	delete pq_l; pq_l = nullptr;

	for (uint8 i = 0; i < m_blocks_info.size() - 1; ++i) {

#ifdef TEST_DEBUG2

		if (m_blocks_info[i].is_multi()) {

			assert(long_l_bwt_its[i] == m_long_l_bwt_seqs[i]->end());
		}
		
#endif
		delete m_long_l_bwt_seqs[i]; m_long_l_bwt_seqs[i] = nullptr; 
	}


	// induce S-type substrs 
	typedef TupleAscCmp3<triple_type> triple_comparator_type2;

	typedef typename ExHeap<triple_type, triple_comparator_type2, MAX_MEM / 4, MAX_ITEM>::heap heap_type2;

	typedef typename heap_type2::block_type heap_block_type2;

	typedef typename stxxl::read_write_pool<heap_block_type2> heap_pool_type2;

	heap_pool_type2 *pq_s_pool = new heap_pool_type2(MAX_MEM / 8 / heap_block_type2::raw_size, MAX_MEM / 8 / heap_block_type2::raw_size);

	heap_type2 *pq_s = new heap_type2(*pq_s_pool);

	m_long_aux_seq = new uint8_vector_type();

	m_long_pos_seq = new offset_vector_type();

	typename alphabet_vector_type::bufreader_reverse_type sorted_l_ch_rev_reader(*sorted_l_ch);

	typename offset_vector_type::bufreader_reverse_type sorted_l_pos_rev_reader(*sorted_l_pos);

	typename uint8_vector_type::bufreader_reverse_type sorted_l_diff_rev_reader(*sorted_l_diff);

	std::vector<typename alphabet_vector_type::const_iterator> long_s_bwt_its(m_blocks_info.size());

	for (uint8 i = 0; i < m_blocks_info.size(); ++i) {

		if (m_blocks_info[i].is_single()) {
		
			--s_rits[i]; // rewind one position, point to the preceding of lefmost L-type character
		}

		if (m_blocks_info[i].is_multi()) {

			long_s_bwt_its[i] = m_long_s_bwt_seqs[i]->begin();			
		}
	}

	{

		alphabet_type cur_bkt = ALPHA_MIN;

		offset_type rank_cnt = m_s->size();

		alphabet_type last_s_ch, last_lms_ch;

		offset_type last_l_rank0, last_l_rank1, last_s_rank0, last_s_rank1, last_lms_rank0;

		last_l_rank0 = last_l_rank1 = last_s_rank0 = last_s_rank1 = last_lms_rank0 = OFFSET_MAX;

		offset_type cur_l_rank0, cur_l_rank1, cur_s_rank0, cur_s_rank1;

		cur_l_rank0 = cur_l_rank1 = cur_s_rank0 = cur_s_rank1 = OFFSET_MAX;

		bool last_l_diff = true;

		last_s_ch = last_lms_ch = ALPHA_MAX; // s_ch != ALPHA_MAX

		if (!sorted_l_ch_rev_reader.empty()) { // may be empty if contains no long

			cur_bkt = *sorted_l_ch_rev_reader;
		}

		while (true) {

			while (!pq_s->empty() && pq_s->top().first == cur_bkt) {

				const triple_type cur_str = pq_s->top();

				cur_s_rank1 = (cur_bkt == last_s_ch && cur_str.second == last_s_rank0) ? last_s_rank1 : rank_cnt;

				rank_cnt = rank_cnt - 1;

				uint8 block_id = getBlockId(cur_str.third);

				if (m_blocks_info[block_id].m_end_pos == cur_str.third) { // leftmost in a single/multi-block

					// in this case, the S*-substr starting with s[cur_str.third] stored in the (block_id - 1)-th block
					uint8 aux = (cur_bkt == last_lms_ch && cur_str.second == last_lms_rank0) ? 0 : 128;

					aux = aux + (block_id - 1);

					m_long_aux_seq->push_back(aux);

					m_long_pos_seq->push_back(cur_str.third);

					last_lms_rank0 = cur_str.second, last_lms_ch = cur_bkt;
				}
				else {

					alphabet_type pre_ch;

					if (m_blocks_info[block_id].is_multi()) {

						pre_ch = *long_s_bwt_its[block_id], ++long_s_bwt_its[block_id];
					}
					else {

						pre_ch = *s_rits[block_id], ++s_rits[block_id];
					}

					if (pre_ch <= cur_bkt) { // preceding is S-type

						pq_s->push(triple_type(pre_ch, cur_s_rank1, cur_str.third - 1));
					}
					else { // current is S*-type

						// in this case, the S*-substr is stored in block_id-th block
						uint8 aux = (cur_bkt == last_lms_ch && cur_str.second == last_lms_rank0) ? 0 : 128;

						aux = aux + block_id;

						m_long_aux_seq->push_back(aux);

						m_long_pos_seq->push_back(cur_str.third);

						last_lms_rank0 = cur_str.second, last_lms_ch = cur_bkt;
					}
				}

				last_s_rank1 = cur_s_rank1, last_s_rank0 = cur_str.second, last_s_ch = cur_bkt;

				pq_s->pop();
			}

			while (!sorted_l_ch_rev_reader.empty() && *sorted_l_ch_rev_reader == cur_bkt) {

				cur_l_rank1 = (last_l_diff == false) ? last_l_rank1 : rank_cnt;

				rank_cnt = rank_cnt - 1;

				uint8 block_id = getBlockId(*sorted_l_pos_rev_reader);

				alphabet_type pre_ch;

				if (m_blocks_info[block_id].is_multi()) {

					pre_ch = *long_s_bwt_its[block_id], ++long_s_bwt_its[block_id];
				}
				else {

					pre_ch = *s_rits[block_id], ++s_rits[block_id];
				}

				pq_s->push(triple_type(pre_ch, cur_l_rank1, *sorted_l_pos_rev_reader - 1));

				last_l_rank1 = cur_l_rank1, last_l_diff = *sorted_l_diff_rev_reader;

				++sorted_l_ch_rev_reader, ++sorted_l_pos_rev_reader, ++sorted_l_diff_rev_reader;
			}

			if (!pq_s->empty()) {

				cur_bkt = pq_s->top().first;

				if (!sorted_l_ch_rev_reader.empty() && *sorted_l_ch_rev_reader > cur_bkt) {

					cur_bkt = *sorted_l_ch_rev_reader;
				}
			}
			else if (!sorted_l_ch_rev_reader.empty()) {

				cur_bkt = *sorted_l_ch_rev_reader;
			}
			else {

				break;
			}
		}	
	}

	// clear
	for (uint8 i = 0; i < m_blocks_info.size() - 1; ++i) {

#ifdef TEST_DEBUG2

		if (m_blocks_info[i].is_multi()) {

			assert(long_s_bwt_its[i] == m_long_s_bwt_seqs[i]->end());
		}

		if (m_blocks_info[i].is_single() && m_blocks_info[i].m_size > D) {

			assert(s_rits[i] == m_s->rbegin() + m_s_len - m_blocks_info[i].m_end_pos + m_blocks_info[i].m_size - 1);
		}
#endif

		delete m_long_s_bwt_seqs[i]; m_long_s_bwt_seqs[i] = nullptr;
	}

#ifdef TEST_DEBUG2

	assert(sorted_l_ch_rev_reader.empty());

	assert(sorted_l_pos_rev_reader.empty());

	assert(sorted_l_diff_rev_reader.empty());
#endif

	delete sorted_l_ch; sorted_l_ch = nullptr;
	
	delete sorted_l_pos; sorted_l_pos = nullptr;

	delete sorted_l_diff; sorted_l_diff = nullptr;

	delete pq_s_pool; pq_s_pool = nullptr;

	delete pq_s; pq_s = nullptr;
}


/// \brief merge sorted short and long S*-substrs using a var-length string sorter
///
template<typename alphabet_type, typename offset_type, uint8 D>
bool DSAComputation<alphabet_type, offset_type, D>::mergeSortedSStarGlobal() {

	std::vector<uint64> blocks_beg_pos;

	for (uint8 i = 0; i < m_blocks_info.size() - 1; ++i) { // exclude the leftmost block (contains no S*-substrs)

		blocks_beg_pos.push_back(m_blocks_info[i].m_beg_pos);
	}

	typedef SubstrSorter<alphabet_type, offset_type, D> str_sorter_type;

	str_sorter_type str_sorter = str_sorter_type(m_short_ch_seqs, m_short_len_seqs, m_short_pos_seqs, 
						m_long_ch_seqs, m_long_pos_seq, m_long_aux_seq, blocks_beg_pos);

	bool is_unique = str_sorter.process(m_s1, m_s->size());

	// clear
	delete m_long_aux_seq; m_long_aux_seq = nullptr;

	delete m_long_pos_seq; m_long_pos_seq = nullptr;

	for (uint8 i = 0; i < m_blocks_info.size() - 1; ++i) {

		delete m_short_ch_seqs[i]; m_short_ch_seqs[i] = nullptr;

		delete m_short_len_seqs[i]; m_short_len_seqs[i] = nullptr;

		delete m_short_pos_seqs[i]; m_short_pos_seqs[i] = nullptr;

		delete m_long_ch_seqs[i]; m_long_ch_seqs[i] = nullptr;	
	}

	return is_unique;
}

/// \brief sort suffix
///
template<typename alphabet_type, typename offset_type, uint8 D>
void DSAComputation<alphabet_type, offset_type, D>::sortSuffixGlobal(offset_vector_type *& _sa1_reverse) {

	uint8 block_num = m_blocks_info.size();

	if (_sa1_reverse == nullptr) {

#ifdef TEST_DEBUG1
		test_output_vector<offset_vector_type>(m_s1);
#endif

		typename offset_vector_type::bufreader_reverse_type s1_rev_reader(*m_s1);

		for (uint8 i = 0; i < block_num; ++i) {

			offset_vector_type *lms_rank_seq = new offset_vector_type();

			if (i == block_num - 1) { // leftmost block
	
				lms_rank_seq->push_back(*s1_rev_reader), ++s1_rev_reader;
			}
			else { // single/multi-block

				for (uint32 j = 0; j < m_blocks_info[i].m_lms_num; ++j, ++s1_rev_reader) {

					lms_rank_seq->push_back(*s1_rev_reader);
				}

			}

			m_lms_rank_seqs.push_back(lms_rank_seq);
		}

#ifdef TEST_DEBUG2

		assert(s1_rev_reader.empty());
#endif

		delete m_s1; m_s1 = nullptr;
	}
	else {

		typedef Pair<offset_type, offset_type> pair_type;

		typedef TupleDscCmp1<pair_type> pair_comparator_type;

		typedef typename ExSorter<pair_type, pair_comparator_type>::sorter sorter_type;

		typename offset_vector_type::bufreader_reverse_type sa1_reader(*_sa1_reverse);

		sorter_type *sorter_lms = new sorter_type(pair_comparator_type(), MAX_MEM);

		for (offset_type i = 0; !sa1_reader.empty(); ++i, ++sa1_reader) { 

			sorter_lms->push(pair_type(*sa1_reader + 1, i)); // sa1_reader + 1 != 0 
		}

		sorter_lms->sort();

		delete _sa1_reverse; _sa1_reverse = nullptr;

		for (uint8 i = 0; i < block_num; ++i) {

			offset_vector_type *lms_rank_seq = new offset_vector_type();

			if (i == block_num - 1) {

				lms_rank_seq->push_back((*sorter_lms)->second); ++(*sorter_lms);
			}
			else {

				for (uint32 j = 0; j < m_blocks_info[i].m_lms_num; ++j, ++(*sorter_lms)) {

					lms_rank_seq->push_back((*sorter_lms)->second);
				}
			}

			m_lms_rank_seqs.push_back(lms_rank_seq);
		}

#ifdef TEST_DEBUG2

		assert(sorter_lms->empty());
#endif

		delete sorter_lms; sorter_lms = nullptr;
	}

	// induce suffixes in each block
	for (uint8 i = 0; i < m_blocks_info.size(); ++i) {

		if (m_blocks_info[i].is_multi()) {

			if (sizeof(alphabet_type) <= sizeof(uint32)) { 
			
				sortSuffixMultiBlock<false>(m_blocks_info[i]); 
			}
			else { 
				
				sortSuffixMultiBlock<true>(m_blocks_info[i]); 
			}
		}
		else if (m_blocks_info[i].is_single()) {

			sortSuffixSingleBlock();
		}
		else {

			sortSuffixLeftmostBlock(m_blocks_info[i]);
		}
	}

	// merge blockwise results 
	mergeSortedSuffixGlobal();

	return;
}

/// \brief sort suffixes in blocks containing one S*-suffixes
///
/// \note directly read bwt from m_s
template<typename alphabet_type, typename offset_type, uint8 D>
void DSAComputation<alphabet_type, offset_type, D>::sortSuffixSingleBlock() {

	m_l_bwt_seqs.push_back(nullptr), m_s_bwt_seqs.push_back(nullptr);

#ifdef TEST_DEBUG2

	m_l_bpos_seqs.push_back(nullptr), m_s_bpos_seqs.push_back(nullptr);
#endif

}

/// \brief sort suffixes in the leftmost block containing no S*-suffixes
///
template<typename alphabet_type, typename offset_type, uint8 D>
void DSAComputation<alphabet_type, offset_type, D>::sortSuffixLeftmostBlock(const BlockInfo & _block_info) {

	alphabet_vector_type *l_bwt_seq = new alphabet_vector_type(), *s_bwt_seq = new alphabet_vector_type();

#ifdef TEST_DEBUG2
	offset_vector_type *l_bpos_seq = new offset_vector_type(), *s_bpos_seq = new offset_vector_type();
#endif

	typename alphabet_vector_type::const_reverse_iterator s_rit = m_s->rbegin() + m_s_len -_block_info.m_end_pos - 1; 

	// induce L-type suffixes
	alphabet_type cur_ch, last_ch;

	uint8 cur_t, last_t;

#ifdef TEST_DEBUG2

	offset_type pos = _block_info.m_end_pos;
#endif

	last_ch = *s_rit, last_t = S_TYPE, ++s_rit; //rightmost is S*-type

#ifdef TEST_DEBUG2
	pos = pos - 1;
#endif
	
	while (s_rit != m_s->rend()) {

		cur_ch = *s_rit;

		cur_t = (cur_ch < last_ch || (cur_ch == last_ch && last_t == S_TYPE)) ? S_TYPE : L_TYPE;

		l_bwt_seq->push_back(cur_ch); // push preceding

#ifdef TEST_DEBUG2
		l_bpos_seq->push_back(pos);
#endif

		if (cur_t == S_TYPE) {

			break;
		}

		last_ch = cur_ch, last_t = cur_t, ++s_rit;	

#ifdef TEST_DEBUG2

		if (0 != pos) pos = pos - 1;
#endif
	}	

	// induce S-type suffixes
	while (s_rit != m_s->rend()) {
		
		s_bwt_seq->push_back(*s_rit); // push preceding

#ifdef TEST_DEBUG2
		s_bpos_seq->push_back(pos);
#endif

		++s_rit;	

#ifdef TEST_DEBUG2
		if(0 != pos) pos = pos - 1;
#endif		
	}

	m_l_bwt_seqs.push_back(l_bwt_seq), m_s_bwt_seqs.push_back(s_bwt_seq);

#ifdef TEST_DEBUG2

	m_l_bpos_seqs.push_back(l_bpos_seq), m_s_bpos_seqs.push_back(s_bpos_seq);

	std::cerr << "block_id: " << _block_info.m_id << " bwt size: " << l_bwt_seq->size() << " " << s_bwt_seq->size() << std::endl;
#endif

	return;
}

/// \brief sort suffixes in blocks containing multiple S*-suffixes
///
template<typename alphabet_type, typename offset_type, uint8 D>
template<bool FORMAT>
void DSAComputation<alphabet_type, offset_type, D>::sortSuffixMultiBlock(const BlockInfo & _block_info) {

	uint32 block_size = _block_info.m_size;

	alphabet_type *s = new alphabet_type[block_size];

	uint32 *fs = nullptr;

	uint32 max_alpha;

	if (FORMAT == true) {

		fs = new uint32[block_size];

		formatMultiBlock(_block_info, s, fs, max_alpha);
	}
	else {

		typename alphabet_vector_type::const_iterator it = m_s->begin() + _block_info.m_beg_pos;

		for (uint32 i = 0; i < block_size; ++i, ++it) {

			s[i] = *it;
		}

		max_alpha = ALPHA_MAX;
	}

	// compute t and sort S*-suffixes
	typedef Pair<offset_type, uint32> pair_type;

	typedef TupleDscCmp1<pair_type> pair_comparator_type;

	uint8 block_id = _block_info.m_id;

	typename offset_vector_type::bufreader_type lms_rank_reader(*m_lms_rank_seqs[block_id]);

	std::vector<pair_type> *container = new std::vector<pair_type>(m_lms_rank_seqs[block_id]->size());

	char *t_buf = new char[block_size / 8 + 1]; BitWrapper t(t_buf);

	t.set(block_size - 1, S_TYPE);

	for (uint32 i = block_size - 2, j = 0; ; --i) {

		t.set(i, (s[i] < s[i + 1] || (s[i] == s[i + 1] && t.get(i + 1) == S_TYPE)) ? S_TYPE : L_TYPE);

		if (t.get(i) == L_TYPE && t.get(i + 1) == S_TYPE) { // s[i + 1] is S*-type
	
			(*container)[j++] = pair_type(*lms_rank_reader, i + 1);

			 ++lms_rank_reader;
		}

		if (i == 0) break;
	}

#ifdef TEST_DEBUG2
	assert(lms_rank_reader.empty());
#endif

	std::sort(container->begin(), container->end(), pair_comparator_type());

#ifdef TEST_DEBUG2

	assert(_block_info.m_lms_num == container->size());

	for (uint32 i = 0; i < _block_info.m_lms_num - 1; ++i) {

		if ((*container)[i].first <= (*container)[i + 1].first) {

			std::cerr << "container[i]: " << (*container)[i].first << " container[i + 1]: " << (*container)[i + 1].first << std::endl;
			
			exit(0);
		}
	}
#endif
	// insert sorted S*-suffixes
	uint32 * sa = new uint32[block_size];

	for (uint32 i = 0; i < block_size; ++i) {

		sa[i] = UINT32_MAX;
	}

	uint32 * bkt = new uint32[max_alpha + 1];

	if (FORMAT == true) {

		getBuckets<uint32>(fs, block_size, bkt, max_alpha + 1, true);
	}
	else {
	
		getBuckets<alphabet_type>(s, block_size, bkt, max_alpha + 1, true);
	}

	for (uint32 i = 0; i < m_lms_rank_seqs[block_id]->size(); ++i) {

		if (FORMAT == true) {

			sa[bkt[fs[(*container)[i].second]]--] = (*container)[i].second;
		}
		else {

			sa[bkt[s[(*container)[i].second]]--] = (*container)[i].second;
		}
	}

	delete container; container = nullptr;

	// induce L-type suffixes
	alphabet_vector_type *l_bwt_seq = new alphabet_vector_type();

#ifdef TEST_DEBUG2
	offset_vector_type *l_bpos_seq = new offset_vector_type();
#endif

	if (FORMAT == true) {

		getBuckets<uint32>(fs, block_size, bkt, max_alpha + 1, false);
	}
	else {

		getBuckets<alphabet_type>(s, block_size, bkt, max_alpha + 1, false);
	}

	for (uint32 i = 0; i < block_size; ++i) {

		if (sa[i] != UINT32_MAX) { // sa[i] != 0

			uint32 j = sa[i] - 1;

			l_bwt_seq->push_back(s[j]);

#ifdef TEST_DEBUG2

			l_bpos_seq->push_back(_block_info.m_beg_pos + j);

			if (_block_info.m_beg_pos + sa[i] >= 205097530 && _block_info.m_beg_pos + sa[i] <= 205097540) {

				std::cerr << "1--sa[i]: " << _block_info.m_beg_pos + sa[i] << " j: " << _block_info.m_beg_pos + j << std::endl;

				std::cerr << "s[sa[i]]: " << s[sa[i]] << " s[j]: " << s[j] << std::endl;

				std::cerr << "t[sa[i]]: " << t.get(sa[i]) << " t[j]: " << t.get(j) << std::endl; 
			}
#endif

			if (t.get(j) == L_TYPE) {
			
				if (FORMAT == true) {
			
					sa[bkt[fs[j]]++] = j;
				}
				else {

					sa[bkt[s[j]]++] = j;
				}
			}

			if (t.get(sa[i]) == S_TYPE) {

				sa[i] = UINT32_MAX; // reserve all L-type
			}	
		}
	}

	// induce S-type suffixes
	alphabet_vector_type *s_bwt_seq = new alphabet_vector_type();

#ifdef TEST_DEBUG2
	offset_vector_type *s_bpos_seq = new offset_vector_type();
#endif

	if (FORMAT == true) {

		getBuckets<uint32>(fs, block_size, bkt, max_alpha + 1, true);
	}
	else {

		getBuckets<alphabet_type>(s, block_size, bkt, max_alpha + 1, true);
	}

	for (uint32 i = block_size - 1; ; --i) {

		if (sa[i] != 0 && sa[i] != UINT32_MAX) {

			uint32 j = sa[i] - 1;

			s_bwt_seq->push_back(s[j]);

#ifdef TEST_DEBUG2
			s_bpos_seq->push_back(_block_info.m_beg_pos + j);

			if (_block_info.m_beg_pos + sa[i] >= 205097530 && _block_info.m_beg_pos + sa[i] <= 205097540) {

				std::cerr << "2--sa[i]: " << _block_info.m_beg_pos + sa[i] << " j: " << _block_info.m_beg_pos + j << std::endl;

				std::cerr << "s[sa[i]]: " << s[sa[i]] << " s[j]: " << s[j] << std::endl;

				std::cerr << "t[sa[i]]: " << t.get(sa[i]) << " t[j]: " << t.get(j) << std::endl; 
			}
#endif

			if (t.get(j) == S_TYPE) {

				if (FORMAT == true) {

					sa[bkt[fs[j]]--] = j;
				}
				else {

					sa[bkt[s[j]]--] = j;
				}
			}
		}

		if (i == 0) break;
	}

	m_l_bwt_seqs.push_back(l_bwt_seq), m_s_bwt_seqs.push_back(s_bwt_seq);

#ifdef TEST_DEBUG2

	m_l_bpos_seqs.push_back(l_bpos_seq), m_s_bpos_seqs.push_back(s_bpos_seq);

	std::cerr << "block_id: " << _block_info.m_id << " bwt size: " << l_bwt_seq->size() << " " << s_bwt_seq->size() << std::endl;
#endif

	delete [] bkt; bkt = nullptr;

	delete [] s; s = nullptr;

	delete [] fs; fs = nullptr;

	delete [] t_buf; t_buf = nullptr;

	delete [] sa; sa = nullptr;
}

/// \brief merge block-wise sorted results
///
template<typename alphabet_type, typename offset_type, uint8 D>
void DSAComputation<alphabet_type, offset_type, D>::mergeSortedSuffixGlobal() {

#ifdef TEST_DEBUG2
	std::cerr << "here5\n";
#endif

	// sort S*-suffixes
	typedef Triple<offset_type, alphabet_type, offset_type> triple_type; // <rank, ch, pos>

	typedef TupleAscCmp1<triple_type> triple_comparator_type;

	typedef typename ExSorter<triple_type, triple_comparator_type>::sorter sorter_type;
	
	sorter_type *sorter_lms = new sorter_type(triple_comparator_type(), MAX_MEM / 2);

	for (uint8 i = 0; i < m_blocks_info.size(); ++i) {

		typename alphabet_vector_type::const_reverse_iterator s_rit = m_s->rbegin() + m_s->size() - 1 - m_blocks_info[i].m_end_pos;
	
		typename offset_vector_type::const_iterator lms_rank_it = m_lms_rank_seqs[i]->begin();

		if (m_blocks_info[i].is_multi()) { // blocks contain multiple S*-suffixes

			alphabet_type cur_ch, last_ch;

			uint8 cur_t, last_t;

			offset_type cur_pos, last_pos;
		
			last_ch = *s_rit, last_t = S_TYPE, last_pos = m_blocks_info[i].m_end_pos, ++s_rit;

			while (lms_rank_it != m_lms_rank_seqs[i]->end()) {

				cur_ch = *s_rit;

				cur_t = (cur_ch < last_ch || (cur_ch == last_ch && last_t == S_TYPE)) ? S_TYPE : L_TYPE;

				if (cur_t == L_TYPE && last_t == S_TYPE) {

					sorter_lms->push(triple_type(*lms_rank_it + 1, last_ch, last_pos)); // lms_rank_it + 1 to avoid 0

					++lms_rank_it;
				} 
			
				last_ch = cur_ch, last_t = cur_t, --last_pos, ++s_rit;
			}
		}
		else { // blocks contain one or no S*-suffixes

			sorter_lms->push(triple_type(*lms_rank_it + 1, *s_rit, m_blocks_info[i].m_end_pos));

			++s_rit, ++lms_rank_it; // lms_rank_it + 1 to avoid 0
		}

		delete m_lms_rank_seqs[i]; m_lms_rank_seqs[i] = nullptr;	
	}

	sorter_lms->sort();


#ifdef TEST_DEBUG2

	std::cerr << "here6";

	std::cerr << "sorter_lms size: " << sorter_lms->size() << std::endl;
#endif

	// induce L-type suffixes
	typedef Triple<alphabet_type, offset_type, offset_type> triple_type2;

	typedef TupleDscCmp2<triple_type2> triple_comparator_type2;

	typedef typename ExHeap<triple_type2, triple_comparator_type2, MAX_MEM / 4, MAX_ITEM>::heap heap_type;

	typedef typename heap_type::block_type heap_block_type;

	typedef typename stxxl::read_write_pool<heap_block_type> heap_pool_type;

	heap_pool_type * pq_l_pool = new heap_pool_type(MAX_MEM / 8 / heap_block_type::raw_size, MAX_MEM / 8 / heap_block_type::raw_size);

	heap_type *pq_l = new heap_type(*pq_l_pool);

	alphabet_vector_type *sorted_l_ch = new alphabet_vector_type();

	offset_vector_type *sorted_l_pos = new offset_vector_type();

	std::vector<typename alphabet_vector_type::const_reverse_iterator> s_rits(m_blocks_info.size());

	std::vector<typename alphabet_vector_type::const_iterator> l_bwt_its(m_blocks_info.size());

#ifdef TEST_DEBUG2
	std::vector<typename offset_vector_type::const_iterator> l_bpos_its(m_blocks_info.size());
#endif

	for (uint8 i = 0; i < m_blocks_info.size(); ++i)  {

		if (m_blocks_info[i].is_single()) {

			s_rits[i] = m_s->rbegin() + m_s->size() - m_blocks_info[i].m_end_pos;
		}
		else {
	
			l_bwt_its[i] = m_l_bwt_seqs[i]->begin();

#ifdef TEST_DEBUG2
			l_bpos_its[i] = m_l_bpos_seqs[i]->begin();
#endif
		}
	}

	{

		alphabet_type cur_bkt = ALPHA_MAX;

		offset_type cur_rank, rank_cnt = 1; // start rank from 1 to avoid <0, 0, s_len - 1>

		cur_bkt = (*sorter_lms)->second; // sorter_lms contains at least the sentinel

		while (true) {

			while (!pq_l->empty() && pq_l->top().first == cur_bkt) {

				const triple_type2 cur_str = pq_l->top();


				cur_rank = rank_cnt;

#ifdef TEST_DEBUG2
				
				if (205097530 <= cur_str.third && 205097540 >= cur_str.third) {

					std::cerr << "lalal--cur_bkt: " << cur_bkt << " last rank: " << cur_str.second << " pos: " << cur_str.third;

					std::cerr << "cur_rank: " << cur_rank << std::endl;
				}
#endif

				rank_cnt = rank_cnt + 1;

				uint8 block_id = getBlockId(cur_str.third);

				if (OFFSET_MIN != cur_str.third) {

					alphabet_type pre_ch;
			
					if (m_blocks_info[block_id].is_single() == false) {

						pre_ch = *l_bwt_its[block_id], ++l_bwt_its[block_id];
					}
					else {

						pre_ch = *s_rits[block_id], ++s_rits[block_id];
					}


#ifdef TEST_DEBUG2
					offset_type pre_pos = *l_bpos_its[block_id];

					++l_bpos_its[block_id];

					if (cur_str.third - 1 != pre_pos) {

						std::cerr << "l-ttype pre_pos correct: " << pre_pos << " wrong: " << cur_str.third - 1 << std::endl;

						exit(0);
					}
#endif

#ifdef TEST_DEBUG2
				
				if (205097530 <= cur_str.third && 205097540 >= cur_str.third) {

					std::cerr << "before--cur_bkt: " << cur_bkt << " last rank: " << cur_str.second << " pos: " << cur_str.third;

					std::cerr << "cur_rank: " << cur_rank << std::endl;
				}
#endif

					if (pre_ch >= cur_bkt) {

						pq_l->push(triple_type2(pre_ch, cur_rank, cur_str.third - 1));
					}

#ifdef TEST_DEBUG2
				
				if (205097530 <= cur_str.third && 205097540 >= cur_str.third) {

					std::cerr << "after--cur_bkt: " << cur_bkt << " last rank: " << cur_str.second << " pos: " << cur_str.third;

					std::cerr << "cur_rank: " << cur_rank << std::endl;
				}
#endif

#ifdef TEST_DEBUG2

					if (205097530 <= cur_str.third && 205097540 >= cur_str.third) {
	
						std::cerr << "3--pos: " << cur_str.third << " pre_pos: " << pre_pos << std::endl;

						std::cerr << "--ch: " << cur_str.first << " pre_ch: " << pre_ch << std::endl;
				
						std::cerr << "--cur_rank: " << cur_rank << " last_rank: " << cur_str.second << std::endl;
					}
#endif
				}

				(*sorted_l_ch).push_back(cur_str.first);

				(*sorted_l_pos).push_back(cur_str.third);

				pq_l->pop();
			}

			while (!sorter_lms->empty() && (*sorter_lms)->second == cur_bkt) {

				const triple_type cur_str = *(*sorter_lms);

				cur_rank = rank_cnt;
				
				rank_cnt = rank_cnt + 1;

				uint8 block_id = getBlockId(cur_str.third);

				alphabet_type pre_ch;

				if (m_blocks_info[block_id].is_single() == false) {

					pre_ch = *l_bwt_its[block_id], ++l_bwt_its[block_id];
				}
				else {

					pre_ch = *s_rits[block_id], ++s_rits[block_id];
				}

#ifdef TEST_DEBUG2
				offset_type pre_pos = *l_bpos_its[block_id];

				++l_bpos_its[block_id];

				if (cur_str.third - 1 != pre_pos) {

					std::cerr << "lms pre_pos correct: " << pre_pos << " wrong: " << cur_str.third - 1 << std::endl;

					exit(0);
				}
#endif

#ifdef TEST_DEBUG2	
				assert(pre_ch > cur_bkt);

				if (205097530 <= cur_str.third && 205097540 >= cur_str.third) {

					std::cerr << "4--pos: " << cur_str.third << " pre_pos: " << pre_pos << std::endl;

					std::cerr << "--ch: " << cur_str.second << " pre_ch: " << pre_ch << std::endl;
			
					std::cerr << "--cur_rank: " << cur_rank << std::endl;
				}
#endif

				pq_l->push(triple_type2(pre_ch, cur_rank, cur_str.third - 1));

				++(*sorter_lms);
			}

			if (!pq_l->empty()) {

				cur_bkt = pq_l->top().first;

				if (!sorter_lms->empty() && (*sorter_lms)->second < cur_bkt) {

					cur_bkt = (*sorter_lms)->second;
				}
			}
			else if (!sorter_lms->empty()) {

				cur_bkt = (*sorter_lms)->second;
			}
			else {

				break;
			}
		}
	}

	// clear
	delete sorter_lms; sorter_lms = nullptr;

	delete pq_l_pool; pq_l_pool = nullptr;

	delete pq_l; pq_l = nullptr;

	for (uint8 i = 0; i < m_blocks_info.size(); ++i) {

#ifdef TEST_DEBUG2

		if (m_blocks_info[i].is_single() == false) {

			assert(l_bwt_its[i] == m_l_bwt_seqs[i]->end());
		}
#endif
		
		delete m_l_bwt_seqs[i]; m_l_bwt_seqs[i] = nullptr; 
	}

#ifdef TEST_DEBUG2

	std::cerr << "here7\n";

	std::cerr << "sorted_l_ch_size: " << sorted_l_ch->size() << std::endl;

	std::cerr << "sorted_l_pos_size: " << sorted_l_pos->size() << std::endl;
#endif

#ifdef TEST_DEBUG1

	std::cerr << "sorted_l_ch: ";
	test_output_vector<alphabet_vector_type>(sorted_l_ch);

	std::cerr << "sorted_l_pos: ";
	test_output_vector<offset_vector_type>(sorted_l_pos);
#endif
	// induce S-type suffixes
	typedef TupleAscCmp2<triple_type2> triple_comparator_type3;

	typedef typename ExHeap<triple_type2, triple_comparator_type3, MAX_MEM / 4, MAX_ITEM>::heap heap_type2;

	typedef typename heap_type2::block_type heap_block_type2;

	typedef typename stxxl::read_write_pool<heap_block_type2> heap_pool_type2;

	heap_pool_type2 *pq_s_pool = new heap_pool_type2(MAX_MEM / 8 / heap_block_type2::raw_size, MAX_MEM / 8 / heap_block_type2::raw_size);

	heap_type2 *pq_s = new heap_type2(*pq_s_pool);

	typename alphabet_vector_type::bufreader_reverse_type sorted_l_ch_rev_reader(*sorted_l_ch);

	typename offset_vector_type::bufreader_reverse_type sorted_l_pos_rev_reader(*sorted_l_pos);

	std::vector<typename alphabet_vector_type::const_iterator> s_bwt_its(m_blocks_info.size());	

#ifdef TEST_DEBUG2
	std::vector<typename offset_vector_type::const_iterator> s_bpos_its(m_blocks_info.size());
#endif

	std::vector<bool> is_l_star(m_blocks_info.size());

#ifdef TEST_DEBUG2

	std::vector<uint64> bkt_s_cnt(m_blocks_info.size(), 0);

	std::vector<uint64> bkt_l_cnt(m_blocks_info.size(), 0);
#endif

	for (uint8 i = 0; i < m_blocks_info.size(); ++i) {

		if (m_blocks_info[i].is_multi()) { // multi-block

			s_bwt_its[i] = m_s_bwt_seqs[i]->begin();

#ifdef TEST_DEBUG2
			s_bpos_its[i] = m_s_bpos_seqs[i]->begin();
#endif

		}
		else if (m_blocks_info[i].is_single()) { // single-block

			--s_rits[i];

			is_l_star[i] = true; // single-blocks contain at most one L*-type suffix
		}
		else { // leftmost block

			s_bwt_its[i] = m_s_bwt_seqs[i]->begin();

#ifdef TEST_DEBUG2
			s_bpos_its[i] = m_s_bpos_seqs[i]->begin();
#endif

			is_l_star[i] = true;
		}
	}

	{

		m_sa_reverse = new offset_vector_type();

		alphabet_type cur_bkt = ALPHA_MIN;

		offset_type cur_rank, rank_cnt = m_s_len + 1;

#ifdef TEST_DEBUG2

		std::cerr << "m_s_len: " << m_s_len << std::endl;
#endif

		if (!sorted_l_ch_rev_reader.empty()) { // may contain no L-type suffixes

			cur_bkt = *sorted_l_ch_rev_reader;
		}

		while (true) {

			while (!pq_s->empty() && pq_s->top().first == cur_bkt) {

				const triple_type2 cur_str = pq_s->top();

				cur_rank = rank_cnt;

				rank_cnt = rank_cnt - 1;

				uint8 block_id = getBlockId(cur_str.third);

				if (m_blocks_info[block_id].m_end_pos != cur_str.third && OFFSET_MIN != cur_str.third) {

					alphabet_type pre_ch;

#ifdef TEST_DEBUG2
					bkt_s_cnt[block_id]++;
#endif

					if (m_blocks_info[block_id].is_single() == false) {

						pre_ch = *s_bwt_its[block_id], ++s_bwt_its[block_id];
					}
					else {

						pre_ch = *s_rits[block_id], ++s_rits[block_id];
					}

#ifdef TEST_DEBUG2
				offset_type pre_pos = *s_bpos_its[block_id];

				++s_bpos_its[block_id];

				if (cur_str.third - 1 == 205097534 || cur_str.third - 1 == 205097533) {

					std::cerr << "cur_pos: " << cur_str.third << " pre_pos: " << cur_str.third - 1;

					std::cerr << "cur_bkt: " << cur_bkt << " pre_ch: " << pre_ch << std::endl;

					std::cerr << "cur_type: " << "s-type" << std::endl;
				}

				if (cur_str.third - 1 != pre_pos) {

					std::cerr << " s-type pre_pos correct: " << pre_pos << " wrong: " << cur_str.third - 1;

					std::cerr << " pre pre pos: " << *s_bpos_its[block_id];

					std::cerr << " cur_rank: " << cur_rank << " cur_bkt: " << (uint64)cur_bkt << std::endl;

					exit(0);
				}

#endif


					if (pre_ch <= cur_bkt) {

						pq_s->push(triple_type2(pre_ch, cur_rank, cur_str.third - 1));
					}
				}

				m_sa_reverse->push_back(cur_str.third);

				pq_s->pop();
			}

			while (!sorted_l_ch_rev_reader.empty() && *sorted_l_ch_rev_reader == cur_bkt) {

				cur_rank = rank_cnt; 

				rank_cnt = rank_cnt - 1;

				uint8 block_id = getBlockId(*sorted_l_pos_rev_reader);

				if (OFFSET_MIN != *sorted_l_pos_rev_reader) {

					alphabet_type pre_ch;

#ifdef TEST_DEBUG2
					bkt_l_cnt[block_id]++;
#endif

					if (m_blocks_info[block_id].is_multi()) {

						pre_ch = *s_bwt_its[block_id], ++s_bwt_its[block_id];
#ifdef TEST_DEBUG2
				offset_type pre_pos = *s_bpos_its[block_id];

				++s_bpos_its[block_id];

				if (*sorted_l_pos_rev_reader - 1 == 205097534 || *sorted_l_pos_rev_reader - 1 == 205097533) {

					std::cerr << "cur_pos: " << *sorted_l_pos_rev_reader << " pre_pos: " << *sorted_l_pos_rev_reader - 1;

					std::cerr << "cur_bkt: " << cur_bkt << " pre_ch: " << pre_ch << std::endl;

					std::cerr << "cur_type: " << "l-type" << std::endl;
				}

				if (*sorted_l_pos_rev_reader - 1 != pre_pos) {

					std::cerr << " l-type pre_pos correct: " << pre_pos << " wrong: " << *sorted_l_pos_rev_reader - 1;

					std::cerr << " pre pre pos: " << *s_bpos_its[block_id];

					std::cerr << " cur_rank: " << cur_rank << " cur_bkt: " << (uint64)cur_bkt << " pre_ch: " << (uint64)pre_ch << std::endl;

					exit(0);
				}
#endif

						if (pre_ch < cur_bkt) {

							pq_s->push(triple_type2(pre_ch, cur_rank, *sorted_l_pos_rev_reader - 1));
						}
					}
					else if (is_l_star[block_id] == true) { // a single or the leftmost block contains at most one L*-type 

						if (m_blocks_info[block_id].is_single()) {

							pre_ch = *s_rits[block_id], ++s_rits[block_id];
						}
						else {

							pre_ch = *s_bwt_its[block_id], ++s_bwt_its[block_id];
						}

#ifdef TEST_DEBUG2
				offset_type pre_pos = *s_bpos_its[block_id];

				++s_bpos_its[block_id];

				if (*sorted_l_pos_rev_reader - 1 == 205097534 || *sorted_l_pos_rev_reader - 1 == 205097533) {

					std::cerr << "cur_pos: " << *sorted_l_pos_rev_reader << " pre_pos: " << *sorted_l_pos_rev_reader - 1;

					std::cerr << " cur_bkt: " << cur_bkt << " pre_ch: " << pre_ch << std::endl;

					std::cerr << "cur_type: " << "l-type" << " cur_rank: " << cur_rank << std::endl;
				}

				if (*sorted_l_pos_rev_reader - 1 != pre_pos) {

					std::cerr << " l-type pre_pos correct: " << pre_pos << " wrong: " << *sorted_l_pos_rev_reader - 1 << std::endl;

					exit(0);
				}
#endif
						if (pre_ch < cur_bkt) {

							pq_s->push(triple_type2(pre_ch, cur_rank, *sorted_l_pos_rev_reader - 1));

							is_l_star[block_id] = false;
						}
					}
				}

				m_sa_reverse->push_back(*sorted_l_pos_rev_reader);

				++sorted_l_ch_rev_reader, ++sorted_l_pos_rev_reader;
			}

			if (!pq_s->empty()) {

				cur_bkt = pq_s->top().first;

				if (!sorted_l_ch_rev_reader.empty() && *sorted_l_ch_rev_reader > cur_bkt) {

					cur_bkt = *sorted_l_ch_rev_reader;
				}
			}
			else if (!sorted_l_ch_rev_reader.empty()) {

				cur_bkt = *sorted_l_ch_rev_reader;
			}
			else {

				break;
			}
		}
	}

	m_sa_reverse->push_back(offset_type(m_s_len - 1)); // push the sentinel

#ifdef TEST_DEBUG2

	std::cerr << "here8\n";

	for (uint8 i = 0; i < m_blocks_info.size(); ++i) {

		std::cerr << "block_id: " << (uint32)i << " bkt_s_cnt: " << bkt_s_cnt[i] << std::endl;

		std::cerr << "block_id: " << (uint32)i << " bkt_l_cnt: " << bkt_l_cnt[i] << std::endl;
	}

	std::cerr << "sa: " << m_sa_reverse->size() << std::endl;

	assert(sorted_l_ch_rev_reader.empty());

	assert(sorted_l_pos_rev_reader.empty());
#endif

	// clear
	for (uint8 i = 0; i < m_blocks_info.size(); ++i) {

#ifdef TEST_DEBUG2

		std::cerr << "block_id: " << (uint32)i << std::endl;


		if (m_blocks_info[i].is_single() == false) {

		//	assert(s_bwt_its[i] == m_s_bwt_seqs[i]->end());

			std::cerr << "s_bwt_its[i] - m_s_bwt_seqs[i]->end(): " << s_bwt_its[i] - m_s_bwt_seqs[i]->end() << std::endl;
		}
		else {

		//	assert(s_rits[i] == m_s->rbegin() + m_s_len - m_blocks_info[i].m_end_pos + m_blocks_info[i].m_size - 1);

			std::cerr << "s_rits[i] - begin: " << s_rits[i] - m_s->rbegin() + m_s_len - m_blocks_info[i].m_end_pos + m_blocks_info[i].m_size - 1 << std::endl;
		}
#endif

		delete m_s_bwt_seqs[i]; m_s_bwt_seqs[i] = nullptr;
	}


	delete sorted_l_ch; sorted_l_ch = nullptr;

	delete sorted_l_pos; sorted_l_pos = nullptr;

	delete pq_s_pool; pq_s_pool = nullptr;

	delete pq_s; pq_s = nullptr;

#ifdef TEST_DEBUG2

	std::cerr << "here9";
#endif

#ifdef TEST_DEBUG1

	test_output_vector<offset_vector_type>(m_sa_reverse);
#endif
}

#endif
