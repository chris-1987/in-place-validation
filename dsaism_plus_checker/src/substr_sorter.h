#ifndef SUBSTR_SORTER_H
#define SUBSTR_SORTER_H


#include "substr.h"
#include "losertree.h"

#define TEST_DEBUG_SUBSTR_SORTER

template<typename alphabet_type, typename offset_type, uint8 D>
class SubstrSorter{

private:

	typedef typename ExVector<alphabet_type>::vector alphabet_vector_type;

	typedef typename ExVector<offset_type>::vector offset_vector_type;

	typedef typename ExVector<uint8>::vector uint8_vector_type;

	typedef typename ExVector<uint32>::vector uint32_vector_type;

public:

	struct StrCompare{

	public:

		std::vector<alphabet_vector_type*> & m_short_ch_seqs;

		std::vector<typename alphabet_vector_type::const_iterator> m_short_ch_its; ///< short lms seqs

		std::vector<uint8_vector_type*> & m_short_len_seqs;

		std::vector<typename uint8_vector_type::const_iterator> m_short_len_its;

		std::vector<alphabet_vector_type*> & m_long_ch_seqs;

		std::vector<typename alphabet_vector_type::const_iterator> m_long_ch_its; ///< long lms seqs

		uint8 m_seq_num; ///< number of seqs

		std::vector<Substr<alphabet_type, D>> m_head; ///< current smallest for each seq

		std::vector<bool> m_exists; ///< existence of each seq

	public:
		/// \brief ctor
		///
		StrCompare(std::vector<alphabet_vector_type*> &_short_ch_seqs, std::vector<uint8_vector_type*> & _short_len_seqs, std::vector<alphabet_vector_type*> &_long_ch_seqs) : m_short_ch_seqs(_short_ch_seqs), m_short_len_seqs(_short_len_seqs), m_long_ch_seqs(_long_ch_seqs) {

			for (uint8 i = 0; i < m_short_ch_seqs.size(); ++i) {

				m_short_ch_its.push_back(m_short_ch_seqs[i]->begin());

				m_short_len_its.push_back(m_short_len_seqs[i]->begin());

				m_long_ch_its.push_back(m_long_ch_seqs[i]->begin());
			}
		
			m_seq_num = m_short_ch_seqs.size() + 1;

			m_head.resize(m_seq_num);

			m_exists.clear();

			m_exists.resize(m_seq_num, false);
		}	

		/// \brief fecth a short lms-substr
		///
		void fetchShort(const uint8 _seq) {

			m_exists[_seq] = m_head[_seq].deserializeShort(m_short_ch_seqs[_seq], m_short_ch_its[_seq], m_short_len_its[_seq]);
		}

		/// \brief fetch a long lms-substr
		///
		/// \note long
		void fetchLong(const uint8 _block_id) {

			m_head[m_seq_num - 1].deserializeLong(m_long_ch_its[_block_id]);
		}

		/// \brief check existence
		///
		bool exists(const uint8 _seq) const {

			return m_exists[_seq];
		}

		/// \brief compare
		///
		int operator()(const uint8 _seqa, const uint8 _seqb) const {

			return m_head[_seqa].cmp(m_head[_seqb]);	
		}		
	};


	StrCompare m_cmp; ///< instance of comparator

	LoserTree3Way<StrCompare> m_ltree; ///< instance of losertree

	std::vector<uint32_vector_type*> & m_short_pos_seqs;

	std::vector<typename uint32_vector_type::const_iterator> m_short_pos_its; ///< pos seqs for short substrs

	offset_vector_type * m_long_pos_seq;

	typename offset_vector_type::const_reverse_iterator m_long_pos_rit;

	uint8_vector_type * m_long_aux_seq;

	typename uint8_vector_type::const_reverse_iterator m_long_aux_rit;
	
	std::vector<uint64> & m_blocks_beg_pos; ///< start position for each block

	uint8 m_seq_num;

	uint8 m_pre_loser;

	uint8 m_cur_loser;

	bool m_pre_long_diff; ///< last and current scanned long substrs are different

	bool m_cur_long_diff; ///< current scanned and next to scanned are different

	/// \brief ctor
	///
	SubstrSorter(std::vector<alphabet_vector_type*>& _short_ch_seqs, std::vector<uint8_vector_type*> & _short_len_seqs, std::vector<uint32_vector_type*> & _short_pos_seqs, std::vector<alphabet_vector_type*>& _long_ch_seqs, offset_vector_type* _long_pos_seq, uint8_vector_type* _long_aux_seq, std::vector<uint64> & _blocks_beg_pos) : m_cmp(_short_ch_seqs, _short_len_seqs, _long_ch_seqs), m_ltree(m_cmp), m_short_pos_seqs(_short_pos_seqs), m_long_pos_seq(_long_pos_seq), m_long_aux_seq(_long_aux_seq), m_blocks_beg_pos(_blocks_beg_pos) { 

		for (uint8 i = 0; i < m_short_pos_seqs.size(); ++i) {

			m_short_pos_its.push_back(m_short_pos_seqs[i]->begin());
		}	

		m_long_pos_rit = _long_pos_seq->rbegin();
	
		m_long_aux_rit = _long_aux_seq->rbegin();

		m_seq_num = m_cmp.m_seq_num;

		m_pre_loser = std::numeric_limits<uint8>::max();

		m_pre_long_diff = true;
	}

	/// \brief check if no more substr to be compared
	///
	bool is_empty() const {

		return m_ltree.done();
	}

	/// \brief get current loser
	///
	Substr<alphabet_type, D>& get_cur_loser() {

		m_cur_loser = m_ltree.top();

		return m_cmp.m_head[m_cur_loser];
	}


	/// \brief get position of current loser
	///
	offset_type get_cur_pos() {

		if (m_cur_loser == m_seq_num - 1) { // long

			offset_type pos = *m_long_pos_rit; 

			++m_long_pos_rit;

			return pos;
		}
		else { // short

			offset_type pos = m_blocks_beg_pos[m_cur_loser] + *m_short_pos_its[m_cur_loser];

			++m_short_pos_its[m_cur_loser];

			return pos; 
		}
	}


	/// \brief two successively retrieved substrs are equal
	///
	/// \note only correct if two are at m_tree[1] and m_tree[2]
	bool is_equal() const {

		return m_ltree.top_equal();
	}
	
	/// \brief
	///
	void forward() {

		m_pre_loser = m_cur_loser;

		if (m_pre_loser == m_seq_num - 1) {

			if (m_long_aux_rit != m_long_aux_seq->rend()) {
			
				m_cmp.m_exists[m_pre_loser] = true;

				uint8 aux = *m_long_aux_rit;

				++m_long_aux_rit;
	
				m_cmp.fetchLong(aux & 127);
					
				m_pre_long_diff = m_cur_long_diff;

				m_cur_long_diff = aux & 128;

				if (m_pre_long_diff == false) { // two identical long substrs

					// no need to replay
				}
				else {

					m_ltree.replay();
				}
			}
			else { // empty

				m_cmp.m_exists[m_pre_loser] = false;

				m_ltree.replay();
			}
		}
		else {

			m_cmp.fetchShort(m_pre_loser);

			m_ltree.replay();
		}
	}

	/// \brief retrieve smallest from each seq
	///
	void start() {

		// retrieve short from each seq
		for (uint8 i = 0; i < m_seq_num - 1; ++i) {

			m_cmp.fetchShort(i);
		}

		// retrieve long from long seq
		if (m_long_aux_rit != m_long_aux_seq->rend()) {
			
			m_cmp.m_exists[m_seq_num - 1] = true;

			uint8 aux = *m_long_aux_rit;

			++m_long_aux_rit;

			m_cmp.fetchLong(aux & 127);

			m_cur_long_diff = aux & 128;
		}
		else {

			m_cmp.m_exists[m_seq_num - 1] = false;
		}

		m_ltree.play_initial(m_seq_num);
	}


	/// \brief 
	///
	bool process(offset_vector_type*& _s1, const uint64 _s_size) {

		bool is_unique = true;
		//
		typedef Pair<offset_type, offset_type> pair_type;

		typedef TupleAscCmp1<pair_type> pair_comparator_type;

		typedef typename ExSorter<pair_type, pair_comparator_type>::sorter sorter_type;
	
		sorter_type *sorter_lms = new sorter_type(pair_comparator_type(), MAX_MEM);

		//
		offset_type name = 0; // for the sentinel

		sorter_lms->push(pair_type(_s_size - 1, name)); // push the sentinel

		++name;

		start(); // start comparing

		Substr<alphabet_type, D> pre_str = get_cur_loser();

		offset_type pos = get_cur_pos();

		sorter_lms->push(pair_type(pos, name)); // push current loser

		forward();

		while (!is_empty()) {

			bool is_diff = true;

			Substr<alphabet_type, D> &cur_str = get_cur_loser();
	
			if (is_equal() == true) { // if is_equal == true, then must be equal

				is_diff = false;
			}
			else { // if is_equal == false

				if (m_pre_loser == m_cur_loser && m_pre_loser == m_seq_num - 1) { // two successive are identical long

					if (m_pre_long_diff == false) {

						is_diff = false;
					}
				}
				else if (pre_str.cmp(cur_str) == 0) { // directly compare them 

					is_diff = false;
				}
			}	

			if (is_diff == false) {

				is_unique = false;
			}
			else {
		
				++name;
			}

			offset_type pos = get_cur_pos();

			sorter_lms->push(pair_type(pos, name));

			pre_str.swap(cur_str);

			forward();
		}

		sorter_lms->sort();

#ifdef TEST_DEBUG_SUBSTR_SORTER
		
		std::cerr << "sorted lms num: " << sorter_lms->size() << std::endl;

		for (uint8 i = 0; i < m_seq_num - 1; ++i) {

			assert(m_cmp.m_short_ch_its[i] == m_cmp.m_short_ch_seqs[i]->end());

			assert(m_cmp.m_short_len_its[i] == m_cmp.m_short_len_seqs[i]->end());

			assert(m_cmp.m_long_ch_its[i] == m_cmp.m_long_ch_seqs[i]->end());

			assert(m_short_pos_its[i] == m_short_pos_seqs[i]->end());
		}

		assert(m_long_pos_rit == m_long_pos_seq->rend());

		assert(m_long_aux_rit == m_long_aux_seq->rend());
#endif

		_s1 = new offset_vector_type();

		for (; !sorter_lms->empty(); ++(*sorter_lms)) {

			_s1->push_back((*sorter_lms)->second);
		}

		delete sorter_lms; sorter_lms = nullptr;

		return is_unique;
	}
};




#endif
