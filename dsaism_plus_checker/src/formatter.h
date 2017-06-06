#ifndef MY_FORMAT_H
#define MY_FORMAT_H

#include "common.h"
#include "io.h"

//#define TEST_DEBUG

/// \brief format the input string
///
/// delete alpha_max, plus each remaining by 1
template<typename alphabet_type>
class Formatter{

private:

	typedef typename ExVector<alphabet_type>::vector alphabet_vector_type;

private:

	std::string m_s_fname;

	std::string m_target_s_fname;

public:
	
	/// \brief ctor
	///
	Formatter(const std::string & _s_fname, const std::string & _target_s_fname) : 
			m_s_fname(_s_fname), m_target_s_fname(_target_s_fname) {}

	/// \brief 
	///
	void run() {

		stxxl::syscall_file *s_file = new stxxl::syscall_file(m_s_fname, stxxl::syscall_file::RDWR | stxxl::syscall_file::DIRECT);

		alphabet_vector_type *s = new alphabet_vector_type(s_file);

		typename alphabet_vector_type::const_iterator it = s->begin();

		stxxl::syscall_file *s_target_file = new 
			stxxl::syscall_file(m_target_s_fname, stxxl::syscall_file::CREAT | stxxl::syscall_file::RDWR | stxxl::syscall_file::DIRECT);

		alphabet_vector_type *s_target = new alphabet_vector_type(s_target_file);

		for (uint64 i = 0; it != s->end(); ++it, ++i) {

			if (*it != std::numeric_limits<alphabet_type>::max()) {

				s_target->push_back(*it + 1); // increment by 1
			}
			else {

				// get rid of it
			}
		}

		delete s; s = nullptr;

		delete s_file; s_file = nullptr;

		delete s_target; s_target = nullptr;

		delete s_target_file; s_target_file = nullptr;
	}
};


#endif
