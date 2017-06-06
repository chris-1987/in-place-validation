#ifndef MY_UTILITY_H
#define MY_UTILITY_H

#include "common.h"


static const uint8 MASK[] = { 0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01 };

/// \brief data wrapper for a bit array in RAM
///
class BitWrapper{

private:

	char * data;

public:

	BitWrapper(char * _data = nullptr) : data(_data) {}

	bool get(size_t _idx) {

		return  (data[_idx / 8] & MASK[_idx % 8]) ? 1 : 0;
	}

	void set(size_t _idx, bool _val) {

		data[_idx / 8] = _val ? (MASK[_idx % 8] | data[_idx / 8]) : ((~MASK[_idx % 8]) & data[_idx / 8]);
	}
};



static const uint8 MASK2[] = { 0xC0, 0x30, 0x0C, 0x03 };

static const uint8 SHIFT2[] = { 6,4,2,0 };

/// \brief data warpper for a 2-bit array in RAM
///
class Bit2Wrapper {
private:

	char * data;

public:

	Bit2Wrapper(char* _data = nullptr) : data(_data) {}

	uint8 get(size_t _idx) {

		return (data[_idx / 4] >> SHIFT2[_idx % 4]) & MASK2[3];
	}

	void set(size_t _idx, uint8 _val) {

		data[_idx / 4] = (_val << SHIFT2[_idx % 4]) | (data[_idx / 4] & ~MASK2[_idx % 4]);
	}
};

#endif //MY_UTILITY_H
