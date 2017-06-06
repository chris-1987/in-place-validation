#ifndef SACHECKERLIGHTWEIGHT_H
#define SACHECKERLIGHTWEIGHT_H


template<typename charT>
class SaCheckerLightWeight{
private:
	const int32_t P = 2047483525;
	const int32_t R = 1732938289;
	
	int32_t fp1;
	int32_t fp2;

public:

	SaCheckerLightWeight(): fp1(0), fp2(0){}

	void computeFP1(charT t) {
		fp1 = ((int64_t)fp1 * R + t + 1) % P;
	}

	void computeFP2(charT t) {
		fp2 = ((int64_t)fp2 * R + t + 1) % P;
	}

	bool validate() {
		return (fp1 == fp2);
	}	
	
	void print() {
		std::cerr << "fp1: " << fp1 << " fp2: " << fp2 << std::endl;	
	}
};


#endif
