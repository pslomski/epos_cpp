#ifndef LIB_EPOS_PEAK_DETECTOR_HPP_INCLUDED
#define LIB_EPOS_PEAK_DETECTOR_HPP_INCLUDED

#include <epos.hpp>

namespace epos
{

	class PeakDetectorCeBr : public PeakDetector {
	private:
		VectorFloat& time;
		VectorFloat& data;
		VectorFloat vec_a; // nachylenie
		Float dc;
		Float dc_delta;
		int impulse_sample_count;
		int curr;
		int stop;
		int insidePeak;
		Float aceil;
		Float yceil;
		Float A;
		Float Amin;
		bool checkLen(int start, int stop);
	public:
		PeakDetectorCeBr(VectorFloat& time_, VectorFloat& data_, Float dc_, Float dc_delta_, int impulse_sample_count_);
		bool next(Peak& peak);
	};

}

#endif // LIB_EPOS_PEAK_DETECTOR_HPP_INCLUDED
