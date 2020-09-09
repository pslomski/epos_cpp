#include "catch.hpp"
#include "epos.hpp"
#include "peak_detector_CeBr.hpp"

using namespace std;
using namespace epos;

TEST_CASE("PeakDetector empy data")
{
	VectorFloat time;
	VectorFloat data;
	Float dc(0.0);
	Float dc_delta(0.0);
	int impulse_sample_count(0);
	PeakDetectorCeBr pd(time, data, dc, dc_delta, impulse_sample_count);
	Peak peak;
	REQUIRE(pd.next(peak) == false);
}

void prepare_peaks()
{

}

TEST_CASE("PeakDetector 1")
{
	VectorFloat time;
	VectorFloat data;
	Float dc(0.0);
	Float dc_delta(0.0);
	int impulse_sample_count(0);
	prepare_peaks();
	PeakDetectorCeBr pd(time, data, dc, dc_delta, impulse_sample_count);
	Peak peak;
	REQUIRE(pd.next(peak) == false);
}
