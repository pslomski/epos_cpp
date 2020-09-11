#include <epos.hpp>
#include <algorithm>
#include "peak_detector_CeBr.hpp"

namespace epos
{
	PeakDetectorCeBr::PeakDetectorCeBr(VectorFloat& time_, VectorFloat& data_, Float dc_, Float dc_delta_, int impulse_sample_count_) :
		time(time_),
		data(data_),
		dc(dc_),
		dc_delta(dc_delta_),
		impulse_sample_count(impulse_sample_count_)
	{
		vec_a.resize(time.size(), 0.0);
		curr = 2;
		stop = static_cast<int>(time.size()) - 2;
		insidePeak = 0;
		aceil = -1.0 * dc_delta / epos::time_resolution; // gorna granica nachylenia
		yceil = 1.0 * (dc - dc_delta); // wartosc graniczna ponizej poziomu tla
		A = 0.0;
		Amin = 0.1;
	}

	bool PeakDetectorCeBr::checkLen(int start, int stop)
	{
		int length = 20;
		int curr = start;
		while ((curr < stop) && (curr - start < length) && (data[curr] < yceil)) {
			curr += 1;
		}
		return curr - start >= length;
	}

	bool PeakDetectorCeBr::next(Peak& peak)
	{
		while (curr < stop) {
			if (insidePeak > 0) {
				insidePeak--;
			}
			else {
				Float a = nachylenie5(data, curr);
				vec_a[curr] = a;
				// rozpatrujemy nachylenia mniejsze(ujemne!) od aceil
				if ((a < aceil) && (data[curr] < yceil)) {
					// punkt odniesienia dla oszacowania amplitudy
					Float y0 = data[curr - 1];
					Float amin = 1000;
					int aimin = -1;
					Float ymin = y0;
					int stop2 = std::min(curr + 10, stop);
					// przegladamy 10 kolejnych punktow w poszukiwaniu maksymalnego nachylenia i maksimum piku
					while (curr < stop2) {
						Float a = nachylenie5(data, curr);
						if (a < amin) {
							amin = a;
							aimin = curr;
						}
						if (data[curr] < ymin) {
							ymin = data[curr];
						}
						curr += 1;
					}

					A = abs(ymin - y0); // amplituda piku
					if (checkLen(curr, stop) && (aimin > -1) && (A > Amin)) {
						// wpasowujemy prosta w zbocze narastania impulsu
						Float a(1), b{ 0 };
						//a, b = fit.LinearFit(self.t[aimin - 2:aimin + 2], self.y[aimin - 2:aimin + 2]);
						// start impulsu wyznacza przeciecie prostej z tlem
						peak.t0 = (dc - b) / a;
						peak.t0idx = Channel::getIdx(peak.t0);
						peak.y0 = a * peak.t0 + b;
						insidePeak = impulse_sample_count; // = 50ns
						curr++; // przed wyjsciem przestawiamy indeks na kolejny punkt dla nastepnej iteracji
						return true;
					}
				}
			}
			curr++;
		}
		return false;
	}
} // namespace epos