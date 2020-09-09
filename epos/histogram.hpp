#ifndef LIB_EPOS_HISTOGRAM_HPP_INCLUDED
#define LIB_EPOS_HISTOGRAM_HPP_INCLUDED

#include "epos.hpp"

namespace epos
{
	class Histogram {
	private:
		double domain_min; // minimum przedzialu
		double domain_max; // maksimum przedzialu
		size_t bins_count; // liczba kanalow
		double bin_width; // szerokosc kanalu
		 // wspolczynniki a i b prostej a*val+b mapujacej wartosc na indeks tablicy
		double a_coeff;
		double b_coeff;
		VectorInt bins; // wartosci histogramu
		size_t hist_max; // max hist value
		size_t hist_max_idx; // index of max hist value
	public:
		Histogram();
		void init(size_t count, Float min, Float max);
		void clear();
		bool add(Float val);
		const VectorInt& getBins() const { return bins; }
		double get_interval(size_t idx) const;
		void saveToFile(StringRef filepath) const;
		size_t getMaxHistIdx() { return hist_max_idx; }
	};
}

#endif // LIB_EPOS_HISTOGRAM_HPP_INCLUDED