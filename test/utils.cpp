#include "utils.hpp"
#include <epos.hpp>

using namespace epos;

void avevar(VectorFloat& data, Float& ave, Float& var)
{
	//Given array data[0..n - 1], returns its mean as aveand its variance as var.
	Float s, ep;
	size_t j, n = data.size();
	ave = 0.0;
	for (j = 0; j < n; j++) ave += data[j];
	ave /= n;
	var = ep = 0.0;
	for (j = 0; j < n; j++) {
		s = data[j] - ave;
		ep += s;
		var += s * s;
	}
	var = (var - ep * ep / n) / (n - 1); // Corrected two - pass formula(14.1.8).
}