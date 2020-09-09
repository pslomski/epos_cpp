#include <stdexcept>
#include "histogram.hpp"

/*
const int NSAMPLES = whatever;
double samples[NSAMPLES] = { 1.0, 3.93, 1e30, ... }; // your data set
const int NBUCKETS = 10; // or whatever
int counts[NBUCKETS] = { 0 };
for (int i = 0; i != NSAMPLES; ++i) {
	counts[TRANSFER(samples[i])]++;
}
*/
using namespace epos;

Histogram::Histogram() :
	domain_min(0),
	domain_max(0.0),
	bins_count(0),
	bin_width(0.0),
	a_coeff(0.0),
	b_coeff(0.0),
	bins(0),
	hist_max(0),
	hist_max_idx(0)
{}

void Histogram::init(size_t count, Float min, Float max)
{
	if (count < 1)
		throw std::invalid_argument("Histogram::init: parameter 'count' should be greater than zero");
	if (min >= max)
		throw std::invalid_argument("Histogram::init: parameter 'min' should be lower than 'max'");
	hist_max = 0;
	hist_max_idx = 0;
	bins_count = count;
	domain_min = min;
	domain_max = max;
	bin_width = (max - min) / count;
	a_coeff = bins_count / (domain_max - domain_min);
	b_coeff = -domain_min * a_coeff;
	bins.resize(count);
	clear();
}

void Histogram::clear()
{
	hist_max = 0;
	hist_max_idx = 0;
	for (auto& val : bins)
		val = 0;
}

bool Histogram::add(Float val)
{
	Float d = a_coeff * val + b_coeff;
	int idx = static_cast<int>(floor(d));
	if ((idx >= 0) && (idx < bins_count)) {
		bins[idx]++;
		if (bins[idx] > hist_max) {
			hist_max = bins[idx];
			hist_max_idx = idx;
		}
		return true;
	}
	return false;
}

double Histogram::get_interval(size_t idx) const
{
	if (idx < bins_count) {
		return domain_min + (idx + 0.5) * bin_width;
	}
	else {
		throw std::out_of_range("Histogram::getInterval - Index out of range");
	}
}

void Histogram::saveToFile(StringRef filepath) const
{
	std::ofstream os(filepath);
	if (os.is_open()) {
		for (size_t i = 0; i < bins.size(); ++i) {
			os << get_interval(i) << "\t" << bins[i] << std::endl;
		}
	}
}
