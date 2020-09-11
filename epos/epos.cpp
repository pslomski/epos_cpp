#define  _CRT_SECURE_NO_WARNINGS
#include "epos.hpp"
#include <sys/types.h>
#include <sys/stat.h>
#include <limits>
#include <algorithm>
#include <math.h>
#include <numeric>

#include "histogram.hpp"

using namespace epos;
using namespace std;

// Mean of [begin, end) range
Float epos::arithmetic_mean(const VectorFloatIter begin, const VectorFloatIter end)
{
	auto count = end - begin;
	if (count == 0)
		return 0.0;

	return std::accumulate(begin, end, 0.0) / count;
}

// Standard deviation of [begin, end) range
Float epos::variance(const VectorFloatIter begin, const VectorFloatIter end)
{
	auto count = end - begin;
	if (count < 2)
		return 0.0;

	Float mean = std::accumulate(begin, end, 0.0) / count;
	Float accum = 0.0;
	std::for_each(begin, end, [&](const Float d) {
		Float delta = d - mean;
		accum += delta * delta;
		});
	return accum / (count -1);
}

void epos::avevar(VectorFloatIter begin, VectorFloatIter end, Float& ave, Float& var)
{
	ave = 0.0;
	var = 0.0;
	auto count = end - begin;
	if (count < 2)
		return;

	ave = std::accumulate(begin, end, 0.0) / count;
	Float accum = 0.0;
	std::for_each(begin, end, [&](const Float d) {
		Float delta = d - ave;
		accum += delta * delta;
		});
	var = accum / (count - 1);
}


Float epos::nachylenie5(VectorFloat& data, int i0)
{
	constexpr Float time_res_inv = 1.0 / epos::time_resolution;
	Float XY = 2.0 * (data[i0 + 2] - data[i0 - 2]) + data[i0 + 1] - data[i0 - 1];
	return time_res_inv * XY / 10.0;
}

void SeriesStatistics::compute_from_range(VectorFloatIter begin, VectorFloatIter end)
{
	arithmetic_mean = 0.0;
	variance = 0.0;
	standard_deviation = 0.0;
	auto count = end - begin;
	if (count > 0) {
		arithmetic_mean = std::accumulate(begin, end, 0.0) / count;
	}
	if (count > 1) {
		Float accum = 0.0;
		std::for_each(begin, end, [&](const double d) {
			Float delta = d - arithmetic_mean;
			accum += delta * delta;
			});
		variance = accum / (count - 1); // with Bessel correction
		standard_deviation = sqrt(variance);
	}
}

// Covariance
template <class Iter> typename Iter::value_type covariance(const Iter& x, const Iter& y)
{
	double sum_x = std::accumulate(std::begin(x), std::end(x), 0.0);
	double sum_y = std::accumulate(std::begin(y), std::end(y), 0.0);

	double mx = sum_x / x.size();
	double my = sum_y / y.size();

	double accum = 0.0;

	for (auto i = 0; i < x.size(); i++)	{
		accum += (x.at(i) - mx) * (y.at(i) - my);
	}

	return accum / (x.size() - 1);
}

long long getFileSize(ifstream& file)
{
	auto current = file.tellg();
	file.seekg(0, ios_base::end);
	auto endpos = file.tellg();
	file.seekg(current, ios_base::beg);
	return endpos;
}

//////////////////////////////////////////////////
// Channel class
//

// Local function to initialize the static Channel::time
static VectorFloat _getTime()
{
	VectorFloat time;
	for (int i = 0; i < epos::NbrSAMP; ++i) {
		time.push_back(Channel::getT(i));
	}
	return time;
}

VectorFloat Channel::time = _getTime();

void Channel::init(size_t dataSize)
{
	data.clear();
	data.reserve(dataSize);
}

Channel::Stats& epos::Channel::computeStats()
{
	stats.compute_statistics(data);
	return stats;
}

// liczy srednia dla pierwszych sample_count wartosci
double Channel::getMean(size_t sampleCount) const
{
	sampleCount = std::min(sampleCount, data.size());
	if (sampleCount > 0) {
		double mean = 0.0;
		for (int i = 0; i < sampleCount; ++i)
			mean += data[i];
		return mean / (Float)sampleCount;
	}
	else
		return 0.0;
}

// Standard deviation of background level.
// W dwa odchylenia standardowe wpada 0,95449989 wartosci.
// return:
// >= 0 - value of Standard deviation
// -1 - error
Float compute_stddev_of_dc(const VectorFloat& data, Float dc)
{
	constexpr Float dc_delta = 0.022; // stala okreslajaca zmiennosc tla oszacowana z wykresu
	Float sigma = 0.0; // odchylenie standardowe
	int count = 0;
	for (Float y : data) {
		if (isInRange(y, dc - dc_delta, dc + dc_delta)) {
			Float delta = y - dc;
			sigma += delta * delta;
			count += 1;
		}
	}
	if (count > 0)
		return sqrt(sigma / count);
	else
		return 0.0;
}

void epos::Channel::Stats::compute_statistics(const VectorFloat& data)
{
	Histogram hist;
	hist.init(1000, -1.0, 1.0);
	min = numeric_limits<Float>::max();
	max = -min;
	idxmin = -1;
	idxmax = -1;
	int i = 0;
	for (auto const& val : data) {
		if (val < min) {
			min = val;
			idxmin = i;
		}
		if (val > max) {
			max = val;
			idxmax = i;
		}
		i++;
		hist.add(val);
	}
	DC = hist.get_interval(hist.getMaxHistIdx());
	DCstddev = compute_stddev_of_dc(data, DC);
}

/*void Channel::addConst(Float c)
{
	for (auto& val : data)
		val += c;
}*/

void Channel::saveToFile(StringRef filePath) const
{
	ofstream os(filePath);
	if (os.is_open()) {
		for (size_t i = 0; i < data.size(); i++) {
			os << getT(i) << "\t" << data[i] << endl;
		}
	}
	else {
		cout << "Brak dostepu do pliku " << filePath << endl;
	}
}

//////////////////////////////////////////////////
// Chanell class

BinChannel::BinChannel()
{
	fFileSize = 0;
}

BinChannel::~BinChannel()
{
	close();
}

bool BinChannel::open(StringRef filePath)
{
	fFile.open(filePath, ios::in | ios::binary);
	if (fFile.is_open()) {
		fFileSize = getFileSize(fFile);
		return fFileSize != -1;
	}
	else {
		cout << "error opening " << filePath;
		fFileSize = 0;
		return false;
	}
}

void BinChannel::close()
{
	fFile.close();
}

bool BinChannel::read(size_t segNr, Channel& seg)
{
	if (segNr >= 0 && segNr < getSegCount()) {
		SegBuf buf(epos::NbrSAMP); // bufor odczytu segmentu
		int bufByteCount = BinChannel::bufByteCount; // liczba bajtow do odczytu
		fFile.seekg(getSegStartPos(segNr));
		fFile.read((char*)buf.data(), bufByteCount);
		if (fFile.gcount() == bufByteCount) {
			seg.init(epos::NbrSAMP);
			constexpr Float scale = 1.0 / numeric_limits<short int>::max();
			for (Sample sample : buf)
				seg.add(scale * sample);

			if (seg.data.size() != epos::NbrSAMP) {
				cerr << "incorrect size for segment: " << segNr << endl;
				return false;
			}
			else
				return true;
		}
		else {
			cerr << "erorr reading segment number: " << segNr << endl;
			return false;
		}
	}
	else {
		cerr << "bad segment number: " << segNr << endl;
		return false;
	}
}

void Segment::saveToFile(StringRef filePath) const
{
	ofstream os(filePath);
	if (os.is_open()) {
		for (size_t i = 0; i < ch1.data.size(); i++) {
			os << Channel::getT(i) << "\t" <<
				ch1.data[i] << "\t" <<
				ch2.data[i] << "\t" <<
				ch3.data[i] << endl;
		}
	}
	else {
		cout << "Brak dostepu do pliku " << filePath << endl;
	}
}

bool DataSetBin::isValid() const
{
	bool res = (fCh1.getSegCount() == fCh2.getSegCount())
		&& (fCh2.getSegCount() == fCh3.getSegCount());
	if (res == false) {
		cerr << "data count mismatch in channels" << endl;
	}
	return res;
}

bool DataSetBin::open(StringRef rootPath, size_t segMin, size_t segMax)
{
	fCurrSeg = segMin;
	bool Res = fCh1.open(rootPath + "ch1.bin") &&
		fCh2.open(rootPath + "ch2.bin") &&
		fCh3.open(rootPath + "ch3.bin") &&
		isValid();
	fSegMax = std::min(segMax, fCh1.getSegCount());
	return Res;
}

bool DataSetBin::read(Segment& data)
{
	if (fCurrSeg <= fSegMax)
		return getData(fCurrSeg++, data);
	else
		return false;
}

bool DataSetBin::getData(size_t segNr, Segment& data)
{
	data.segNr = segNr;
	return fCh1.read(segNr, data.ch1) &&
		fCh2.read(segNr, data.ch2) &&
		fCh3.read(segNr, data.ch3);
}

DataSetTxt::DataSetTxt(Float scalingFactor) :
	fScalingFactor(scalingFactor),
	fRootPath(""),
	fStop(-1),
	fCurrent(-1),
	fSegNr(NbrSEG)
{
}

bool DataSetTxt::open(StringRef rootPath, int fileMin, int fileMax)
{
	fRootPath = rootPath;
	fCurrent = fileMin - 1;
	fStop = fileMax;
	fSegNr = NbrSEG;
	return true;
}

bool DataSetTxt::openFile(int fileNr)
{
	bool result = false;
	stringstream fileName;
	fileName << fRootPath << fileNr * NbrSEG << "-" << (fileNr + 1) * NbrSEG << ".txt";
	fDataFile.close();
	fDataFile.open(fileName.str());
	result = fDataFile.is_open();
	if (result) {
		// Skip 4 header lines
		string line;
		for (size_t i = 0; i < 4; i++)
			getline(fDataFile, line);
	}
	else
		cerr << "error opening " << fileName.str() << endl;
	return result;
}

bool isClipped(Sample val)
{
	constexpr short int clipMin = numeric_limits<short int>::min();
	constexpr short int clipMax = numeric_limits<short int>::max();
	return (val <= clipMin) || (val >= clipMax);
}

bool DataSetTxt::read(Segment& data)
{
	if (fSegNr >= NbrSEG) {
		// Start new file
		fSegNr = 0;
		fCurrent++;
		if ((fCurrent > fStop) || (openFile(fCurrent) == false))
			return false;
	}
	data.clipped = false;
	data.segNr = fSegNr;
	data.fileNr = fCurrent;
	data.ch1.init(NbrSAMP);
	data.ch2.init(NbrSAMP);
	data.ch3.init(NbrSAMP);
	for (int i = 0; i < NbrSAMP; i++) {
		string line;
		getline(fDataFile, line);
		char* p = strtok((char*)line.c_str(), "\t");

		Sample ch1 = atoi(p);
		p = strtok(NULL, "\t");
		Sample ch2 = atoi(p);
		p = strtok(NULL, "\t");
		Sample ch3 = atoi(p);

		data.clipped = data.clipped || isClipped(ch1) || isClipped(ch2) || isClipped(ch3);
		data.ch1.add(ch1 * fScalingFactor);
		data.ch2.add(ch2 * fScalingFactor);
		data.ch3.add(ch3 * fScalingFactor);
	}
	fSegNr++;
	return data.ch1.data.size() == NbrSAMP;
}

void Highpass::Init(double RC, double dt, double startVal)
{
	fAlpha = RC / (RC + dt);
	fXprev = startVal;
	fYprev = startVal;
}

double Highpass::Get(double x)
{
	//double y = fAlpha * fYprev + fAlpha * (x - fXprev);
	double y = fAlpha * (fYprev + x - fXprev);
	fXprev = x;
	fYprev = y;
	return y;
}

void Highpass::FilterBuf(double RC, double dt, const VectorFloat& in, VectorFloat& out)
{

	double alfa = RC / (RC + dt);
	out[0] = in[0];
	for (size_t i = 1; i < in.size(); i++) {
		out[i] = alfa * (out[i - 1] + in[i] - in[i - 1]);
	}
}

void Lowpass::Init(double RC, double dt, double startVal)
{
	fAlpha = RC / (RC + dt);
	fXprev = startVal;
	fYprev = fAlpha * startVal;
}

double Lowpass::Get(double x)
{
	double y = fYprev + fAlpha * (x - fYprev);
	fXprev = x;
	fYprev = y;
	return y;
}

void Lowpass::FilterBuf(double RC, double dt, const VectorFloat& in, VectorFloat& out)
{
	double alfa = dt / (RC + dt);
	out[0] = alfa * in[0];
	for (size_t i = 1; i < in.size(); i++) {
		out[i] = out[i - 1] + alfa * (in[i] - out[i - 1]);
	}
}

