#include <algorithm>
#include <math.h>
#include "nr.h"
#include "nrutil.h"
#include "step1.hpp"
#include <impulse_fun.hpp>
#include <peak_detector_CeBr.hpp>

using namespace std;
using namespace epos;

Step1::Step1()
{
	init("");
}

Step1::~Step1()
{
}

Peaks get_all_peaks(PeakDetector& peak_detector)
{
	constexpr Float Amin = 0.1;
	Peaks peaks;
	Peak peak;
	while (peak_detector.next(peak)) {
		if (peak.A > Amin) {
			peaks.push_back(peak);
		}
	}
	return peaks;
}

void Step1::init(StringRef RootPath)
{
	fCountClipped = 0;
	fCountFailed = 0;
	fCountOk = 0;
	fRootPath = RootPath;
	for (size_t i = 0; i < epos::NbrSAMP; i++) {
		fX.push_back(Channel::getT(i));
	}
}

void Step1::process_segment(epos::Segment& seg)
{
	if (seg.clipped) {
		fCountClipped++;
	}
	else {
		if (on_segment(seg)) {
			fCountOk++;
		}
		else
			fCountFailed++;
	}
}

bool Step1::on_segment(epos::Segment& seg)
{
	bool res = false;
	if (check_background(seg.ch1, 2.5, 2.5)
		&& check_background(seg.ch2, 2.5, 2.5)
		&& check_background(seg.ch3, 10.0, 10.0)
		&& process_ch1(seg.ch1)
		&& process_ch2(seg.ch2)
		&& process_ch3(seg.ch3)) {
		write_results_to_file(seg);
		res = true;
	}
	return res;
}

// Checks if the signal starts at the background level
bool Step1::check_background(epos::Channel& channel, epos::Float beginDCcoeff, epos::Float endDCcoeff)
{
	channel.stats.compute_statistics(channel.data);

	Float DCDelta = beginDCcoeff * channel.stats.DCstddev;
	//*
	// czy srednia wartosc 15 poczatkowych sampli miesci sie w granicach poziomu tla
	Float yStart = epos::arithmetic_mean(channel.data.begin(), std::next(channel.data.begin(), 15));
	if (!epos::isInRange(yStart, channel.stats.DC - DCDelta, channel.stats.DC + DCDelta)) {
		//self._whyFailed = channel.ID + ': poczatek segmentu poza tlem'
		return false;
	}
	// czy srednia wartosc 15 koncowych sampli miesci sie w granicach poziomu tla
	DCDelta = endDCcoeff * channel.stats.DCstddev;
	Float yStop = epos::arithmetic_mean(std::next(channel.data.end(), -15), channel.data.end());
	if (!epos::isInRange(yStop, channel.stats.DC - DCDelta, channel.stats.DC + DCDelta)) {
		//self._whyFailed = channel.ID + ": koniec segmentu poza tlem";
		return false;
	}
	//*/
	return true;
}

bool Step1::get_fit_params(Channel& channel, int t0idx, Float Ampl, FitParams& fitPar, Float s, Float tau, int DCcount, int impulseLen)
{
	return false;
}

bool Step1::fit(FitParams& fitPar, Channel::Result& chRes, Float maxChisq, Float maxSigmaS, Float maxSigmaTau)
{
	return false;
}

bool Step1::process_ch1(Channel& channel)
{
	bool res = false;
	Float DCdelta = 2 * channel.stats.DCstddev;
	int impulseSampleCount = 70 * 2; // artitrary chosen impulse length
	FitParams fitPar;
	PeakDetectorCeBr peakDetector = PeakDetectorCeBr(channel.time, channel.data, channel.stats.DC, DCdelta, impulseSampleCount);
	Peaks peaks = get_all_peaks(peakDetector);
	if (peaks.size() == 1) {
		Peak& peak = peaks[0];
		// s(0.9) i tau(20.0) dobrze okreslaja ksztalt impulsu CeBr
		Float s = 0.9;
		Float tau = 20.0;
		int DCcount = 20;
		Float maxChisq = 10.0;
		Float maxSigmaS = 2.0;
		Float maxSigmaTau = 4.0;
		if (get_fit_params(channel, peak.t0idx, peak.A, fitPar, s, tau, DCcount, impulseSampleCount) &&
			fit(fitPar, channel.res, maxChisq, maxSigmaS, maxSigmaTau))
		{
			// # dt - s; (s, tau, dt, J, C)
			channel.res.t0 = channel.res.a[2] - channel.res.a[0];
			res = true;
		}

	}
	else {
		//self._whyFailed = channel.ID + ': liczba pikow: ' + str(len(peaks))
	}
	return res;
}

bool Step1::process_ch2(epos::Channel& channel)
{

	return false;
}

bool Step1::process_ch3(epos::Channel& channel)
{
	return false;
}

void Step1::write_results_to_file(epos::Segment& seg)
{
}

bool Step1::processCeBr(int buf[1600], double DC, double zeroTime, epos::VectorFloat& a)
{
	const size_t ParamCount = 5;
	bool res = true;
	const double SigmaSError = 1.0;//dobrac doswiadczalnie
	VectorFloat x, y;
	size_t zeroIndex = Channel::getIdx(zeroTime);
	size_t startIndex = std::max(size_t(0), zeroIndex - 40);
	//pierwornie bylo 'startIndex + 200', ale dla PorousSilica trzeba 300, zeby mrqmin nie wywalalo bledu alokacji. Zbyt mala liczba sampli od max. w prawo powoduje zrzuty pamieci 
	size_t stopIndex = std::min(startIndex + 300, epos::NbrSAMP);

	//z eksperymentow w Pythonie
	//[s=0.88259  tau=20.3255
	//a[1] - s
	//a[2] - tau
	//a[3] - dt przesuniecie
	//a[4] - J amplituda
	//a[5] - C stala
	a.resize(ParamCount + 1);
	a[1] = 0.88;//2.5 charakterystyczne dla CeBr
	a[2] = 20.3;// wartosc 24 jest charakterystyczne dla CeBr a 100 dla SiPM i dobre ustawienie na poczatku jest bardzo wazne!
	a[3] = zeroTime;
	a[4] = 1;// wstepnie 1 dla oszacowania J
	a[5] = 0;

	double yfun_min = numeric_limits<double>::max();
	double ydat_min = yfun_min;
	array<double, ParamCount + 1> da;
	for (size_t i = startIndex; i < stopIndex; i++)
	{
		double yfun;
		double xdat = Channel::getT(i);
		double ydat = (buf[i] - DC) / numeric_limits<short int>::max();
		x.push_back(xdat);
		y.push_back(ydat);
		//dla oszacowania J szukamy minimow danych y i f(x)
		fnImpulse(xdat, a.data(), &yfun, da.data(), ParamCount);
		if (yfun < yfun_min)
			yfun_min = yfun;
		if (ydat < ydat_min)
			ydat_min = ydat;
	}
	a[4] = 0.9 * ydat_min / yfun_min;

	VectorFloat sig_a(a.size(), 0.0);
	res = fit(x, y, a, sig_a);
	if (res)
	{
		if (sig_a[1] > SigmaSError) {
			res = false; // zbyt duzy blad s
		}
	}
	return res;
}

VectorFloat Step1::processSiPM(int buf[1600], double DC, double zeroTime)
{
	VectorFloat x, y, a(5 + 1);
	size_t zeroIndex = Channel::getIdx(zeroTime);
	size_t startIndex = std::max(size_t(0), zeroIndex - 40);
	size_t stopIndex = std::min(startIndex + 300, epos::NbrSAMP);

	//a[1] - s
	//a[2] - tau
	//a[3] - dt przesuniecie
	//a[4] - J amplituda
	//a[5] - C stala
	a[1] = 1.5;
	a[2] = 100.0;//24.0;. 24 jest charakterystyczne dla CeBr a 100 dla SiPM i dobre ustawienie na poczatku jest bardzo wazne!
	a[3] = zeroTime;
	a[4] = 300.0;//100.0
	a[5] = 0;

	for (size_t i = startIndex; i < stopIndex; i++) {
		x.push_back(Channel::getT(i));
		y.push_back((buf[i] - DC) / numeric_limits<short int>::max());
	}

	VectorFloat sig_a(a.size(), 0.0);
	fit(x, y, a, sig_a);

	return a;
}

double ComputeSigma(double chisq, double covar, int ndata, int ma)
{
	return sqrt(chisq * covar / (double)(ndata - ma));
}

// x, y - set of data points x[1..ndata], y[1..ndata]
// funcs - nonlinear function dependent on coefficients a[1..ma]
// chisq - dokladnosc dopasowania
bool Step1::fit(VectorFloat& x, VectorFloat& y, VectorFloat& a, VectorFloat& sig_a)
{
	bool res = true;
	const size_t MaxSteps = 100; // maksymalna liczba krokow fitowania

	VectorFloat sig; // Odchylenia standardowe danych y. Przyjmujemy = 1
	int ndata;
	VectorInt ia;
	int ma;// coefficient count

	//standardowe odchyl. kazdego z data (przyjmuje domyslnie 1)
	for (size_t i = 0; i < x.size(); i++) {
		sig.push_back(1.0);
	}

	// wartosc=!=0 oznacza ze parametr ma byc dopasowywany. Wartosc 0 oznacza, ze ma byc zafiksowany na wart. pocz.
	for (size_t i = 0; i < a.size(); i++)
		ia.push_back(1);

	ma = static_cast<int>(a.size() - 1);
	ndata = static_cast<int>(x.size() - 1);
	double alamda = -1; // set alamda<0 for initialization
	double** covar = dmatrix(1, ma, 1, ma);
	double** alpha = dmatrix(1, ma, 1, ma);

	int step = 0;
	double lambda = -1.0; // set alamda<0 for initialization
	double chisq = 0.0;

	mrqmin(x.data(), y.data(), sig.data(), ndata, a.data(), ia.data(), ma, covar, alpha, &chisq, fnImpulse, &lambda);

	std::cout << "#step " << step << ", chisq=" << chisq << ", lambda=" << lambda << std::endl;
	double ochisq = chisq + 2.0;

	/* Iterating with mrqmin until relative change in chisq is less than 0.001 */
	while ((fabs((ochisq - chisq) / ochisq) > 0.00001) || (lambda >= 0.1))
	{
		step++;
		ochisq = chisq;
		mrqmin(x.data(), y.data(), sig.data(), ndata, a.data(), ia.data(), ma, covar, alpha, &chisq, fnImpulse, &lambda);

		std::cout << "#step " << step << ", chisq=" << chisq << ", lambda=" << lambda << std::endl;
		if (step > MaxSteps) {
			res = false; // brak zbieznosci w zadanej liczbie krokow
			break;
		}
		if (lambda > 2) {
			res = false; // rozbieznosc
			break;
		}
	}

	if (res)
	{
		/* Here we run mrqmin one more time with alambda = 0.0 to get covariant matrix, in case we want the errors */
		lambda = 0.0;//last call with 0
		mrqmin(x.data(), y.data(), sig.data(), ndata, a.data(), ia.data(), ma, covar, alpha, &chisq, fnImpulse, &lambda);

		a[0] = chisq;//wyrzucane jest na zewnatrz Chi bedace surowa suma dla wszystkich punktow

		sig_a.resize(a.size());
		sig_a[0] = chisq;
		for (size_t i = 1; i < a.size(); i++)
		{
			//obliczanie bledu dopasowania parametrow a
			sig_a[i] = ComputeSigma(chisq, covar[i][i], ndata, ma);
			cout << "#a[" << i << "]=" << a[i] << ", sig=" << sig_a[i] << endl;
		}
	}
	free_dmatrix(covar, 1, ma, 1, ma);
	free_dmatrix(alpha, 1, ma, 1, ma);

	return res;
}

