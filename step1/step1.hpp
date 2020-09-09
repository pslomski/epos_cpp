#ifndef EPOS_STEP1_H_INCLUDED
#define EPOS_STEP1_H_INCLUDED

#include "epos.hpp"
#include "fit.hpp"

class Step1
{
	friend class UnitTestStep1Processor;
private:
	size_t fCountClipped;
	size_t fCountFailed;
	size_t fCountOk;
	std::string fRootPath;
	epos::VectorFloat fX;
	bool on_segment(epos::Segment& seg);
	bool check_background(epos::Channel& channel, epos::Float beginDCcoeff, epos::Float endDCcoeff);
	bool process_ch1(epos::Channel& channel);
	bool process_ch2(epos::Channel& channel);
	bool process_ch3(epos::Channel& channel);
	bool get_fit_params(epos::Channel& channel, int t0idx, epos::Float Ampl, epos::FitParams& fitPar, epos::Float s, epos::Float tau, int DCcount, int impulseLen);
	bool fit(epos::FitParams& fitPar, epos::Channel::Result& chRes, epos::Float maxChisq, epos::Float maxSigmaS, epos::Float maxSigmaTau);
	void write_results_to_file(epos::Segment& seg);
	static bool fit(epos::VectorFloat& x, epos::VectorFloat& y, epos::VectorFloat& a, epos::VectorFloat & sig_a);
public:
	Step1();
	~Step1();
	void init(epos::StringRef RootPath);
	void process_segment(epos::Segment& seg);
	static bool processCeBr(int buf[1600], double DC, double zeroTime, epos::VectorFloat& a);
	epos::VectorFloat processSiPM(int buf[1600], double DC, double zeroTime);
};

#endif // EPOS_STEP1_H_INCLUDED
