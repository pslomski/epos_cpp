// epos_ftest.cpp : Ten plik zawiera funkcję „main”. W nim rozpoczyna się i kończy wykonywanie programu.
//

#include <iostream>
#include <limits>
#include "../epos/nr.h"
#include "../epos/epos.h"
#include "../epos/fnImpulse.h"
#include "../epos/ConvolutionFunction.h"
#include "Stopwatch.hpp"
#include "TestThreads.h"

using namespace std;
using namespace epos;
using namespace stopwatch;

void timeOfReading1000Segmenets()
{
	String path = "d:/proj/epos/";
	DataSetTxtRaw dataSet;
	dataSet.open(path + "POS18201279-Octacosane-C28-B/", 0, 0);
	Segment data;
	Stopwatch stopwatch;
	while (dataSet.read(data)) {

	}
	auto time = stopwatch.elapsed();
	std::cout << "time = " << time;
}

void runInpulseFunction(VectorFloat& time, fitfun_type fitfun)
{
	constexpr int na = 6;
	double a[na];
	a[1] = 1.0;
	a[2] = 1.0;
	a[3] = 0.0;
	a[4] = 1.0;
	a[5] = 0.0;
	double dyda[na];
	for (double t : time) {
		double y;
		fitfun(t, a, &y, dyda, na);
	}
}

void runInpulseFunctionInLoop(fitfun_type fitfun, int loopCount, const string& funName) {
	Stopwatch stopwatch;
	for (int i = 0; i < loopCount; ++i) {
		runInpulseFunction(Channel::time, fitfun);
	}
	auto time = stopwatch.elapsed();
	std::cout << funName << " time = " << time << endl;
}

typedef void AveVarFun(const VectorFloat&, Float&, Float&);

void runAveVarInLoop(AveVarFun avevar_fun, int loopCount, const string& funName) {
	VectorFloat data(epos::NbrSAMP, 0);
	Float ave, var;
	Stopwatch stopwatch;
	for (int i = 0; i < loopCount; ++i) {
		avevar_fun(data, ave, var);
	}
	auto time = stopwatch.elapsed();
	std::cout << funName << " time = " << time << endl;
}

template<class T>
void runInLoop(T fun, int loopCount, const string& funName) {
	VectorFloat data(epos::NbrSAMP, 0);
	Float ave, var;
	Stopwatch stopwatch;
	for (int i = 0; i < loopCount; ++i) {
		fun(data, ave, var);
	}
	auto time = stopwatch.elapsed();
	std::cout << funName << " time = " << time << endl;
}

template <class _Fn, class... _Args>
void runInLoop2(_Fn&& _Fx, _Args&&... _Ax)
{
	forward<_Fn>(_Fx)(forward<_Args>(_Ax)...);
}

void timeOfImpulseFunction()
{
	constexpr int loopCount(10000);
	runInpulseFunctionInLoop(ConvolutionFunction, loopCount, "ConvolutionFunction");
	runInpulseFunctionInLoop(fnImpulse, loopCount, "fnImpulse");
	runInpulseFunctionInLoop(fnImpulse, loopCount, "fnImpulse");
	runInpulseFunctionInLoop(ConvolutionFunction, loopCount, "ConvolutionFunction");
}

void f1(double val)
{

}

void f2(const vector<double>& v, double* res)
{
	*res = 0.0;
	unique_lock<mutex> lck{ TestThreads::m };
	for (auto& x : v)
		*res += x;
}

void timeOfAveVar()
{
	Float ave, var;
	//runInLoop2(epos::avevar, ave, var);

	//double res1;
	//vector<double> vec1{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	//runInLoop2(f1, res1);
	//runInLoop2(f2, vec1, &res1);

	constexpr int loopCount(10000);
	runInLoop<AveVarFun>(avevar, loopCount, "avevar");
	runInLoop<AveVarFun>(avevar, loopCount, "avevar");
}

int main()
{
	//timeOfImpulseFunction();
	//TestThreads::run();
	timeOfAveVar();
}