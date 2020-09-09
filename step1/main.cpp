#define _USE_MATH_DEFINES

#include <assert.h>
#include <cmath>
#include <string>
#include <ostream>
#include <fstream>
#include "epos.hpp"
#include "cxxopts.hpp"
#include "step1.hpp"

using namespace std;
using namespace epos;

void step1(DataSet& dataSet)
{
	Step1 step1Processor;
	Segment seg;
	while (dataSet.read(seg)) {
		cout << "fileNr" << seg.fileNr << " fileSeg" << seg.segNr << endl;
		step1Processor.process_segment(seg);
	}
}

int main(int argc, char* argv[])
{
	string rootDir = "D:/proj/epos/POS18201279-Octacosane-C28-B/";
	DataSetBin dataSet;
	if (dataSet.open(rootDir, 0, 0)) {
		step1(dataSet);
		return 0;
	}
	else
		return -1;
}