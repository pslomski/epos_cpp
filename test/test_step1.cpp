#include "catch.hpp"
#include "epos.hpp"
#include "step1.hpp"

using namespace epos;
using namespace std;

class UnitTestStep1Processor: private Step1 {
public:
	bool test_checkBackground(epos::Channel& channel, epos::Float beginDCcoeff, epos::Float endDCcoeff) {
		return check_background(channel, beginDCcoeff, endDCcoeff);
	}
};

TEST_CASE("Step1Processor: all zero data")
{
	// Prapare artificial data of all zeros
	Channel ch;
	ch.init(epos::NbrSAMP);
	for (int i = 0; i < epos::NbrSAMP; ++i)
		ch.add(0.0);

	UnitTestStep1Processor step1Processor;
	REQUIRE(step1Processor.test_checkBackground(ch, 2.5, 2.5));
}

TEST_CASE("Step1Processor: begin of data lower than background")
{
	// Prapare artificial data
	Channel ch;
	ch.init(epos::NbrSAMP);
	for (int i = 0; i < epos::NbrSAMP; ++i) {
		if (i < 50)
			ch.add(-0.2);
		else
			ch.add(0.0);
	}

	UnitTestStep1Processor step1Processor;
	REQUIRE(step1Processor.test_checkBackground(ch, 2.5, 2.5) == false);
}

TEST_CASE("Step1Processor: end of data lower than background")
{
	// Prapare artificial data
	Channel ch;
	ch.init(epos::NbrSAMP);
	for (int i = 0; i < epos::NbrSAMP; ++i) {
		if (i > epos::NbrSAMP - 50)
			ch.add(-0.2);
		else
			ch.add(0.0);
	}

	UnitTestStep1Processor step1Processor;
	REQUIRE(step1Processor.test_checkBackground(ch, 2.5, 2.5) == false);
}

TEST_CASE("Step1Processor: real valid data")
{
	string rootDir = "D:/proj/epos/POS18201279-Octacosane-C28-B/";
	DataSetBin dataSet;
	if (dataSet.open(rootDir, 0, 0)) {
		Segment seg;
		dataSet.read(seg);
		UnitTestStep1Processor step1Processor;
		REQUIRE(step1Processor.test_checkBackground(seg.ch1, 2.5, 2.5));
	}
	else
		REQUIRE(false);
}