#include "catch.hpp"
#include "epos.hpp"

using namespace std;
using namespace epos;

TEST_CASE("DataSetTxt open")
{
	bool res;
	String path = "d:/proj/epos/";
	DataSetTxtRaw dataSet;
	res = dataSet.open(path + "POS18201279-Octacosane-C28-B/", 0, 0);
	REQUIRE(res);
	if (res) {
		Segment data;
		bool res = dataSet.read(data);
		REQUIRE(res == true);
		REQUIRE(data.ch1.data.size() == epos::NbrSAMP);
	}
}