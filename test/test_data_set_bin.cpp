#include "catch.hpp"
#include "epos.hpp"

using namespace std;
using namespace epos;

TEST_CASE("DataSetBin open")
{
	bool res;
	String path = "d:/proj/epos/";
	DataSetBin dataSet;
	res = dataSet.open(path + "POS18201279-Octacosane-C28-B/", 0, 2);
	REQUIRE(res);
	if (res) {
		int count = 0;
		Segment data;
		while (dataSet.read(data)) {
			count++;
		}
		REQUIRE(count == 3);
	}
}