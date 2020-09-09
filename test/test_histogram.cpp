#include "catch.hpp"
#include "epos.hpp"
#include "histogram.hpp"

using namespace epos;
using namespace std;

TEST_CASE("Histogram")
{
	size_t count = 2;
	Histogram hist;
	hist.init(count, 0, 10);
	const VectorInt& bins = hist.getBins();
	REQUIRE(bins.size() == count);
	
	REQUIRE(hist.add(-0.00001) == false);
	REQUIRE(hist.add(-0.001) == false);
	REQUIRE(hist.add(-1) == false);
	REQUIRE(hist.add(10) == false);
	REQUIRE(hist.add(10.00001) == false);
	REQUIRE(hist.add(10.001) == false);
	REQUIRE(hist.add(11.0) == false);

	REQUIRE(hist.add(0.0) == true);
	REQUIRE(bins[0] == 1);
	REQUIRE(bins[1] == 0);

	REQUIRE(hist.add(5.0) == true);
	REQUIRE(bins[0] == 1);
	REQUIRE(bins[1] == 1);
	
	REQUIRE(hist.add(4.9999) == true);
	REQUIRE(bins[0] == 2);
	REQUIRE(bins[1] == 1);

	REQUIRE(hist.add(5.0001) == true);
	REQUIRE(bins[0] == 2);
	REQUIRE(bins[1] == 2);

	REQUIRE(hist.add(9.9999999) == true);
	REQUIRE(bins[0] == 2);
	REQUIRE(bins[1] == 3);

}

TEST_CASE("Histogram2")
{
	size_t count = 10;
	Histogram hist;
	hist.init(count, 0, 10);
	const VectorInt& bins = hist.getBins();
	REQUIRE(bins.size() == count);
	hist.add(0.1);
	//REQUIRE(1 == hist.getMaxHist());
	REQUIRE(0 == hist.getMaxHistIdx());
	hist.add(1.1);
	//REQUIRE(1 == hist.getMaxHist());
	REQUIRE(0 == hist.getMaxHistIdx());
	hist.add(1.1);
	//REQUIRE(2 == hist.getMaxHist());
	REQUIRE(1 == hist.getMaxHistIdx());
	hist.add(9.1);
	hist.add(9.1);
	hist.add(9.1);
	//REQUIRE(3 == hist.getMaxHist());
	REQUIRE(9 == hist.getMaxHistIdx());
}

TEST_CASE("Histogram3")
{
	Histogram hist;
	REQUIRE_THROWS(hist.init(0, 0.0, 1.0));
	REQUIRE_THROWS(hist.init(10, 0.0, 0.0));
	REQUIRE_THROWS(hist.init(10, 1.0, -1.0));
}