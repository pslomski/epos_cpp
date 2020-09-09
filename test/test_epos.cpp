#include "catch.hpp"
#include <epos.hpp>
#include "utils.hpp"

using namespace epos;
using namespace std;

TEST_CASE("sgn")
{
	// for T = int
	REQUIRE(sgn(0) == 0);

	REQUIRE(sgn(-1) == -1);
	REQUIRE(sgn(numeric_limits<int>::min()) == -1);

	REQUIRE(sgn(1) == 1);
	REQUIRE(sgn(numeric_limits<int>::max()) == 1);

	// for T = double
	REQUIRE(sgn(0.0) == 0);

	REQUIRE(sgn(-numeric_limits<double>::min()) == -1); // negative of minimum normalized positive value
	REQUIRE(sgn(-1.0) == -1);
	REQUIRE(sgn(-numeric_limits<double>::max()) == -1);

	REQUIRE(sgn(numeric_limits<double>::min()) == 1); // minimum normalized positive value
	REQUIRE(sgn(1.0) == 1);
	REQUIRE(sgn(numeric_limits<double>::max()) == 1);
}

TEST_CASE("inRange")
{
	REQUIRE(isInRange(0, 0, 0) == true);
	REQUIRE(isInRange(0, -1, 1) == true);
	REQUIRE(isInRange(0, 1, 2) == false);
	REQUIRE(isInRange(0, -2, -1) == false);
}

TEST_CASE("SeriesStatistics 1")
{
	VectorFloat v = { 1, 2, 3 };
	SeriesStatistics ss;
	ss.compute_from_range(v.begin(), v.end());
	REQUIRE(are_equal(ss.arithmetic_mean, 2.0));
	REQUIRE(are_equal(ss.variance, 1.0));
	REQUIRE(are_equal(ss.standard_deviation, 1.0));
}

TEST_CASE("SeriesStatistics 2")
{
	VectorFloat v = { 1.1, 2, 2.9 };
	SeriesStatistics ss;
	ss.compute_from_range(v.begin(), v.end());
	REQUIRE(are_equal(ss.arithmetic_mean, 2.0));
	REQUIRE(are_equal(ss.variance, 0.81));
	REQUIRE(are_equal(ss.standard_deviation, 0.9));
}

//auto it = std::next(v.begin(), index);
TEST_CASE("stddev")
{
	Float m, sdev;
	VectorFloat v;

	m = arithmetic_mean(v.begin(), v.end());
	REQUIRE(are_equal(m, 0.0));
	sdev = variance(v.begin(), v.end());
	REQUIRE(are_equal(sdev, 0.0));

	v = { 1, 2, 3 };
	m = arithmetic_mean(v.begin(), v.end());
	REQUIRE(are_equal(m, 2.0));
	sdev = variance(v.begin(), v.end());
	REQUIRE(are_equal(sdev, 1.0));

	v = { -1, -2, 1, 2, 3 };
	m = arithmetic_mean(std::next(v.begin(), 2), std::next(v.begin(), 2 + 3));
	REQUIRE(are_equal(m, 2.0));
	sdev = variance(std::next(v.begin(), 2), std::next(v.begin(), 2 + 3));
	REQUIRE(are_equal(sdev, 1.0));

	v = { -1, -2, 1, 2, 3 };
	m = arithmetic_mean(std::next(v.end(), -3), v.end());
	REQUIRE(are_equal(m, 2.0));
	sdev = variance(std::next(v.end(), -3), v.end());
	REQUIRE(are_equal(sdev, 1.0));

	v = { -1, -2, 1, 2, 3, -1, -2 };
	m = arithmetic_mean(std::next(v.begin(), 2), std::next(v.begin(), 2 + 3));
	REQUIRE(are_equal(m, 2.0));
	sdev = variance(std::next(v.begin(), 2), std::next(v.begin(), 2 + 3));
	REQUIRE(are_equal(sdev, 1.0));

	m = arithmetic_mean(v.begin(), v.end());
	REQUIRE(are_equal(m, 0.0));
	sdev = variance(v.begin(), v.end());
	REQUIRE(are_equal(sdev, 4.0));

	avevar(v.begin(), v.end(), m, sdev);
	REQUIRE(are_equal(m, 0.0));
	REQUIRE(are_equal(sdev, 4.0));
}

TEST_CASE("arithmetic mean")
{
	VectorFloat v = { 1, 2, 3 };
	Float mean = epos::arithmetic_mean(v.begin(), v.end());
	REQUIRE(are_equal(mean, 2.0));
	REQUIRE(v.end() - v.begin() == 3);
	v.clear();
	REQUIRE(v.end() - v.begin() == 0);
}

/*TEST_CASE("Lookup Table"){
	LookupTable ltexp(Channel::time.size() * 2, Channel::time.front(), Channel::time.back(), exp);
	for(Float t: Channel::time)
		REQUIRE(ltexp.flut(t) == exp(t));
}*/
