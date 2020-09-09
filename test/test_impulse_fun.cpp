#include <catch.hpp>
#include <epos.hpp>
#include "impulse_fun.hpp"

using namespace std;
using namespace epos;

TEST_CASE("fnImpulse")
{
	constexpr int na = 6;
	double a[na];
	//a[1] - s
	//a[2] - tau
	//a[3] - dt - przesuniecie
	//a[4] - J - amplituda
	//a[5] - C - sta³a
	//dyda[1] - WartoscPochodnejPoS
	//dyda[2] - WartoscPochodnejPoTau
	//dyda[3] - WartoscPochodnejPoDT
	//dyda[4] - WartoscPochodnejPoJ
	//dyda[5] - WartoscPochodnejPoC
	a[1] = 1.0;
	a[2] = 1.0;
	a[3] = 0.0;
	a[4] = 1.0;
	a[5] = 0.0;
	double dyda1[na], dyda2[na];
	double arrx[3] = { -1.0, 0.0, 1.0 };
	for (double x : arrx) {
		double y1, y2;
		fnImpulse(x, a, &y1, dyda1, na);
		ConvolutionFunction(x, a, &y2, dyda2, na);
		REQUIRE(y1 == y2);
		constexpr double prec = 1e-6;
		REQUIRE(are_equal(dyda1[1], dyda2[1], prec));
		REQUIRE(are_equal(dyda1[2], dyda2[2], prec));
		REQUIRE(are_equal(dyda1[3], dyda2[3], prec));
		REQUIRE(are_equal(dyda1[4], dyda2[4], prec));
		REQUIRE(are_equal(dyda1[5], dyda2[5], prec));
	}
}