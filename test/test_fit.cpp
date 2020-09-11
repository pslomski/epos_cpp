#include <catch.hpp>
#include <nr.hpp>
#include <nrutil.hpp>
#include <epos.hpp>
#include <fit.hpp>

using namespace epos;

void linear_fun(double x, double a[], double* y, double dyda[], int na)
{
	*y = a[1] * x + a[2];
	dyda[1] = a[1];
	dyda[2] = 1;
}

void fill_matrix(double** alpha, int ma)
{
	for (int j = 1; j <= ma; j++) {
		for (int k = 1; k <= ma; k++)
			alpha[j][k] = j * k;
	}
}

TEST_CASE("fill_matrix")
{
	int ma = 2;
	Float** m = dmatrix(1, ma, 1, ma);
	fill_matrix(m, ma);
	REQUIRE(m[1][1] == 1.0);
	REQUIRE(m[1][2] == 2.0);
	REQUIRE(m[2][1] == 2.0);
	REQUIRE(m[2][2] == 4.0);
	free_dmatrix(m, 1, ma, 1, ma);
}

TEST_CASE("mrqmin")
{
	//return;
	Float x[3] = { 0.0, 0.0, 1.0 };
	Float y[3] = { 0.0, 0.0, 1.0 };
	Float sig[3] = { 1.0, 1.0, 1.0 };
	int ndata = 2;
	int ma = 2;
	Float a[3] = { 0.0, 1.1, 0.1 };
	int ia[3] = { 1, 1, 1 };
	double** covar = dmatrix(1, ma, 1, ma);
	double** alpha = dmatrix(1, ma, 1, ma);
	Float chisq = 0;
	Float lambda = -1.0;
	mrqmin(x, y, sig, ndata, a, ia, ma, covar, alpha, &chisq, linear_fun, &lambda);
	mrqmin(x, y, sig, ndata, a, ia, ma, covar, alpha, &chisq, linear_fun, &lambda);
	lambda = 0.0;
	mrqmin(x, y, sig, ndata, a, ia, ma, covar, alpha, &chisq, linear_fun, &lambda);
	free_dmatrix(alpha, 1, ma, 1, ma);
	free_dmatrix(covar, 1, ma, 1, ma);
}

TEST_CASE("linear_fun")
{
	Float x, y;
	VectorFloat a{ 0, 1.0, 0.0 };
	VectorFloat dyda{ 0, 0, 0 };

	x = 0.0;
	linear_fun(x, a.data(), &y, dyda.data(), 2);
	REQUIRE(y == Approx(0.0));
	REQUIRE(dyda[1] == Approx(a[1]));
	REQUIRE(dyda[2] == Approx(1.0));

	x = 1.0;
	linear_fun(x, a.data(), &y, dyda.data(), 2);
	REQUIRE(y == Approx(1.0));
	REQUIRE(dyda[1] == Approx(a[1]));
	REQUIRE(dyda[2] == Approx(1.0));

	a = { 0, -1.0, 0.0 };
	x = 0.0;
	linear_fun(x, a.data(), &y, dyda.data(), 2);
	REQUIRE(y == Approx(0.0));
	REQUIRE(dyda[1] == Approx(a[1]));
	REQUIRE(dyda[2] == Approx(1.0));

	x = 1.0;
	linear_fun(x, a.data(), &y, dyda.data(), 2);
	REQUIRE(y == Approx(-1.0));
	REQUIRE(dyda[1] == Approx(a[1]));
	REQUIRE(dyda[2] == Approx(1.0));
}

TEST_CASE("Fit linear_fun")
{
	return;
	VectorFloat x{ 0, 1, 2, 3 };
	VectorFloat y{ 0, 1, 2, 3 };
	VectorFloat a{ 0, 1.1, 0.2 };
	VectorFloat sig_a{ 0, 0, 0 };

	REQUIRE(epos::fit_fun(x, y, a, sig_a, linear_fun) == true);
	REQUIRE(a[0] == 1.0);
	REQUIRE(a[1] == 0.0);
}