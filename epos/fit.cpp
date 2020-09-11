/* Driver for routine mrqmin */

#include <stdio.h>
#include <math.h>
#include <nr.hpp>
#include <nrutil.hpp>
#include <epos.hpp>
#include "fit.hpp"

#define NPT 100
#define MA 6
#define SPREAD 0.001

// using namespace std;
using namespace epos;

/* For a given value of x, and input parameters a[], this function
returns the resulting value of y and the vector dyda[] which
contains the derivatives of y with respect to all of the a[]
parameters */
void fgauss(Float x, Float a[], Float* y, Float dyda[], int na)
{
	int i;
	Float fac, ex, arg;
	*y = 0.0;
	for (i = 1; i <= na - 1; i += 3) {
		arg = (x - a[i + 1]) / a[i + 2];
		ex = exp(-arg * arg);
		fac = a[i] * ex * 2.0 * arg;
		*y += a[i] * ex;
		dyda[i] = ex;
		dyda[i + 1] = fac / a[i + 2];
		dyda[i + 2] = fac * arg / a[i + 2];
	}
	return;
}

int fit_()
{
	long idum = (-911);
	int i, * ia, iter, itst, j, k, mfit = MA;
	double alamda, chisq, ochisq, * x, * y, * sig, ** covar, ** alpha;
	static double a[MA + 1] = { 0.0,5.0,2.0,3.0,2.0,5.0,3.0 };
	static double gues[MA + 1] = { 0.0,4.5,2.2,2.8,2.5,4.9,2.8 };

	ia = ivector(1, MA);
	x = dvector(1, NPT);
	y = dvector(1, NPT);
	sig = dvector(1, NPT);
	covar = dmatrix(1, MA, 1, MA);
	alpha = dmatrix(1, MA, 1, MA);
	/* First try a sum of two Gaussians */
	for (i = 1; i <= NPT; i++) {
		x[i] = 0.1 * i;
		y[i] = 0.0;
		for (j = 1; j <= MA; j += 3) {
			y[i] += a[j] * exp(-DSQR((x[i] - a[j + 1]) / a[j + 2]));
		}
		y[i] *= (1.0 + SPREAD * 0 /*gasdev(&idum)*/);
		sig[i] = SPREAD * y[i];
	}
	for (i = 1; i <= mfit; i++) ia[i] = 1;
	for (i = 1; i <= MA; i++) a[i] = gues[i];
	for (iter = 1; iter <= 2; iter++) {
		alamda = -1;
		mrqmin(x, y, sig, NPT, a, ia, MA, covar, alpha, &chisq, fgauss, &alamda);
		k = 1;
		itst = 0;
		for (;;) {
			printf("\n%s %2d %17s %10.4f %10s %9.2e\n",
				"Iteration #", k, "chi-squared:", chisq, "alamda:", alamda);
			printf("%8s %8s %8s %8s %8s %8s\n",
				"a[1]", "a[2]", "a[3]", "a[4]", "a[5]", "a[6]");
			for (i = 1; i <= 6; i++) printf("%9.4f", a[i]);
			printf("\n");
			k++;
			ochisq = chisq;
			mrqmin(x, y, sig, NPT, a, ia, MA, covar, alpha, &chisq, fgauss, &alamda);
			if (chisq > ochisq)
				itst = 0;
			else if (fabs(ochisq - chisq) < 0.1)
				itst++;
			if (itst < 4) continue;
			alamda = 0.0;
			mrqmin(x, y, sig, NPT, a, ia, MA, covar, alpha, &chisq, fgauss, &alamda);
			printf("\nUncertainties:\n");
			for (i = 1; i <= 6; i++) printf("%9.4f", sqrt(covar[i][i]));
			printf("\n");
			printf("\nExpected results:\n");
			printf(" %7.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n",
				5.0, 2.0, 3.0, 2.0, 5.0, 3.0);
			break;
		}
		if (iter == 1) {
			printf("press return to continue with constraint\n");
			(void)getchar();
			printf("holding a[2] and a[5] constant\n");
			for (j = 1; j <= MA; j++) a[j] += 0.1;
			a[2] = 2.0;
			ia[2] = 0;
			a[5] = 5.0;
			ia[5] = 0;
		}
	}
	free_dmatrix(alpha, 1, MA, 1, MA);
	free_dmatrix(covar, 1, MA, 1, MA);
	free_dvector(sig, 1, NPT);
	free_dvector(y, 1, NPT);
	free_dvector(x, 1, NPT);
	free_ivector(ia, 1, MA);
	return 0;
}

double compute_sigma(double chisq, double covar, int ndata, int ma)
{
	return sqrt(chisq * covar / (double)(ndata - ma));
}

bool fit1(VectorFloat& x, VectorFloat& y, VectorFloat& a, VectorFloat& sig_a, fitfun_type function)
{
	bool res = true;
	const size_t MaxSteps = 100; // maksymalna liczba krokow fitowania
	VectorFloat sig(x.size(), 1.0); // Odchylenia standardowe danych y. Przyjmujemy = 1
	VectorInt ia(x.size(), 1);
	int ndata = static_cast<int>(x.size() - 1);
	int ma = static_cast<int>(a.size() - 1); // coefficient count
	double** covar = dmatrix(1, ma, 1, ma);
	double** alpha = dmatrix(1, ma, 1, ma);
	int step = 0;
	double lambda = -1.0; // set alamda<0 for initialization
	double chisq = 0.0;

	mrqmin(x.data(), y.data(), sig.data(), ndata, a.data(), ia.data(), ma, covar, alpha, &chisq, function, &lambda);

	std::cout << "#step " << step << ", chisq=" << chisq << ", lambda=" << lambda << std::endl;
	double ochisq = chisq + 2.0;

	/* Iterating with mrqmin until relative change in chisq is less than 0.001 */
	while ((fabs((ochisq - chisq) / ochisq) > 0.00001) || (lambda >= 0.1)) {
		step++;
		ochisq = chisq;
		mrqmin(x.data(), y.data(), sig.data(), ndata, a.data(), ia.data(), ma, covar, alpha, &chisq, function, &lambda);

		//std::cout << "#step " << step << ", chisq=" << chisq << ", lambda=" << lambda << std::endl;
		if (step > MaxSteps) {
			res = false; // brak zbieznosci w zadanej liczbie krokow
			break;
		}
		if (lambda > 2) {
			res = false; // rozbieznosc
			break;
		}
	}

	if (res) {
		/* Here we run mrqmin one more time with alambda = 0.0 to get covariant matrix, in case we want the errors */
		lambda = 0.0;//last call with 0
		mrqmin(x.data(), y.data(), sig.data(), ndata, a.data(), ia.data(), ma, covar, alpha, &chisq, function, &lambda);

		a[0] = chisq;//wyrzucane jest na zewnatrz Chi bedace surowa suma dla wszystkich punktow

		sig_a.resize(a.size());
		sig_a[0] = chisq;
		for (size_t i = 1; i < a.size(); i++) {
			//obliczanie bledu dopasowania parametrow a
			sig_a[i] = compute_sigma(chisq, covar[i][i], ndata, ma);
			//cout << "#a[" << i << "]=" << a[i] << ", sig=" << sig_a[i] << endl;
		}
	}
	free_dmatrix(covar, 1, ma, 1, ma);
	free_dmatrix(alpha, 1, ma, 1, ma);

	return res;
}

bool epos::fit_fun(VectorFloat& x, VectorFloat& y, VectorFloat& a, VectorFloat& sig_a, fitfun_type function)
{
	bool res = true;
	const size_t MaxSteps = 100; // maksymalna liczba krokow fitowania
	VectorFloat sig(x.size(), 1.0); // Odchylenia standardowe danych y. Przyjmujemy = 1
	VectorInt ia(x.size(), 1);
	int ndata = static_cast<int>(x.size() - 1);
	int ma = static_cast<int>(a.size() - 1); // coefficient count
	double** covar = dmatrix(1, ma, 1, ma);
	double** alpha = dmatrix(1, ma, 1, ma);
	int step = 0;
	double lambda = -1.0; // set alamda<0 for initialization
	double chisq = 0.0;

	mrqmin(x.data(), y.data(), sig.data(), ndata, a.data(), ia.data(), ma, covar, alpha, &chisq, function, &lambda);

	//std::cout << "#step " << step << ", chisq=" << chisq << ", lambda=" << lambda << std::endl;
	double ochisq = chisq + 2.0;

	/* Iterating with mrqmin until relative change in chisq is less than 0.001 */
	while ((fabs((ochisq - chisq) / ochisq) > 0.00001) || (lambda >= 0.1)) {
		step++;
		ochisq = chisq;
		mrqmin(x.data(), y.data(), sig.data(), ndata, a.data(), ia.data(), ma, covar, alpha, &chisq, function, &lambda);

		//std::cout << "#step " << step << ", chisq=" << chisq << ", lambda=" << lambda << std::endl;
		if (step > MaxSteps) {
			res = false; // brak zbieznosci w zadanej liczbie krokow
			break;
		}
		if (lambda > 2) {
			res = false; // rozbieznosc
			break;
		}
	}

	if (res) {
		/* Here we run mrqmin one more time with lambda = 0.0 to get covariant matrix, in case we want the errors */
		lambda = 0.0;//last call with 0
		mrqmin(x.data(), y.data(), sig.data(), ndata, a.data(), ia.data(), ma, covar, alpha, &chisq, function, &lambda);

		a[0] = chisq;//wyrzucane jest na zewnatrz Chi bedace surowa suma dla wszystkich punktow

		sig_a.resize(a.size());
		sig_a[0] = chisq;
		for (size_t i = 1; i < a.size(); i++) {
			//obliczanie bledu dopasowania parametrow a
			sig_a[i] = compute_sigma(chisq, covar[i][i], ndata, ma);
			//cout << "#a[" << i << "]=" << a[i] << ", sig=" << sig_a[i] << endl;
		}
	}
	free_dmatrix(covar, 1, ma, 1, ma);
	free_dmatrix(alpha, 1, ma, 1, ma);

	return res;
}