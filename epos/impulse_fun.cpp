#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>

using namespace std;

//wartosc funkcji i pochodnej w danym punkcie dla danych parametrow funkcji
struct DaneWyjscioweDoFitu {
	double WartoscFunkcji;
	double WartoscPochodnejPoTau;
	double WartoscPochodnejPoS;
	double WartoscPochodnejPoDT;
	double WartoscPochodnejPoJ;
	double WartoscPochodnejPoC;
	DaneWyjscioweDoFitu() { DaneWyjscioweDoFitu(0.0, 0.0, 0.0, 0.0, 0.0, 0.0); }
	DaneWyjscioweDoFitu(double wf, double wppt, double wpps, double wppdt, double wppj, double wppc) :
		WartoscFunkcji(wf), WartoscPochodnejPoTau(wppt), WartoscPochodnejPoS(wpps),
		WartoscPochodnejPoDT(wppdt), WartoscPochodnejPoJ(wppj), WartoscPochodnejPoC(wppc) {}
};

static double Y(double x, double s, double tau, double dt)
{
	return exp((s * s) / (4 * (tau * tau))) * erfc(s / (2 * tau) - (x - dt) / s) * exp(-(x - dt) / tau);
}

static double impulse(double x, double s, double tau, double dt, double J, double C)
{
	constexpr double deltax(0.5);
	return -0.5 * J * (
		Y(x, s, tau, dt) + erfc((x - dt) / s) -
		Y(x + deltax, s, tau, dt) - erfc((x + deltax - dt) / s))
		+ C; // uwzgledniony MINUS przed funkcja bo impulsy sa ujemne
}

static DaneWyjscioweDoFitu compute_all(double x, double s, double tau, double dt, double J, double C)
{
	constexpr double delta = 1e-6;
	DaneWyjscioweDoFitu data;
	data.WartoscFunkcji = impulse(x, s, tau, dt, J, C);
	data.WartoscPochodnejPoC = (impulse(x, s, tau, dt, J, C + delta) - impulse(x, s, tau, dt, J, C)) / delta;
	data.WartoscPochodnejPoDT = (impulse(x, s, tau, dt + delta, J, C) - impulse(x, s, tau, dt, J, C)) / delta;
	data.WartoscPochodnejPoJ = (impulse(x, s, tau, dt, J + delta, C) - impulse(x, s, tau, dt, J, C)) / delta;
	data.WartoscPochodnejPoS = (impulse(x, s + delta, tau, dt, J, C) - impulse(x, s, tau, dt, J, C)) / delta;
	data.WartoscPochodnejPoTau = (impulse(x, s, tau + delta, dt, J, C) - impulse(x, s, tau, dt, J, C)) / delta;
	return data;
}

static DaneWyjscioweDoFitu WszystkoDoFitu(double x, double s, double tau, double dt, double J, double C)
{
	constexpr double sqrtPi(1.7724538509055160272981674833411); // sqrt(M_PI);
	constexpr double deltax(0.5);
	double Funkcja, PochodnaPoTau, PochodnaPoS, PochodnaPoDT, PochodnaPoJ, PochodnaPoC; // wartosc funkcji i pochodnych w danym punkcie
	double tausq = tau * tau;
	double ssq = s * s;

	Funkcja = -0.5 * J * (
		Y(x, s, tau, dt) + erfc((x - dt) / s) -
		Y(x + deltax, s, tau, dt) - erfc((x + deltax - dt) / s))
		+ C; // uwzgledniony MINUS przed funkcja bo impulsy sa ujemne

	PochodnaPoTau = -0.5 * J * (
		s / (sqrt(M_PI) * pow(tau, 2)) * exp(pow(s, 2) / (4 * pow(tau, 2))) * exp(-pow(s / (2 * tau) - (x - dt) / s, 2)) * exp(-(x - dt) / tau) + (x - dt) / (pow(tau, 2)) * Y(x, s, tau, dt) - pow(s, 2) / (2 * pow(tau, 3)) * Y(x, s, tau, dt) -
		s / (sqrt(M_PI) * pow(tau, 2)) * exp(pow(s, 2) / (4 * pow(tau, 2))) * exp(-pow(s / (2 * tau) - (x + deltax - dt) / s, 2)) * exp(-(x + deltax - dt) / tau) - (x + deltax - dt) / (pow(tau, 2)) * Y(x + deltax, s, tau, dt) + pow(s, 2) / (2 * pow(tau, 3)) * Y(x + deltax, s, tau, dt)
		);

	PochodnaPoS = -0.5 * J * (
		s / (2 * pow(tau, 2)) * Y(x, s, tau, dt) - 2 / sqrt(M_PI) * ((x - dt) / pow(s, 2) + 1 / (2 * tau)) * exp(-pow(s / (2 * tau) - (x - dt) / s, 2)) * exp(pow(s, 2) / (4 * pow(tau, 2))) * exp(-(x - dt) / tau) -
		s / (2 * pow(tau, 2)) * Y(x + deltax, s, tau, dt) + 2 / sqrt(M_PI) * ((x + deltax - dt) / pow(s, 2) + 1 / (2 * tau)) * exp(-pow(s / (2 * tau) - (x + deltax - dt) / s, 2)) * exp(pow(s, 2) / (4 * pow(tau, 2))) * exp(-(x + deltax - dt) / tau) +
		2 * (x - dt) / (sqrt(M_PI) * pow(s, 2)) * exp(-pow(x - dt, 2) / pow(s, 2)) -
		2 * (x + deltax - dt) / (sqrt(M_PI) * pow(s, 2)) * exp(-pow(x + deltax - dt, 2) / pow(s, 2))
		);

	PochodnaPoDT = -0.5 * J * (
		1 / tau * Y(x, s, tau, dt) - 2 / (s * sqrt(M_PI)) * exp(pow(s, 2) / (4 * pow(tau, 2))) * exp(-(x - dt) / tau) * exp(-pow(s / (2 * tau) - (x - dt) / s, 2))
		- 1 / tau * Y(x + deltax, s, tau, dt) + 2 / (s * sqrt(M_PI)) * exp(pow(s, 2) / (4 * pow(tau, 2))) * exp(-(x + deltax - dt) / tau) * exp(-pow(s / (2 * tau) - (x + deltax - dt) / s, 2)) +
		2 / (s * sqrt(M_PI)) * exp(-pow(x - dt, 2) / pow(s, 2)) -
		2 / (s * sqrt(M_PI)) * exp(-pow(x + deltax - dt, 2) / pow(s, 2))
		);

	PochodnaPoJ = -0.5 * (Y(x, s, tau, dt) - Y(x + deltax, s, tau, dt) + erfc((x - dt) / s) - erfc((x + deltax - dt) / s));

	PochodnaPoC = 1;

	return DaneWyjscioweDoFitu(Funkcja, PochodnaPoTau, PochodnaPoS, PochodnaPoDT, PochodnaPoJ, PochodnaPoC);
}

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
void fnImpulse(double x, double a[], double* y, double dyda[], int na)
{
	DaneWyjscioweDoFitu Oblicz = compute_all(x, a[1], a[2], a[3], a[4], a[5]);

	*y = Oblicz.WartoscFunkcji;
	dyda[1] = Oblicz.WartoscPochodnejPoS;
	dyda[2] = Oblicz.WartoscPochodnejPoTau;
	dyda[3] = Oblicz.WartoscPochodnejPoDT;
	dyda[4] = Oblicz.WartoscPochodnejPoJ;
	dyda[5] = Oblicz.WartoscPochodnejPoC;
}

void ConvolutionFunction(double x, double a[], double* y, double dyda[], int na)
{
	DaneWyjscioweDoFitu Oblicz = WszystkoDoFitu(x, a[1], a[2], a[3], a[4], a[5]);

	*y = Oblicz.WartoscFunkcji;
	dyda[1] = Oblicz.WartoscPochodnejPoS;
	dyda[2] = Oblicz.WartoscPochodnejPoTau;
	dyda[3] = Oblicz.WartoscPochodnejPoDT;
	dyda[4] = Oblicz.WartoscPochodnejPoJ;
	dyda[5] = Oblicz.WartoscPochodnejPoC;
}
