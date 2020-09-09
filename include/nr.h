#ifndef _NR_H_
#define _NR_H_

using Float = double;

#ifndef _FCOMPLEX_DECLARE_T_
typedef struct FCOMPLEX { Float r, i; } fcomplex;
#define _FCOMPLEX_DECLARE_T_
#endif /* _FCOMPLEX_DECLARE_T_ */

#ifndef _ARITHCODE_DECLARE_T_
typedef struct {
	unsigned long* ilob, * iupb, * ncumfq, jdif, nc, minint, nch, ncum, nrad;
} arithcode;
#define _ARITHCODE_DECLARE_T_
#endif /* _ARITHCODE_DECLARE_T_ */

#ifndef _HUFFCODE_DECLARE_T_
typedef struct {
	unsigned long* icod, * ncod, * left, * right, nch, nodemax;
} huffcode;
#define _HUFFCODE_DECLARE_T_
#endif /* _HUFFCODE_DECLARE_T_ */

#include <stdio.h>

void addint(double** uf, double** uc, double** res, int nf);
void airy(Float x, Float* ai, Float* bi, Float* aip, Float* bip);
void amebsa(Float** p, Float y[], int ndim, Float pb[], Float* yb, Float ftol, Float(*funk)(Float[]), int* iter, Float temptr);
void amoeba(Float** p, Float y[], int ndim, Float ftol, Float(*funk)(Float[]), int* iter);
Float amotry(Float** p, Float y[], Float psum[], int ndim, Float(*funk)(Float[]), int ihi, Float fac);
Float amotsa(Float** p, Float y[], Float psum[], int ndim, Float pb[], Float* yb, Float(*funk)(Float[]), int ihi, Float* yhi, Float fac);
void anneal(Float x[], Float y[], int iorder[], int ncity);
double anorm2(double** a, int n);
void arcmak(unsigned long nfreq[], unsigned long nchh, unsigned long nradd, arithcode* acode);
void arcode(unsigned long* ich, unsigned char** codep, unsigned long* lcode, unsigned long* lcd, int isign, arithcode* acode);
void arcsum(unsigned long iin[], unsigned long iout[], unsigned long ja, int nwk, unsigned long nrad, unsigned long nc);
void asolve(unsigned long n, double b[], double x[], int itrnsp);
void atimes(unsigned long n, double x[], double r[], int itrnsp);
void avevar(Float data[], unsigned long n, Float* ave, Float* var);
void balanc(Float** a, int n);
void banbks(Float** a, unsigned long n, int m1, int m2, Float** al, unsigned long indx[], Float b[]);
void bandec(Float** a, unsigned long n, int m1, int m2, Float** al, unsigned long indx[], Float* d);
void banmul(Float** a, unsigned long n, int m1, int m2, Float x[], Float b[]);
void bcucof(Float y[], Float y1[], Float y2[], Float y12[], Float d1, Float d2, Float** c);
void bcuint(Float y[], Float y1[], Float y2[], Float y12[], Float x1l, Float x1u, Float x2l, Float x2u, Float x1,
	Float x2, Float* ansy, Float* ansy1, Float* ansy2);
void beschb(double x, double* gam1, double* gam2, double* gampl, double* gammi);
Float bessi(int n, Float x);
Float bessi0(Float x);
Float bessi1(Float x);
void bessik(Float x, Float xnu, Float* ri, Float* rk, Float* rip, Float* rkp);
Float bessj(int n, Float x);
Float bessj0(Float x);
Float bessj1(Float x);
void bessjy(Float x, Float xnu, Float* rj, Float* ry, Float* rjp, Float* ryp);
Float bessk(int n, Float x);
Float bessk0(Float x);
Float bessk1(Float x);
Float bessy(int n, Float x);
Float bessy0(Float x);
Float bessy1(Float x);
Float beta(Float z, Float w);
Float betacf(Float a, Float b, Float x);
Float betai(Float a, Float b, Float x);
Float bico(int n, int k);
void bksub(int ne, int nb, int jf, int k1, int k2, Float*** c);
Float bnldev(Float pp, int n, long* idum);
Float brent(Float ax, Float bx, Float cx, Float(*f)(Float), Float tol, Float* xmin);
void broydn(Float x[], int n, int* check, void (*vecfunc)(int, Float[], Float[]));
void bsstep(Float y[], Float dydx[], int nv, Float* xx, Float htry,
	Float eps, Float yscal[], Float* hdid, Float* hnext,
	void (*derivs)(Float, Float[], Float[]));
void caldat(long julian, int* mm, int* id, int* iyyy);
void chder(Float a, Float b, Float c[], Float cder[], int n);
Float chebev(Float a, Float b, Float c[], int m, Float x);
void chebft(Float a, Float b, Float c[], int n, Float(*func)(Float));
void chebpc(Float c[], Float d[], int n);
void chint(Float a, Float b, Float c[], Float cint[], int n);
Float chixy(Float bang);
void choldc(Float** a, int n, Float p[]);
void cholsl(Float** a, int n, Float p[], Float b[], Float x[]);
void chsone(Float bins[], Float ebins[], int nbins, int knstrn,
	Float* df, Float* chsq, Float* prob);
void chstwo(Float bins1[], Float bins2[], int nbins, int knstrn,
	Float* df, Float* chsq, Float* prob);
void cisi(Float x, Float* ci, Float* si);
void cntab1(int** nn, int ni, int nj, Float* chisq,
	Float* df, Float* prob, Float* cramrv, Float* ccc);
void cntab2(int** nn, int ni, int nj, Float* h, Float* hx, Float* hy,
	Float* hygx, Float* hxgy, Float* uygx, Float* uxgy, Float* uxy);
void convlv(Float data[], unsigned long n, Float respns[], unsigned long m, int isign, Float ans[]);
void copy(double** aout, double** ain, int n);
void correl(Float data1[], Float data2[], unsigned long n, Float ans[]);
void cosft(Float y[], int n, int isign);
void cosft1(Float y[], int n);
void cosft2(Float y[], int n, int isign);
void covsrt(Float** covar, int ma, int ia[], int mfit);
void crank(unsigned long n, Float w[], Float* s);
void cyclic(Float a[], Float b[], Float c[], Float alpha, Float beta, Float r[], Float x[], unsigned long n);
void daub4(Float a[], unsigned long n, int isign);
Float dawson(Float x);
Float dbrent(Float ax, Float bx, Float cx, Float(*f)(Float), Float(*df)(Float), Float tol, Float* xmin);
void ddpoly(Float c[], int nc, Float x, Float pd[], int nd);
int decchk(char string[], int n, char* ch);
void derivs(Float x, Float y[], Float dydx[]);
Float df1dim(Float x);
void dfour1(double data[], unsigned long nn, int isign);
void dfpmin(Float p[], int n, Float gtol, int* iter, Float* fret, Float(*func)(Float[]), void (*dfunc)(Float[], Float[]));
Float dfridr(Float(*func)(Float), Float x, Float h, Float* err);
void dftcor(Float w, Float delta, Float a, Float b, Float endpts[],
	Float* corre, Float* corim, Float* corfac);
void dftint(Float(*func)(Float), Float a, Float b, Float w, Float* cosint, Float* sinint);
void difeq(int k, int k1, int k2, int jsf, int is1, int isf, int indexv[], int ne, Float** s, Float** y);
void dlinmin(Float p[], Float xi[], int n, Float* fret, Float(*func)(Float[]), void (*dfunc)(Float[], Float[]));
double dpythag(double a, double b);
void drealft(double data[], unsigned long n, int isign);
void dsprsax(double sa[], unsigned long ija[], double x[], double b[], unsigned long n);
void dsprstx(double sa[], unsigned long ija[], double x[], double b[], unsigned long n);
void dsvbksb(double** u, double w[], double** v, int m, int n, double b[], double x[]);
void dsvdcmp(double** a, int m, int n, double w[], double** v);
void eclass(int nf[], int n, int lista[], int listb[], int m);
void eclazz(int nf[], int n, int (*equiv)(int, int));
Float ei(Float x);
void eigsrt(Float d[], Float** v, int n);
Float elle(Float phi, Float ak);
Float ellf(Float phi, Float ak);
Float ellpi(Float phi, Float en, Float ak);
void elmhes(Float** a, int n);
Float erfcc(Float x);
Float erff(Float x);
Float erffc(Float x);
void eulsum(Float* sum, Float term, int jterm, Float wksp[]);
Float evlmem(Float fdt, Float d[], int m, Float xms);
Float expdev(long* idum);
Float expint(int n, Float x);
Float f1(Float x);
Float f1dim(Float x);
Float f2(Float y);
Float f3(Float z);
Float factln(int n);
Float factrl(int n);
void fasper(Float x[], Float y[], unsigned long n, Float ofac, Float hifac,
	Float wk1[], Float wk2[], unsigned long nwk, unsigned long* nout,
	unsigned long* jmax, Float* prob);
void fdjac(int n, Float x[], Float fvec[], Float** df,
	void (*vecfunc)(int, Float[], Float[]));
void fgauss(Float x, Float a[], Float* y, Float dyda[], int na);
void fill0(double** u, int n);
void fit(Float x[], Float y[], int ndata, Float sig[], int mwt,
	Float* a, Float* b, Float* siga, Float* sigb, Float* chi2, Float* q);
void fitexy(Float x[], Float y[], int ndat, Float sigx[], Float sigy[],
	Float* a, Float* b, Float* siga, Float* sigb, Float* chi2, Float* q);
void fixrts(Float d[], int m);
void fleg(Float x, Float pl[], int nl);
void flmoon(int n, int nph, long* jd, Float* frac);
Float fmin(Float x[]);
void four1(Float data[], unsigned long nn, int isign);
void fourew(FILE* file[5], int* na, int* nb, int* nc, int* nd);
void fourfs(FILE* file[5], unsigned long nn[], int ndim, int isign);
void fourn(Float data[], unsigned long nn[], int ndim, int isign);
void fpoly(Float x, Float p[], int np);
void fred2(int n, Float a, Float b, Float t[], Float f[], Float w[],
	Float(*g)(Float), Float(*ak)(Float, Float));
Float fredin(Float x, int n, Float a, Float b, Float t[], Float f[], Float w[],
	Float(*g)(Float), Float(*ak)(Float, Float));
void frenel(Float x, Float* s, Float* c);
void frprmn(Float p[], int n, Float ftol, int* iter, Float* fret,
	Float(*func)(Float[]), void (*dfunc)(Float[], Float[]));
void ftest(Float data1[], unsigned long n1, Float data2[], unsigned long n2, Float* f, Float* prob);
Float gamdev(int ia, long* idum);
Float gammln(Float xx);
Float gammp(Float a, Float x);
Float gammq(Float a, Float x);
Float gasdev(long* idum);
void gaucof(int n, Float a[], Float b[], Float amu0, Float x[], Float w[]);
void gauher(Float x[], Float w[], int n);
void gaujac(Float x[], Float w[], int n, Float alf, Float bet);
void gaulag(Float x[], Float w[], int n, Float alf);
void gauleg(Float x1, Float x2, Float x[], Float w[], int n);
void gaussj(Float** a, int n, Float** b, int m);
void gcf(Float* gammcf, Float a, Float x, Float* gln);
Float golden(Float ax, Float bx, Float cx, Float(*f)(Float), Float tol,
	Float* xmin);
void gser(Float* gamser, Float a, Float x, Float* gln);
void hpsel(unsigned long m, unsigned long n, Float arr[], Float heap[]);
void hpsort(unsigned long n, Float ra[]);
void hqr(Float** a, int n, Float wr[], Float wi[]);
void hufapp(unsigned long index[], unsigned long nprob[], unsigned long n, unsigned long i);
void hufdec(unsigned long* ich, unsigned char* code, unsigned long lcode, unsigned long* nb, huffcode* hcode);
void hufenc(unsigned long ich, unsigned char** codep, unsigned long* lcode, unsigned long* nb, huffcode* hcode);
void hufmak(unsigned long nfreq[], unsigned long nchin, unsigned long* ilong, unsigned long* nlong, huffcode* hcode);
void hunt(Float xx[], unsigned long n, Float x, unsigned long* jlo);
void hypdrv(Float s, Float yy[], Float dyyds[]);
fcomplex hypgeo(fcomplex a, fcomplex b, fcomplex c, fcomplex z);
void hypser(fcomplex a, fcomplex b, fcomplex c, fcomplex z, fcomplex* series, fcomplex* deriv);
unsigned short icrc(unsigned short crc, unsigned char* bufptr, unsigned long len, short jinit, int jrev);
unsigned short icrc1(unsigned short crc, unsigned char onech);
unsigned long igray(unsigned long n, int is);
void iindexx(unsigned long n, long arr[], unsigned long indx[]);
void indexx(unsigned long n, Float arr[], unsigned long indx[]);
void interp(double** uf, double** uc, int nf);
int irbit1(unsigned long* iseed);
int irbit2(unsigned long* iseed);
void jacobi(Float** a, int n, Float d[], Float** v, int* nrot);
void jacobn(Float x, Float y[], Float dfdx[], Float** dfdy, int n);
long julday(int mm, int id, int iyyy);
void kendl1(Float data1[], Float data2[], unsigned long n, Float* tau, Float* z, Float* prob);
void kendl2(Float** tab, int i, int j, Float* tau, Float* z, Float* prob);
void kermom(double w[], double y, int m);
void ks2d1s(Float x1[], Float y1[], unsigned long n1, void (*quadvl)(Float, Float, Float*, Float*, Float*, Float*), Float* d1, Float* prob);
void ks2d2s(Float x1[], Float y1[], unsigned long n1, Float x2[], Float y2[], unsigned long n2, Float* d, Float* prob);
void ksone(Float data[], unsigned long n, Float(*func)(Float), Float* d, Float* prob);
void kstwo(Float data1[], unsigned long n1, Float data2[], unsigned long n2, Float* d, Float* prob);
void laguer(fcomplex a[], int m, fcomplex* x, int* its);
void lfit(Float x[], Float y[], Float sig[], int ndat, Float a[], int ia[], int ma, Float** covar, Float* chisq, void (*funcs)(Float, Float[], int));
void linbcg(unsigned long n, double b[], double x[], int itol, double tol, int itmax, int* iter, double* err);
void linmin(Float p[], Float xi[], int n, Float* fret, Float(*func)(Float[]));
void lnsrch(int n, Float xold[], Float fold, Float g[], Float p[], Float x[], Float* f, Float stpmax, int* check, Float(*func)(Float[]));
void load(Float x1, Float v[], Float y[]);
void load1(Float x1, Float v1[], Float y[]);
void load2(Float x2, Float v2[], Float y[]);
void locate(Float xx[], unsigned long n, Float x, unsigned long* j);
void lop(double** out, double** u, int n);
void lubksb(Float** a, int n, int* indx, Float b[]);
void ludcmp(Float** a, int n, int* indx, Float* d);
void machar(int* ibeta, int* it, int* irnd, int* ngrd, int* machep, int* negep, int* iexp, int* minexp, int* maxexp, Float* eps, Float* epsneg, Float* xmin, Float* xmax);
void matadd(double** a, double** b, double** c, int n);
void matsub(double** a, double** b, double** c, int n);
void medfit(Float x[], Float y[], int ndata, Float* a, Float* b, Float* abdev);
void memcof(Float data[], int n, int m, Float* xms, Float d[]);
int metrop(Float de, Float t);
void mgfas(double** u, int n, int maxcyc);
void mglin(double** u, int n, int ncycle);
Float midexp(Float(*funk)(Float), Float aa, Float bb, int n);
Float midinf(Float(*funk)(Float), Float aa, Float bb, int n);
Float midpnt(Float(*func)(Float), Float a, Float b, int n);
Float midsql(Float(*funk)(Float), Float aa, Float bb, int n);
Float midsqu(Float(*funk)(Float), Float aa, Float bb, int n);
void miser(Float(*func)(Float[]), Float regn[], int ndim, unsigned long npts, Float dith, Float* ave, Float* var);
void mmid(Float y[], Float dydx[], int nvar, Float xs, Float htot, int nstep, Float yout[], void (*derivs)(Float, Float[], Float[]));
void mnbrak(Float* ax, Float* bx, Float* cx, Float* fa, Float* fb, Float* fc, Float(*func)(Float));
void mnewt(int ntrial, Float x[], int n, Float tolx, Float tolf);
void moment(Float data[], int n, Float* ave, Float* adev, Float* sdev, Float* var, Float* skew, Float* curt);
void mp2dfr(unsigned char a[], unsigned char s[], int n, int* m);
void mpadd(unsigned char w[], unsigned char u[], unsigned char v[], int n);
void mpdiv(unsigned char q[], unsigned char r[], unsigned char u[],
	unsigned char v[], int n, int m);
void mpinv(unsigned char u[], unsigned char v[], int n, int m);
void mplsh(unsigned char u[], int n);
void mpmov(unsigned char u[], unsigned char v[], int n);
void mpmul(unsigned char w[], unsigned char u[], unsigned char v[], int n, int m);
void mpneg(unsigned char u[], int n);
void mppi(int n);
void mprove(Float** a, Float** alud, int n, int indx[], Float b[], Float x[]);
void mpsad(unsigned char w[], unsigned char u[], int n, int iv);
void mpsdv(unsigned char w[], unsigned char u[], int n, int iv, int* ir);
void mpsmu(unsigned char w[], unsigned char u[], int n, int iv);
void mpsqrt(unsigned char w[], unsigned char u[], unsigned char v[], int n, int m);
void mpsub(int* is, unsigned char w[], unsigned char u[], unsigned char v[], int n);

typedef void (*fitfun_type)(double, double[], double*, double[], int);

void mrqcof(Float x[], Float y[], Float sig[], int ndata, Float a[],
	int ia[], int ma, Float** alpha, Float beta[], Float* chisq,
	fitfun_type funcs
	//void (*funcs)(Float, Float[], Float*, Float[], int)
);

void mrqmin(Float x[], Float y[], Float sig[], int ndata, Float a[],
	int ia[], int ma, Float** covar, Float** alpha, Float* chisq,
	fitfun_type funcs,
	//void (*funcs)(Float, Float [], Float *, Float [], int),
	Float* alamda);
void newt(Float x[], int n, int* check,
	void (*vecfunc)(int, Float[], Float[]));
void odeint(Float ystart[], int nvar, Float x1, Float x2,
	Float eps, Float h1, Float hmin, int* nok, int* nbad,
	void (*derivs)(Float, Float[], Float[]),
	void (*rkqs)(Float[], Float[], int, Float*, Float, Float,
		Float[], Float*, Float*, void (*)(Float, Float[], Float[])));
void orthog(int n, Float anu[], Float alpha[], Float beta[], Float a[],
	Float b[]);
void pade(double cof[], int n, Float* resid);
void pccheb(Float d[], Float c[], int n);
void pcshft(Float a, Float b, Float d[], int n);
void pearsn(Float x[], Float y[], unsigned long n, Float* r, Float* prob,
	Float* z);
void period(Float x[], Float y[], int n, Float ofac, Float hifac,
	Float px[], Float py[], int np, int* nout, int* jmax, Float* prob);
void piksr2(int n, Float arr[], Float brr[]);
void piksrt(int n, Float arr[]);
void pinvs(int ie1, int ie2, int je1, int jsf, int jc1, int k,
	Float*** c, Float** s);
Float plgndr(int l, int m, Float x);
Float poidev(Float xm, long* idum);
void polcoe(Float x[], Float y[], int n, Float cof[]);
void polcof(Float xa[], Float ya[], int n, Float cof[]);
void poldiv(Float u[], int n, Float v[], int nv, Float q[], Float r[]);
void polin2(Float x1a[], Float x2a[], Float** ya, int m, int n,
	Float x1, Float x2, Float* y, Float* dy);
void polint(Float xa[], Float ya[], int n, Float x, Float* y, Float* dy);
void powell(Float p[], Float** xi, int n, Float ftol, int* iter, Float* fret,
	Float(*func)(Float[]));
void predic(Float data[], int ndata, Float d[], int m, Float future[], int nfut);
Float probks(Float alam);
void psdes(unsigned long* lword, unsigned long* irword);
void pwt(Float a[], unsigned long n, int isign);
void pwtset(int n);
Float pythag(Float a, Float b);
void pzextr(int iest, Float xest, Float yest[], Float yz[], Float dy[],
	int nv);
Float qgaus(Float(*func)(Float), Float a, Float b);
void qrdcmp(Float** a, int n, Float* c, Float* d, int* sing);
Float qromb(Float(*func)(Float), Float a, Float b);
Float qromo(Float(*func)(Float), Float a, Float b,
	Float(*choose)(Float(*)(Float), Float, Float, int));
void qroot(Float p[], int n, Float* b, Float* c, Float eps);
void qrsolv(Float** a, int n, Float c[], Float d[], Float b[]);
void qrupdt(Float** r, Float** qt, int n, Float u[], Float v[]);
Float qsimp(Float(*func)(Float), Float a, Float b);
Float qtrap(Float(*func)(Float), Float a, Float b);
Float quad3d(Float(*func)(Float, Float, Float), Float x1, Float x2);
void quadct(Float x, Float y, Float xx[], Float yy[], unsigned long nn,
	Float* fa, Float* fb, Float* fc, Float* fd);
void quadmx(Float** a, int n);
void quadvl(Float x, Float y, Float* fa, Float* fb, Float* fc, Float* fd);
Float ran0(long* idum);
Float ran1(long* idum);
Float ran2(long* idum);
Float ran3(long* idum);
Float ran4(long* idum);
void rank(unsigned long n, unsigned long indx[], unsigned long irank[]);
void ranpt(Float pt[], Float regn[], int n);
void ratint(Float xa[], Float ya[], int n, Float x, Float* y, Float* dy);
void ratlsq(double (*fn)(double), double a, double b, int mm, int kk,
	double cof[], double* dev);
double ratval(double x, double cof[], int mm, int kk);
Float rc(Float x, Float y);
Float rd(Float x, Float y, Float z);
void realft(Float data[], unsigned long n, int isign);
void rebin(Float rc, int nd, Float r[], Float xin[], Float xi[]);
void red(int iz1, int iz2, int jz1, int jz2, int jm1, int jm2, int jmf,
	int ic1, int jc1, int jcf, int kc, Float*** c, Float** s);
void relax(double** u, double** rhs, int n);
void relax2(double** u, double** rhs, int n);
void resid(double** res, double** u, double** rhs, int n);
Float revcst(Float x[], Float y[], int iorder[], int ncity, int n[]);
void reverse(int iorder[], int ncity, int n[]);
Float rf(Float x, Float y, Float z);
Float rj(Float x, Float y, Float z, Float p);
void rk4(Float y[], Float dydx[], int n, Float x, Float h, Float yout[],
	void (*derivs)(Float, Float[], Float[]));
void rkck(Float y[], Float dydx[], int n, Float x, Float h,
	Float yout[], Float yerr[], void (*derivs)(Float, Float[], Float[]));
void rkdumb(Float vstart[], int nvar, Float x1, Float x2, int nstep,
	void (*derivs)(Float, Float[], Float[]));
void rkqs(Float y[], Float dydx[], int n, Float* x,
	Float htry, Float eps, Float yscal[], Float* hdid, Float* hnext,
	void (*derivs)(Float, Float[], Float[]));
void rlft3(Float*** data, Float** speq, unsigned long nn1,
	unsigned long nn2, unsigned long nn3, int isign);
Float rofunc(Float b);
void rotate(Float** r, Float** qt, int n, int i, Float a, Float b);
void rsolv(Float** a, int n, Float d[], Float b[]);
void rstrct(double** uc, double** uf, int nc);
Float rtbis(Float(*func)(Float), Float x1, Float x2, Float xacc);
Float rtflsp(Float(*func)(Float), Float x1, Float x2, Float xacc);
Float rtnewt(void (*funcd)(Float, Float*, Float*), Float x1, Float x2,
	Float xacc);
Float rtsafe(void (*funcd)(Float, Float*, Float*), Float x1, Float x2,
	Float xacc);
Float rtsec(Float(*func)(Float), Float x1, Float x2, Float xacc);
void rzextr(int iest, Float xest, Float yest[], Float yz[], Float dy[], int nv);
void savgol(Float c[], int np, int nl, int nr, int ld, int m);
void score(Float xf, Float y[], Float f[]);
void scrsho(Float(*fx)(Float));
/* Float select(unsigned long k, unsigned long n, Float arr[]); */
Float selip(unsigned long k, unsigned long n, Float arr[]);
void shell(unsigned long n, Float a[]);
void shoot(int n, Float v[], Float f[]);
void shootf(int n, Float v[], Float f[]);
void simp1(Float** a, int mm, int ll[], int nll, int iabf, int* kp,
	Float* bmax);
void simp2(Float** a, int n, int l2[], int nl2, int* ip, int kp, Float* q1);
void simp3(Float** a, int i1, int k1, int ip, int kp);
void simplx(Float** a, int m, int n, int m1, int m2, int m3, int* icase,
	int izrov[], int iposv[]);
void simpr(Float y[], Float dydx[], Float dfdx[], Float** dfdy,
	int n, Float xs, Float htot, int nstep, Float yout[],
	void (*derivs)(Float, Float[], Float[]));
void sinft(Float y[], int n);
void slvsm2(double** u, double** rhs);
void slvsml(double** u, double** rhs);
void sncndn(Float uu, Float emmc, Float* sn, Float* cn, Float* dn);
double snrm(unsigned long n, double sx[], int itol);
void sobseq(int* n, Float x[]);
void solvde(int itmax, Float conv, Float slowc, Float scalv[],
	int indexv[], int ne, int nb, int m, Float** y, Float*** c, Float** s);
void sor(double** a, double** b, double** c, double** d, double** e,
	double** f, double** u, int jmax, double rjac);
void sort(unsigned long n, Float arr[]);
void sort2(unsigned long n, Float arr[], Float brr[]);
void sort3(unsigned long n, Float ra[], Float rb[], Float rc[]);
void spctrm(FILE* fp, Float p[], int m, int k, int ovrlap);
void spear(Float data1[], Float data2[], unsigned long n, Float* d, Float* zd,
	Float* probd, Float* rs, Float* probrs);
void sphbes(int n, Float x, Float* sj, Float* sy, Float* sjp, Float* syp);
void splie2(Float x1a[], Float x2a[], Float** ya, int m, int n, Float** y2a);
void splin2(Float x1a[], Float x2a[], Float** ya, Float** y2a, int m, int n,
	Float x1, Float x2, Float* y);
void spline(Float x[], Float y[], int n, Float yp1, Float ypn, Float y2[]);
void splint(Float xa[], Float ya[], Float y2a[], int n, Float x, Float* y);
void spread(Float y, Float yy[], unsigned long n, Float x, int m);
void sprsax(Float sa[], unsigned long ija[], Float x[], Float b[],
	unsigned long n);
void sprsin(Float** a, int n, Float thresh, unsigned long nmax, Float sa[],
	unsigned long ija[]);
void sprspm(Float sa[], unsigned long ija[], Float sb[], unsigned long ijb[],
	Float sc[], unsigned long ijc[]);
void sprstm(Float sa[], unsigned long ija[], Float sb[], unsigned long ijb[],
	Float thresh, unsigned long nmax, Float sc[], unsigned long ijc[]);
void sprstp(Float sa[], unsigned long ija[], Float sb[], unsigned long ijb[]);
void sprstx(Float sa[], unsigned long ija[], Float x[], Float b[],
	unsigned long n);
void stifbs(Float y[], Float dydx[], int nv, Float* xx,
	Float htry, Float eps, Float yscal[], Float* hdid, Float* hnext,
	void (*derivs)(Float, Float[], Float[]));
void stiff(Float y[], Float dydx[], int n, Float* x,
	Float htry, Float eps, Float yscal[], Float* hdid, Float* hnext,
	void (*derivs)(Float, Float[], Float[]));
void stoerm(Float y[], Float d2y[], int nv, Float xs,
	Float htot, int nstep, Float yout[],
	void (*derivs)(Float, Float[], Float[]));
void svbksb(Float** u, Float w[], Float** v, int m, int n, Float b[],
	Float x[]);
void svdcmp(Float** a, int m, int n, Float w[], Float** v);
void svdfit(Float x[], Float y[], Float sig[], int ndata, Float a[],
	int ma, Float** u, Float** v, Float w[], Float* chisq,
	void (*funcs)(Float, Float[], int));
void svdvar(Float** v, int ma, Float w[], Float** cvm);
void toeplz(Float r[], Float x[], Float y[], int n);
void tptest(Float data1[], Float data2[], unsigned long n, Float* t, Float* prob);
void tqli(Float d[], Float e[], int n, Float** z);
Float trapzd(Float(*func)(Float), Float a, Float b, int n);
void tred2(Float** a, int n, Float d[], Float e[]);
void tridag(Float a[], Float b[], Float c[], Float r[], Float u[],
	unsigned long n);
Float trncst(Float x[], Float y[], int iorder[], int ncity, int n[]);
void trnspt(int iorder[], int ncity, int n[]);
void ttest(Float data1[], unsigned long n1, Float data2[], unsigned long n2,
	Float* t, Float* prob);
void tutest(Float data1[], unsigned long n1, Float data2[], unsigned long n2,
	Float* t, Float* prob);
void twofft(Float data1[], Float data2[], Float fft1[], Float fft2[],
	unsigned long n);
void vander(double x[], double w[], double q[], int n);
void vegas(Float regn[], int ndim, Float(*fxn)(Float[], Float), int init,
	unsigned long ncall, int itmx, int nprn, Float* tgral, Float* sd,
	Float* chi2a);
void voltra(int n, int m, Float t0, Float h, Float* t, Float** f,
	Float(*g)(int, Float), Float(*ak)(int, int, Float, Float));
void wt1(Float a[], unsigned long n, int isign,
	void (*wtstep)(Float[], unsigned long, int));
void wtn(Float a[], unsigned long nn[], int ndim, int isign,
	void (*wtstep)(Float[], unsigned long, int));
void wwghts(Float wghts[], int n, Float h,
	void (*kermom)(double[], double, int));
int zbrac(Float(*func)(Float), Float* x1, Float* x2);
void zbrak(Float(*fx)(Float), Float x1, Float x2, int n, Float xb1[],
	Float xb2[], int* nb);
Float zbrent(Float(*func)(Float), Float x1, Float x2, Float tol);
void zrhqr(Float a[], int m, Float rtr[], Float rti[]);
Float zriddr(Float(*func)(Float), Float x1, Float x2, Float xacc);
void zroots(fcomplex a[], int m, fcomplex roots[], int polish);

#endif /* _NR_H_ */
