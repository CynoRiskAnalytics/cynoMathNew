#include "pch.h"
#include "framework.h"
#include "math.h"
#include "specialFunctions.h"

const double pi = 3.141592653589793;
//Purpose
//		Standard normal inverse function
//Parameters
//		x	input probability
//Return Values
//		double
//Remarks
//	Converted from VBA function	
//Function InverseCDF(ByVal p As Double) As Double
//  Adapted for Microsoft Visual Basic from Peter Acklam's
//  "An algorithm for computing the inverse normal cumulative distribution function"
//  (http://home.online.no/~pjacklam/notes/invnorm/)
//  by John Herrero (3-Jan-03)
	

//Function InverseCDF(ByVal p As Double) As Double
CYNOMATHUTILITY_API double __stdcall cyno_StdNormInv(double p)
{

	//  Adapted for Microsoft Visual Basic from Peter Acklam's
	//  "An algorithm for computing the inverse normal cumulative distribution function"
	//  (http://home.online.no/~pjacklam/notes/invnorm/)
	//  by John Herrero (3-Jan-03)


	const double  a1 = -39.6968302866538;
	const double  a2 = 220.946098424521;
	const double  a3 = -275.928510446969;
	const double  a4 = 138.357751867269;
	const double  a5 = -30.6647980661472;
	const double  a6 = 2.50662827745924;

	const double  b1 = -54.4760987982241;
	const double  b2 = 161.585836858041;
	const double  b3 = -155.698979859887;
	const double  b4 = 66.8013118877197;
	const double  b5 = -13.2806815528857;

	const double  c1 = -7.78489400243029E-03;
	const double  c2 = -0.322396458041136;
	const double  c3 = -2.40075827716184;
	const double  c4 = -2.54973253934373;
	const double  c5 = 4.37466414146497;
	const double  c6 = 2.93816398269878;

	const double  d1 = 7.78469570904146E-03;
	const double  d2 = 0.32246712907004;
	const double  d3 = 2.445134137143;
	const double  d4 = 3.75440866190742;



	const double p_low = 0.02425;
	const double p_high = 1 - p_low;

	//Define work variables
	//Dim q As Double, r As Double
	double q, r, InverseCDF;

	//If argument out of bounds, raise error
	//If p <= 0 Or p >= 1 Then err.Raise 5
	if (isnan(p))
		return NAN;
	if (p == 0.0)
		return -INFINITY;
	if (p == 1.0)
		return INFINITY;
	if (p < 0.0 || p > 1.0)
		return NAN;


	if (p < p_low) {
		//Rational approximation for lower region
		q = sqrt(-2 * log(p));
		InverseCDF = (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) /
			((((d1 * q + d2) * q + d3) * q + d4) * q + 1);
	}
	else if (p <= p_high) {
		//Rational approximation for lower region
		q = p - 0.5;
		r = q * q;
		InverseCDF = (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * q /
			(((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + 1);
	}
	else if (p < 1) {
		//Rational approximation for upper region
		q = sqrt(-2 * log(1 - p));
		InverseCDF = -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) /
			((((d1 * q + d2) * q + d3) * q + d4) * q + 1);
	}

	return InverseCDF;
	//End Function
}

//Purpose
//		Regular normal inverse function
//Parameters
//		p	input	probability
//		xMean	input	Mean of normal distribution
//		xStdev	input	Standard deviation of normal distribution
//Return Values
//		double
//Remarks

CYNOMATHUTILITY_API double __stdcall cyno_NormInv(double p, double xMean, double xStdev) {
	if (isnan(p) || isnan(xMean) || isnan(xStdev) || xStdev <= 0.0)
		return NAN;
	return xMean + cyno_StdNormInv(p) * xStdev;
}

//Purpose
//		Standard normal CDF
//Parameters
//		x	input	z value
//Return Values
//		double
//Remarks
//		Converted from VBA function
//Function cumnormhart(x As Double) As Double
//computing the inverse normal cumulative distribution
//function where the relative error has an absolute value
//less than 1.15 � 10-9 in the entire region
//by Peter John Acklam
//http://home.online.no/~pjacklam/notes/invnorm/
//Change history

CYNOMATHUTILITY_API double __stdcall cyno_StdNormCDF(double x)
{

	//computing the inverse normal cumulative distribution
	//function where the relative error has an absolute value
	//less than 1.15 � 10-9 in the entire region
	//by Peter John Acklam
	//http://home.online.no/~pjacklam/notes/invnorm/

	
	double XAbs, result, Exponential, Build;
	XAbs = fabs(x);
	if (XAbs > 37)
		result = 0;
	else {
	
		Exponential = exp(-pow(XAbs, 2) / 2);
		if (XAbs < 7.07106781186547) {
			Build = 3.52624965998911E-02 * XAbs + 0.700383064443688;
			Build = Build * XAbs + 6.37396220353165;
			Build = Build * XAbs + 33.912866078383;
			Build = Build * XAbs + 112.079291497871;
			Build = Build * XAbs + 221.213596169931;
			Build = Build * XAbs + 220.206867912376;
			result = Exponential * Build;
			Build = 8.83883476483184E-02 * XAbs + 1.75566716318264;
			Build = Build * XAbs + 16.064177579207;
			Build = Build * XAbs + 86.7807322029461;
			Build = Build * XAbs + 296.564248779674;
			Build = Build * XAbs + 637.333633378831;
			Build = Build * XAbs + 793.826512519948;
			Build = Build * XAbs + 440.413735824752;
			result = result / Build;
		}
		else {
			Build = XAbs + 0.65;
			Build = XAbs + 4 / Build;
			Build = XAbs + 3 / Build;
			Build = XAbs + 2 / Build;
			Build = XAbs + 1 / Build;
			result = Exponential / Build / 2.506628274631;
		}
	}
	
	if (x > 0)
		result = 1 - result;
	return result;
}
//End Function

//Purpose
//		Regular normal CDF
//Parameters
//		x		input	z value
//		xMean	input	Mean of normal distribution
//		xStdev	input	Standard deviation of normal distribution
//Return Values
//		double
//Remarks

CYNOMATHUTILITY_API double __stdcall cyno_NormCDF(double x, double xMean, double xStdev)
{
	if (isnan(x) || isnan(xMean) || isnan(xStdev) || xStdev <= 0.0)
		return NAN;
	return cyno_StdNormCDF((x - xMean) / xStdev);
}

//Purpose
//		Standard normal density function
//Parameters
//		x		input	z value
//Return Values
//		double
//Remarks

CYNOMATHUTILITY_API double __stdcall cyno_StdNormDen(double x)
{
	return exp(-pow(x, 2) / 2) / sqrt(2 * pi);
}

//Purpose
//		Regular normal density function
//Parameters
//		x		in	z value
//		xMean	in	Mean of normal distribution
//		xStdev	in	Standard deviation of normal distribution
//Return Values
//		double
//Remarks

CYNOMATHUTILITY_API double __stdcall cyno_NormDen(double x, double xMean, double xStdev){
	if (isnan(x) || isnan(xMean) || isnan(xStdev) || xStdev <= 0.0)
		return NAN;
	return exp(-pow(x - xMean, 2) / (2 * pow(xStdev, 2))) / sqrt(2 * pi) / xStdev;
}



CYNOMATHUTILITY_API double __stdcall cyno_BetaDen(const double a, const double b, double x){
	if (isnan(a) || isnan(b) || isnan(x) || a <= 0.0 || b <= 0.0 || x < 0.0 || x > 1.0)
		return NAN;
	return pow(x, a - 1) * pow(1 - x, b - 1) / beta(a, b);
}


CYNOMATHUTILITY_API double _stdcall cyno_BetaCDF(const double a, const double b, double x){
	if (isnan(a) || isnan(b) || isnan(x) || a <= 0.0 || b <= 0.0 || x < 0.0 || x > 1.0)
		return NAN;
	return beta_inc(a, b, x);
}

CYNOMATHUTILITY_API double _stdcall cyno_BetaInv(const double a, const double b, double p){
	if (isnan(a) || isnan(b) || isnan(p) || a <= 0.0 || b <= 0.0 || p < 0.0 || p > 1.0)
		return NAN;

	return cyno_Beta_inv(a, b, p);
}


double gammln(const double xx) {
	int j;
	double x, tmp, y, ser;
	static const double cof[14] = { 57.1562356658629235,-59.5979603554754912,
	14.1360979747417471,-0.491913816097620199,.339946499848118887e-4,
	.465236289270485756e-4,-.983744753048795646e-4,.158088703224912494e-3,
	-.210264441724104883e-3,.217439618115212643e-3,-.164318106536763890e-3,
	.844182239838527433e-4,-.261908384015814087e-4,.368991826595316234e-5 };
	if (xx <= 0) return NAN;
	y = x = xx;
	tmp = x + 5.24218750000000000;
	tmp = (x + 0.5) * log(tmp) - tmp;
	ser = 0.999999999999997092;
	for (j = 0; j < 14; j++) ser += cof[j] / ++y;
	return tmp + log(2.5066282746310005 * ser / x);
}
double factrl(const int n) {
	static double a[171];
	static bool init = true;
	if (init) {
		init = false;
		a[0] = 1.;
		for (int i = 1; i < 171; i++) a[i] = i * a[i - 1];
	}
	if (n < 0 || n > 170) return NAN;
	return a[n];
}
double factln(const int n) {
	enum { NTOP = 2000 };
	static double a[NTOP];
	static bool init = true;
	if (init) {
		init = false;
		for (int i = 0; i < NTOP; i++) a[i] = gammln(i + 1.);
	}
	if (n < 0) return NAN;
	if (n < NTOP) return a[n];
	return gammln(n + 1.);
}
double bico(const int n, const int k) {
	if (n < 0 || k<0 || k>n) return NAN;
	if (n < 171) return floor(0.5 + factrl(n) / (factrl(k) * factrl(n - k)));
	return floor(0.5 + exp(factln(n) - factln(k) - factln(n - k)));
}
double beta(const double a, const double b) {
	return exp(gammln(a) + gammln(b) - gammln(a + b));
}


double beta_cf(const double a, const double b, double x) {

	const int itmax = 150;
	const double fpmin = 2.2250738585072014e-308;
	const double eps = 2.2204460492503131e-16;

	double qab, qap, qam, c, d, h, aa, del;
	int m, m2;

	qab = a + b;
	qap = a + 1;
	qam = a - 1;
	c = 1;
	d = 1 - qab * x / qap;
	if (fabs(d) < fpmin) d = fpmin;
	d = 1 / d;
	h = d;

	for (m = 1; m <= itmax; ++m) {
		m2 = 2 * m;
		aa = m * (b - m) * x / ((qam + m2) * (a + m2));
		d = 1 + aa * d;
		if (fabs(d) < fpmin) d = fpmin;
		c = 1 + aa / c;
		if (fabs(c) < fpmin) c = fpmin;
		d = 1 / d;
		h *= d * c;
		aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
		d = 1 + aa * d;
		if (fabs(d) < fpmin) d = fpmin;
		c = 1 + aa / c;
		if (fabs(c) < fpmin) c = fpmin;
		d = 1 / d;
		del = d * c;
		h *= del;
		if (fabs(del - 1) < eps) break;
	}
	return h;
}

double _stdcall cyno_Beta_inv(const double a, const double b, double p)
{
	if (isnan(a) || isnan(b) || isnan(p) || a <= 0.0 || b <= 0.0 || p < 0.0 || p > 1.0)
		return NAN;

	if (p == 0 || p == 1) return p;
	else {
		const unsigned int MAXIT = 50;
		const double eps = 2.2204460492503131e-16;

		double  x, dx, bt, dc, am, bm, v, u, vx, nx, lb, ub;
		int j;

		x = a / (a + b), lb = 0, ub = 1;
		if (p > 0.5) {
			x = 2 * (1 - p) * x + 2 * p - 1;
		}
		else if (p < 0.5) {
			x = 2 * p * x;
		}
		dx = 1;
		bt = 1 / beta(a, b);
		dc = (a + 1) / (a + b + 2);
		am = a - 1;
		bm = b - 1;
		j = 0;
		while (fabs(dx) > 5 * eps && j < MAXIT) {
			v = pow(x, am);
			u = pow(1 - x, bm);
			dx = v * x * (1 - x) * u * bt;
			if (x < dc) dx *= beta_cf(a, b, x) / a;
			else {
				dx *= beta_cf(b, a, 1 - x) / b;
				dx = 1 - dx;
			}
			vx = dx;
			dx = (p - dx) / (v * u * bt);
			nx = x + dx;
			if (nx < lb || nx > ub) {
				if (vx > p) ub = x;
				else lb = x;
				x = (lb + ub) / 2;
				dx = (ub - lb) / 2;
			}
			else x = nx;
			++j;
		}
		return x;
	}
}
double beta_inc(const double a, const double b, double x) {
	if (isnan(a) || isnan(b) || isnan(x) || a <= 0.0 || b <= 0.0 || x < 0.0 || x > 1.0)
		return NAN;

	double bt;
	if (x == 0) return 0;
	else if (x == 1) return 1;
	else {
		bt = exp(gammln(a + b) - gammln(a) - gammln(b) + a * log(x) + b * log(1.0 - x));
		if (x < (a + 1) / (a + b + 2)) return bt * beta_cf(a, b, x) / a;
		else return 1.0 - bt * beta_cf(b, a, 1.0 - x) / b;
	}
}
