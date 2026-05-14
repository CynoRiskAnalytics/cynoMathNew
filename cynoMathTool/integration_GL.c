#include "pch.h"
#include "framework.h"
#include "math.h"
#include "specialFunctions.h"

//Purpose
//		Half Infinite Domain - Gaussian Laguerre
//		adapted from Numeric Receipe 3 gaulag
//Parameters
//		x	out array
//		w	out array
//		inputAlpha
//		npoints	input number of points
//Return Values
//		none
//Remarks
//		1) Adapted from Numeric Receipe 3 gaulag
//		2) arrays x and w must be sized properly in caller
//		3) Since it was adapted from C++ directly,
//		if it is called in VBA, lower bound of array x and array w don't need to be 0.

CYNOMATHUTILITY_API void __stdcall cyno_IntegrationGL(double x[], double w[], double inputAlpha, short npoints)
{
	const int MAXIT = 10;
	const double EPS = 1.0e-14;
	long i, its, j;
	double ai, p1, p2, p3, pp, z, z1;
	short n;

	n = npoints;

	for (i = 0; i < n; i++) {
		if (i == 0) {
			z = (1.0 + inputAlpha) * (3.0 + 0.92 * inputAlpha) / (1.0 + 2.4 * n + 1.8 * inputAlpha);
		}
		else if (i == 1) {
			z += (15.0 + 6.25 * inputAlpha) / (1.0 + 0.9 * inputAlpha + 2.5 * n);
		}
		else {
			ai = i - 1.0;
			z += ((1.0 + 2.55 * ai) / (1.9 * ai) + 1.26 * ai * inputAlpha /
				(1.0 + 3.5 * ai)) * (z - x[i - 2]) / (1.0 + 0.3 * inputAlpha);
		}
		for (its = 0; its < MAXIT; its++) {
			p1 = 1.0;
			p2 = 0.0;
			for (j = 0; j < n; j++) {
				p3 = p2;
				p2 = p1;
				p1 = ((2.0 * j + 1 + inputAlpha - z) * p2 - (j + inputAlpha) * p3) / (j + 1.0);
			}
			pp = (n * p1 - (n + inputAlpha) * p2) / z;
			z1 = z;
			z = z1 - p1 / pp;
			if (fabs(z - z1) <= EPS) break;
		}
		if (its >= MAXIT) return;
		x[i] = z;
		w[i] = -exp(gammln(inputAlpha + n) - gammln((double)n)) / (pp * n * p2);
	}
}
