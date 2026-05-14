/*
 optimizer.cpp

 Feb 8, 2020 - created by Jesse Wang

 line search using gradient

 */

#include "pch.h"
#include "framework.h"
#include "optimizer.h"

CYNOMATHUTILITY_API bool __stdcall cyno_lineSearch(long func,  int numOfDim, double *xold, double* fold, double* g, double* p, double* x, double* f, double max_step){

//Feb 9, 2020 - created by Jesse Wang - linear search with gradient decendent
//Feb 9, 2020 -newly crated not tested yet
// func -is target function for which roots are desired
//numofDim - dimension of line search
//xold- intinital guess of the solution (vector)
//fold- guess of the target function solution
//g - gradient of the function at xold (vector)
//p newton direction at xold (vector)
//x solution of the target function (vector)
//f value of the target function at x
//stpmax - maximum steps

// Declare the function pointer
	FUNCPTRARR vbFunc;
	vbFunc = (FUNCPTRARR)func;

	int n = numOfDim;

	double sum = 0;
	for (int i = 0; i < n; i++) sum += p[i] * p[i];
	sum = sqrt(sum);

	if (sum > max_step) for (int i = 0; i < n; i++) p[i] *= max_step/sum;

	double slope = 0;
	for (int i = 0; i < n; i++) slope += g[i] * p[i];
	

	double test = 0;
	for (int i = 0; i < n; i++) {
		double xold_abs = fabs(xold[i]);
		double temp = fabs(p[i]) / ((xold_abs > 1.0) ? xold_abs : 1.0);
		if (temp > test) test = temp;
	}


	double alamin = LSTOLX;
	if (test > 0) alamin = LSTOLX / test;

	double alam = 1, f2 = 0, alam2 = 0, tmplam;

	while (true) {
		for (int i = 0; i < n; i++) x[i] = xold[i] + alam * p[i];
		*f = vbFunc(x);
		if (alam < alamin) {
			for (int i = 0; i < n; i++) xold[i] = x[i];
			return true;
		}
		else if (*f <= *fold + LSALF * alam * slope) return true;
		else {
			if (alam == 1) tmplam = -slope / (2 * (*f - *fold - slope));
			else {
				double rhs1 = *f - *fold - alam * slope;
				double rhs2 = f2 - *fold - alam2 * slope;
				double a = (rhs1 / (alam * alam) - rhs2 / (alam2 * alam2)) / (alam - alam2);
				double b = (-alam2 * rhs1 / (alam * alam) + alam * rhs2 / (alam2 * alam2)) / (alam - alam2);
				if (a == 0) tmplam = -slope / (2 * b);
				else {
					double disc = b * b - 3 * a * slope;
					if (disc < 0) tmplam = 0.5 * alam;
					else if (b <= 0) tmplam = (-b + sqrt(disc)) / (3 * a);
					else tmplam = -slope / (b + sqrt(disc));
				}
				if (tmplam > 0.5 * alam) tmplam = 0.5 * alam;
			}
		}
		alam2 = alam;
		f2 = *f;
		alam = (tmplam > 0.1 * alam) ? tmplam : (0.1 * alam);
	}
}

