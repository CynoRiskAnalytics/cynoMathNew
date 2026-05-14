#include "pch.h"
#include "framework.h"
#include "math.h"


double D1mach(short i)
{

	static double D[6];
	static short iSetArray;

	if (iSetArray != 1) {
		D[1] = 2.2250738585072E-308;
		D[2] = 1.79769313486227E+308;
		D[3] = 1.11022302462516E-16;
		D[4] = 2.22044404925031E-16;
		D[5] = 0.301029995663981;
		iSetArray = 1;
	}
	return D[i];
}



CYNOMATHUTILITY_API void __stdcall cyno_Root(FUNCPTR vbFunc, double* x1, double* x2, double xGuess, double RelError, double ABSERROR, short* iFlag)
{
	//
	// Note: in this VBA version F is passed as an untyped object, so will be subject
	//   to late binding, with reduced efficiency.  C or VB programmers should
	//   change this to ensure early binding (eg, in VB, by using the Interface
	//   instruction
	//

	//***BEGIN PROLOGUE  DFZERO
	//***PURPOSE  Search for a zero of a function F.eval(X) in a given interval
	//            (x1,x2).  It is designed primarily for problems where F.eval(x1)
	//            and F.eval(x2) have opposite signs.
	//***LIBRARY   SLATEC
	//***CATEGORY  F1B
	//***TYPE      DOUBLE PRECISION (FZERO-S, DFZERO-D)
	//***KEYWORDS  BISECTION, NONLINEAR, ROOTS, ZEROS
	//***AUTHOR  Shampine, L. F., (SNLA)
	//           Watts, H. A., (SNLA)
	//***DESCRIPTION
	//
	//     DFZERO searches for a zero of a DOUBLE PRECISION function F.eval(X)
	//     between the given DOUBLE PRECISION values x1 and x2 until the width
	//     of the interval (x1,x2) has collapsed to within a tolerance
	//     specified by the stopping criterion,
	//        ABS(x1-x2) <= 2.*(RW*ABS(x1)+AE).
	//     The method used is an efficient combination of bisection and the
	//     secant rule and is due to T. J. Dekker.
	//
	//     Description Of Arguments
	//
	//   F     :EXT   - Name of the DOUBLE PRECISION external function.  This
	//                  name must be in an EXTERNAL statement in the calling
	//                  program.  F must be a function of one DOUBLE
	//                  PRECISION argument.
	//
	//   x1     :INOUT - One end of the DOUBLE PRECISION interval (x1,x2).  The
	//                  value returned for x1 usually is the better
	//                  approximation to a zero of F.
	//
	//   x2     :INOUT - The other end of the DOUBLE PRECISION interval (x1,x2)
	//
	//   xGuess     :IN    - A (better) DOUBLE PRECISION guess of a zero of F
	//                  which could help in speeding up convergence.  If F(x1)
	//                  and F(xGuess) have opposite signs, a root will be found in
	//                  the interval (x1,xGuess);  if not, but F(xGuess) and F(x2) have
	//                  opposite signs, a root will be found in the interval
	//                  (xGuess,x2);  otherwise, the interval (x1,x2) will be
	//                  searched for a possible root.  When no better guess
	//                  is known, it is recommended that xGuess be set to x1 or x2,
	//                  since if xGuess is not interior to the interval (x1,x2), it
	//                  will be ignored.
	//
	//   RelError    :IN    - Relative error used for RW in the stopping criterion.
	//                  If the requested RelError is less than machine precision,
	//                  then RW is set to approximately machine precision.
	//
	//   AbsError    :IN    - Absolute error used in the stopping criterion.  If
	//                  the given interval (x1,x2) contains the origin, then a
	//                  nonzero value should be chosen for AbsError.
	//
	//   IFLAG :OUT   - A status code.  User must check IFLAG after each
	//                  call.  Control returns to the user from DFZERO in all
	//                  cases.
	//
	//                1  x1 is within the requested tolerance of a zero.
	//                   The interval (x1,x2) collapsed to the requested
	//                   tolerance, the function changes sign in (x1,x2), and
	//                   F(X) decreased in magnitude as (x1,x2) collapsed.
	//
	//                2  F(x1) = 0.  However, the interval (x1,x2) may not have
	//                   collapsed to the requested tolerance.
	//
	//                3  x1 may be near a singular point of F(X).
	//                   The interval (x1,x2) collapsed to the requested tol-
	//                   erance and the function changes sign in (x1,x2), but
	//                   F(X) increased in magnitude as (x1,x2) collapsed, i.e.
	//                     ABS(F(x1 out)) > MAX(ABS(F(x1 in)),ABS(F(x2 in)))
	//
	//                4  No change in sign of F(X) was found although the
	//                   interval (x1,x2) collapsed to the requested tolerance.
	//                   The user must examine this case and decide whether
	//                   x1 is near a local minimum of F(X), or x1 is near a
	//                   zero of even multiplicity, or neither of these.
	//
	//                5  Too many (> 500) function evaluations used.
	//
	//***REFERENCES  L. F. Shampine and H. A. Watts, FZERO, a root-solving
	//                 code, Report SC-TM-70-631, Sandia Laboratories,
	//                 September 1970.
	//               T. J. Dekker, Finding a zero by means of successive
	//                 linear interpolation, Constructive Aspects of the
	//                 Fundamental Theorem of Algebra, edited by x1. Dejon
	//                 and P. Henrici, Wiley-Interscience, 1969.
	//***ROUTINES CALLED  D1MACH
	//***REVISION HISTORY  (YYMMDD)
	//   700901  DATE WRITTEN
	//   890531  Changed all specific intrinsics to generic.  (WRB)
	//   890531  REVISION DATE from Version 3.2
	//   891214  Prologue converted to Version 4.0 format.  (BAB)
	//   920501  Reformatted the REFERENCES section.  (WRB)
	//***END PROLOGUE  DFZERO
	//
	//   F is now an object with a method .eval(x)
	//   F is now a function pointer
	//*****************************

	double a, ACBS, ACMB, AW, CMB;
	double ER, fa, fb, fc, fx, FZ;
	double p, q, RW, t, tol, z;
	short IC, KOUNT;


	// Declare the function pointer
	//FUNCPTR vbFunc;

	// Point the function pointer at the passed-in address.
	//vbFunc = f; //(FUNCPTR)f;

	//chaged to use pointer directly


	ER = 2 * D1mach(4);
	//
	//   Initialize.
	//
		//z = xGuess
	z = xGuess;
	//If (xGuess < x1 And xGuess < x2) Or (xGuess > x1 And xGuess > x2) Then z = x2
	if ((xGuess < *x1 && xGuess < *x2) || (xGuess > * x1&& xGuess > * x2))
		z = *x2;
	//RW = RelError: If RW < ER Then RW = ER
	//AW = ABSERROR: If AW < 0# Then AW = 0#
	RW = RelError;
	if (RW < ER)
		RW = ER;
	AW = ABSERROR;
	if (AW < 0)
		AW = 0;



	IC = 0;
	t = z;
	FZ = vbFunc(t);
	fc = FZ;
	t = *x1;
	fb = vbFunc(t);
	KOUNT = 2;



	if (((FZ > 0) - (FZ < 0)) == ((fb > 0) - (fb < 0)))
		goto l1;
	*x2 = z;
	goto l2;
l1:
	if (z == *x2)
		goto l2;
	t = *x2;
	fc = vbFunc(t);
	KOUNT = 3;
	if (((FZ > 0) - (FZ < 0)) == ((fc > 0) - (fc < 0)))
		goto l2;
	*x1 = z;
	fb = FZ;
l2:
	a = *x2;
	fa = fc;
	ACBS = fabs(*x1 - *x2);
	fx = fabs(fb);
	if (fx < fabs(fc))
		fx = fabs(fc);


l3:
	if (fabs(fc) >= fabs(fb))
		goto l4;
	a = *x1;
	fa = fb;
	*x1 = *x2;
	fb = fc;
	*x2 = a;
	fc = fa;
l4:
	CMB = 0.5 * (*x2 - *x1);
	ACMB = fabs(CMB);
	tol = RW * fabs(*x1) + AW;


	if (ACMB <= tol)
		goto l10;
	if (fb == (double)0)
		goto l11;
	if (KOUNT >= 500)
		goto l14;


	p = (*x1 - a) * fb;
	q = fa - fb;
	if (p >= (double)0)
		goto l5;
	p = -p;
	q = -q;

	
l5:
	a = *x1;
	fa = fb;
	IC = IC + 1;
	if (IC < 4)
		goto l6;
	if ((double)8 * ACMB >= ACBS)
		goto l8;
	IC = 0;
	ACBS = ACMB;
	
l6:
	if (p > fabs(q)* tol)
		goto l7;
	
	if (CMB < 0)
		*x1 = *x1 - fabs(tol);
	else
		*x1 = *x1 + fabs(tol);
	goto l9;
	
l7:
	if (p >= CMB * q)
		goto l8;
	
	*x1 = *x1 + p / q;
	goto l9;

l8:
	*x1 = *x1 + CMB;
	
l9:
	t = *x1;
	fb = vbFunc(t);
	KOUNT = KOUNT + 1;
	
	if (((fb > 0) - (fb < 0)) != ((fc > 0) - (fc < 0)))
		goto l3;
	*x2 = a;
	fc = fa;
	goto l3;


l10:
	if (((fb > 0) - (fb < 0)) == ((fc > 0) - (fc < 0)))
		goto l13;
	if (fabs(fb) > fx)
		goto l12;
	*iFlag = 1;
	return;
l11:
	*iFlag = 2;
	return;
l12:
	*iFlag = 3;
	return;
l13:
	*iFlag = 4;
	return;
l14:
	*iFlag = 5;
	return;
}



