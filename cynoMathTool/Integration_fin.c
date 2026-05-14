#include "pch.h"
#include "framework.h"
#include "math.h"


#define MAXLIMIT 601
#define MAXPOINTS 301

void dqelg(short n, double EPSTAB[], double* result, double* ABSERR,
    double RES3LA[], short* NRES)
{
  
    //***SUBSIDIARY
    //***PURPOSE  The routine determines the limit of a given sequence of
    //            approximations, by means of the Epsilon algorithm of
    //            P.Wynn. An estimate of the absolute error is also given.
    //            The condensed Epsilon table is computed. Only those
    //            elements needed for the computation of the next diagonal
    //            are preserved.
    //***LIBRARY   SLATEC
    //***TYPE      DOUBLE PRECISION (QELG-S, SLATEC_DQELG-D)
    //***KEYWORDS  CONVERGENCE ACCELERATION, EPSILON ALGORITHM, EXTRAPOLATION
    //***AUTHOR  Piessens, Robert
    //             Applied Mathematics and Programming Division
    //             K. U. Leuven
    //           de Doncker, Elise
    //             Applied Mathematics and Programming Division
    //             K. U. Leuven
    //***DESCRIPTION
    //
    //           Epsilon algorithm
    //           Standard fortran subroutine
    //           Double precision version
    //
    //           PARAMETERS
    //              N      - Integer
    //                       EPSTAB(N) contains the new element in the
    //                       first column of the epsilon table.
    //
    //              EPSTAB - Double precision
    //                       Vector of dimension 52 containing the elements
    //                       of the two lower diagonals of the triangular
    //                       epsilon table. The elements are numbered
    //                       starting at the right-hand corner of the
    //                       triangle.
    //
    //              RESULT - Double precision
    //                       Resulting approximation to the integral
    //
    //              ABSERR - Double precision
    //                       Estimate of the absolute error computed from
    //                       RESULT and the 3 previous results
    //
    //              RES3LA - Double precision
    //                       Vector of dimension 3 containing the last 3
    //                       results
    //
    //              NRES   - Integer
    //                       Number of calls to the routine
    //                       (should be zero at first call)
    //
    // Feb 1, 2010 -created by Jesse from an older version
  

    double DELTA1, DELTA2, DELTA3;
    double EPMACH, EPSINF, error;
    double ERR1, ERR2, ERR3, E0, E1;
    double E1ABS, E2, E3, OFLOW;
    double res, SS, tol1, tol2, TOL3;
    double x;

    short i, ib, IB2, IE, INDX;
    short K1, K2, K3, LIMEXP;
    short NEWELM, num;
    
   
    EPMACH = D1mach(4);
    OFLOW = D1mach(2);
    *NRES = *NRES + 1;
    *ABSERR = OFLOW;
    *result = EPSTAB[n];
    if (n < 3)
        goto l100;
    LIMEXP = 50;
    EPSTAB[n + 2] = EPSTAB[n];
    NEWELM = (n - 1) / 2; //Potential bug!
    EPSTAB[n] = OFLOW;
    num = n;
    K1 = n;

  

    for (i = 1; i <= NEWELM; i++) { //40
        K2 = K1 - 1;
        K3 = K1 - 2;
        res = EPSTAB[K1 + 2];
        E0 = EPSTAB[K3];
        E1 = EPSTAB[K2];
        E2 = res;
        E1ABS = fabs(E1);
        DELTA2 = E2 - E1;
        ERR2 = fabs(DELTA2);
        tol2 = fabs(E2) * EPMACH;
        if (E1ABS > fabs(E2))
            tol2 = E1ABS * EPMACH;
        DELTA3 = E1 - E0;
        ERR3 = fabs(DELTA3);
        TOL3 = fabs(E0) * EPMACH;
        if (E1ABS > fabs(E0))
            TOL3 = E1ABS * EPMACH;
        if (ERR2 > tol2 || ERR3 > TOL3)
            goto l10;

       

        * result = res;
        *ABSERR = ERR2 + ERR3;
        // ***JUMP OUT OF DO-LOOP
        goto l100;
    l10:
        E3 = EPSTAB[K1];
        EPSTAB[K1] = E1;
        DELTA1 = E1 - E3;
        ERR1 = fabs(DELTA1);
        tol1 = fabs(E3) * EPMACH;
        if (E1ABS > fabs(E3))
            tol1 = E1ABS * EPMACH;

    
        if (ERR1 <= tol1 || ERR2 <= tol2 || ERR3 <= TOL3)
            goto l20;
        SS = (double)1 / DELTA1 + (double)1 / DELTA2 - (double)1 / DELTA3;
        EPSINF = fabs(SS * E1);
 
        if (EPSINF > 0.0001)
            goto l30;
    l20:
        n = i + i - 1;
        // ***JUMP OUT OF DO-LOOP
        goto l50;
     
    l30:
        res = E1 + (double)1 / SS;
        EPSTAB[K1] = res;
        K1 = K1 - 2;
        error = ERR2 + fabs(res - E2) + ERR3;
        if (error > * ABSERR)
            goto l40;
        *ABSERR = error;
        *result = res;
    l40:
        ;
    }
    
  

l50:
    if (n == LIMEXP)
        n = 2 * (LIMEXP / 2) - 1; 
    ib = 1;
    if ((num / 2) * 2 == num)
        ib = 2;
    IE = NEWELM + 1;
    for (i = 1; i <= IE; i++) {
        IB2 = ib + 2;
        EPSTAB[ib] = EPSTAB[IB2];
        ib = IB2;
    }
    if (num == n)
        goto l80;
    INDX = num - n + 1;
    for (i = 1; i <= n; i++) {
        EPSTAB[i] = EPSTAB[INDX];
        INDX = INDX + 1;
    }
l80:
    if (*NRES >= 4)
        goto l90;
    RES3LA[*NRES] = *result;
    *ABSERR = OFLOW;
    goto l100;

    
l90:
    *ABSERR = fabs(*result - RES3LA[3]) + fabs(*result - RES3LA[2]) + fabs(*result - RES3LA[1]);
    RES3LA[1] = RES3LA[2];
    RES3LA[2] = RES3LA[3];
    RES3LA[3] = *result;
l100:
    x = (double)5 * EPMACH * fabs(*result);
    if (x > * ABSERR)
        *ABSERR = x;
}



void dqk21(FUNCPTR vbFunc, double a, double b, double* result, double* ABSERR,
    double* RESABS, double* RESASC)
{
  
    //***PURPOSE  To compute I = Integral of F over (A,B), with error
    //                           estimate
    //                       J = Integral of ABS(F) over (A,B)
    //***LIBRARY   SLATEC (QUADPACK)
    //***CATEGORY  H2A1A2
    //***TYPE      DOUBLE PRECISION (QK21-S, SLATEC_DQK21-D)
    //***KEYWORDS  21-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE
    //***AUTHOR  Piessens, Robert
    //             Applied Mathematics and Programming Division
    //             K. U. Leuven
    //           de Doncker, Elise
    //             Applied Mathematics and Programming Division
    //             K. U. Leuven
    //***DESCRIPTION
    //
    //           Integration rules
    //           Standard fortran subroutine
    //           Double precision version
    //
    //           PARAMETERS
    //            ON ENTRY
    //              F      - Double precision
    //                       Function subprogram defining the integrand
    //                       FUNCTION F(X). The actual name for F needs to be
    //                       Declared E X T E R N A L in the driver program.
    //
    //              A      - Double precision
    //                       Lower limit of integration
    //
    //              B      - Double precision
    //                       Upper limit of integration
    //
    //            ON RETURN
    //              RESULT - Double precision
    //                       Approximation to the integral I
    //                       RESULT is computed by applying the 21-POINT
    //                       KRONROD RULE (RESK) obtained by optimal addition
    //                       of abscissae to the 10-POINT GAUSS RULE (RESG).
    //
    //              ABSERR - Double precision
    //                       Estimate of the modulus of the absolute error,
    //                       which should not exceed ABS(I-RESULT)
    //
    //              RESABS - Double precision
    //                       Approximation to the integral J
    //
    //              RESASC - Double precision
    //                       Approximation to the integral of ABS(F-I/(B-A))
    //                       over (A,B)
    //
    //***REFERENCES  (NONE)
    //***ROUTINES CALLED  SLATEC_D1MACH
  


    double ABSC, CENTR, DHLGTH, EPMACH, fc;
    double FSUM, FVAL1, FVAL2, FV1[11], FV2[11];
    double HLGTH, RESG, RESK, RESKH, UFLOW, z;
    static double WG[6], WGK[12], XGK[12];
    static short iSetArray;
    short j, JTW, JTWM1;

    // Declare the function pointer
    //FUNCPTR vbFunc;

    // Point the function pointer at the passed-in address.
    //vbFunc = (FUNCPTR)f;
    //
    //
    //           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
    //           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
    //           CORRESPONDING WEIGHTS ARE GIVEN.
    //
    //           XGK    - ABSCISSAE OF THE 21-POINT KRONROD RULE
    //                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 10-POINT
    //                    GAUSS RULE
    //                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
    //                    ADDED TO THE 10-POINT GAUSS RULE
    //
    //           WGK    - WEIGHTS OF THE 21-POINT KRONROD RULE
    //
    //           WG     - WEIGHTS OF THE 10-POINT GAUSS RULE
    //
    //
    // GAUSS QUADRATURE WEIGHTS AND KRONROD QUADRATURE ABSCISSAE AND WEIGHTS
    // AS EVALUATED WITH 80 DECIMAL DIGIT ARITHMETIC BY L. W. FULLERTON,
    // BELL LABS, NOV. 1981.
    //
    if (iSetArray != 1) {
        WG[1] = 6.66713443086881E-02;
        WG[2] = 0.149451349150581;
        WG[3] = 0.219086362515982;
        WG[4] = 0.269266719309996;
        WG[5] = 0.295524224714753;
        //
        XGK[1] = 0.995657163025808;
        XGK[2] = 0.973906528517172;
        XGK[3] = 0.930157491355708;
        XGK[4] = 0.865063366688985;
        XGK[5] = 0.780817726586417;
        XGK[6] = 0.679409568299024;
        XGK[7] = 0.562757134668605;
        XGK[8] = 0.433395394129247;
        XGK[9] = 0.29439286270146;
        XGK[10] = 0.148874338981631;
        XGK[11] = (double)0;
        //
        WGK[1] = 1.16946388673719E-02;
        WGK[2] = 3.25581623079647E-02;
        WGK[3] = 0.054755896574352;
        WGK[4] = 0.07503967481092;
        WGK[5] = 9.31254545836976E-02;
        WGK[6] = 0.109387158802298;
        WGK[7] = 0.123491976262066;
        WGK[8] = 0.134709217311473;
        WGK[9] = 0.14277593857706;
        WGK[10] = 0.147739104901338;
        WGK[11] = 0.149445554002917;

        iSetArray = 1;
    }
    
    //           LIST OF MAJOR VARIABLES
    //           -----------------------
    //
    //           CENTR  - MID POINT OF THE INTERVAL
    //           HLGTH  - HALF-LENGTH OF THE INTERVAL
    //           ABSC   - ABSCISSA
    //           FVAL*  - FUNCTION VALUE
    //           RESG   - RESULT OF THE 10-POINT GAUSS FORMULA
    //           RESK   - RESULT OF THE 21-POINT KRONROD FORMULA
    //           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B),
    //                    I.E. TO I/(B-A)
    //
    //
    //           MACHINE DEPENDENT CONSTANTS
    //           ---------------------------
    //
    //           EPMACH IS THE LARGEST RELATIVE SPACING.
    //           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
    //
    //***FIRST EXECUTABLE STATEMENT  SLATEC_DQK21
    //    EPMACH = Slatec_D1mach(4)
    //    UFLOW = Slatec_D1mach(1)
    ////
    //    CENTR = 0.5 * (a + b)
    //    HLGTH = 0.5 * (b - a)
    //    DHLGTH = Abs(HLGTH)

    EPMACH = D1mach(4);
    UFLOW = D1mach(1);
    //
    CENTR = 0.5 * (a + b);
    HLGTH = 0.5 * (b - a);
    DHLGTH = fabs(HLGTH);


    RESG = (double)0;
    fc = vbFunc(CENTR);
    RESK = WGK[11] * fc;
    *RESABS = fabs(RESK);
    for (j = 1; j <= 5; j++) {
        JTW = 2 * j;
        ABSC = HLGTH * XGK[JTW];
        FVAL1 = vbFunc(CENTR - ABSC);
        FVAL2 = vbFunc(CENTR + ABSC);
        FV1[JTW] = FVAL1;
        FV2[JTW] = FVAL2;
        FSUM = FVAL1 + FVAL2;
        RESG = RESG + WG[j] * FSUM;
        RESK = RESK + WGK[JTW] * FSUM;
        *RESABS = *RESABS + WGK[JTW] * (fabs(FVAL1) + fabs(FVAL2));
    }



    for (j = 1; j <= 5; j++) {
        JTWM1 = 2 * j - 1;
        ABSC = HLGTH * XGK[JTWM1];
        FVAL1 = vbFunc(CENTR - ABSC);
        FVAL2 = vbFunc(CENTR + ABSC);
        FV1[JTWM1] = FVAL1;
        FV2[JTWM1] = FVAL2;
        FSUM = FVAL1 + FVAL2;
        RESK = RESK + WGK[JTWM1] * FSUM;
        *RESABS = *RESABS + WGK[JTWM1] * (fabs(FVAL1) + fabs(FVAL2));
    }

 
    RESKH = RESK * 0.5;
    *RESASC = WGK[11] * fabs(fc - RESKH);
    for (j = 1; j <= 10; j++)
        *RESASC = *RESASC + WGK[j] * (fabs(FV1[j] - RESKH) + fabs(FV2[j] - RESKH));
    *result = RESK * HLGTH;
    *RESABS = *RESABS * DHLGTH;
    *RESASC = *RESASC * DHLGTH;
    *ABSERR = fabs((RESK - RESG) * HLGTH);
    if (*RESASC != (double)0 && *ABSERR != (double)0) {
        z = pow(((double)200 * (*ABSERR) / (*RESASC)), 1.5);
        if (z > (double)1)
            *ABSERR = *RESASC;
        else
            *ABSERR = *RESASC * z;
    }
    if (*RESABS > UFLOW / ((double)50 * EPMACH)) {
        z = (EPMACH * (double)50) * *RESABS;
        if (z > * ABSERR)
            *ABSERR = z;
    }

}







void dqpsrt(short LIMIT, short LAST, short* MAXERR, double* ERMAX,
    double ELIST[], short IORD[], short* NRMAX)
{
    //***BEGIN PROLOGUE  SLATEC_DQPSRT
    //***SUBSIDIARY
    //***PURPOSE  This routine maintains the descending ordering in the
    //            list of the local error estimated resulting from the
    //            interval subdivision process. At each call two error
    //            estimates are inserted using the sequential search
    //            method, top-down for the largest error estimate and
    //            bottom-up for the smallest error estimate.
    //***TYPE      DOUBLE PRECISION (QPSRT-S, SLATEC_DQPSRT-D)
    //***KEYWORDS  SEQUENTIAL SORTING
    //***AUTHOR  Piessens, Robert
    //             Applied Mathematics and Programming Division
    //             K. U. Leuven
    //           de Doncker, Elise
    //             Applied Mathematics and Programming Division
    //             K. U. Leuven
    //***DESCRIPTION
    //
    //           Ordering routine
    //           Standard fortran subroutine
    //           Double precision version
    //
    //           PARAMETERS (MEANING AT OUTPUT)
    //              LIMIT  - Integer
    //                       Maximum number of error estimates the list
    //                       can contain
    //
    //              LAST   - Integer
    //                       Number of error estimates currently in the list
    //
    //              MAXERR - Integer
    //                       MAXERR points to the NRMAX-th largest error
    //                       estimate currently in the list
    //
    //              ERMAX  - Double precision
    //                       NRMAX-th largest error estimate
    //                       ERMAX = ELIST(MAXERR)
    //
    //              ELIST  - Double precision
    //                       Vector of dimension LAST containing
    //                       the error estimates
    //
    //              IORD   - Integer
    //                       Vector of dimension LAST, the first K elements
    //                       of which contain pointers to the error
    //                       estimates, such that
    //                       ELIST(IORD(1)),...,  ELIST(IORD(K))
    //                       form a decreasing sequence, with
    //                       K = LAST if LAST<=(LIMIT/2+2), and
    //                       K = LIMIT+1-LAST otherwise
    //
    //              NRMAX  - Integer
    //                       MAXERR = IORD(NRMAX)
    //
    //***SEE ALSO  DQAGE, DQAGIE, DQAGPE, DQAWSE

 

    double ERRMAX, ERRMIN;
    short i, IBEG, IDO, ISUCC, j, JBND, JUPBN, k;

  

    if (LAST > 2)
        goto l10;
    IORD[1] = 1;
    IORD[2] = 2;
    goto l90;

  

l10:
    ERRMAX = ELIST[*MAXERR];
    if (*NRMAX == 1)
        goto l30;
    IDO = *NRMAX - 1;
    for (i = 1; i <= IDO; i++) {
        ISUCC = IORD[*NRMAX - 1];
        // ***JUMP OUT OF DO-LOOP
        if (ERRMAX <= ELIST[ISUCC])
            goto l30;
        IORD[*NRMAX] = ISUCC;
        *NRMAX = *NRMAX - 1;
    }

   

l30:
    JUPBN = LAST;
    if (LAST > (LIMIT / 2 + 2)) 
        JUPBN = LIMIT + 3 - LAST;
    ERRMIN = ELIST[LAST];
 

    JBND = JUPBN - 1;
    IBEG = *NRMAX + 1;
    if (IBEG > JBND)
        goto l50;
    for (i = IBEG; i <= JBND; i++) {
        ISUCC = IORD[i];
        // ***JUMP OUT OF DO-LOOP
        if (ERRMAX >= ELIST[ISUCC])
            goto l60;
        IORD[i - 1] = ISUCC;
    }
l50:
    IORD[JBND] = *MAXERR;
    IORD[JUPBN] = LAST;
    goto l90;


  

l60:
    IORD[i - 1] = *MAXERR;
    k = JBND;
    for (j = i; j <= JBND; j++) {
        ISUCC = IORD[k];
        // ***JUMP OUT OF DO-LOOP
        if (ERRMIN < ELIST[ISUCC])
            goto l80;
        IORD[k + 1] = ISUCC;
        k = k - 1;
    }
    IORD[i] = LAST;
    goto l90;
l80:
    IORD[k + 1] = LAST;

    //
    //           SET MAXERR AND ERMAX.
    //
l90:
    *MAXERR = IORD[*NRMAX];
    *ERMAX = ELIST[*MAXERR];
}




CYNOMATHUTILITY_API void __stdcall cyno_IntegrationFin(FUNCPTR vbFunc, short MinIndex, short MaxIndex, double* element, double a, double b,
    double EPSABS, double EPSREL, short LIMIT,
    double* result, double* ABSERR, short* neval, short* ier)
{
  

    //***BEGIN PROLOGUE  DQAGPE
    //***PURPOSE  Approximate a given definite integral I = Integral of F
    //            over (A,B), hopefully satisfying the accuracy claim:
    //                 ABS(I-RESULT)<=MAX(EPSABS,EPSREL*ABS(I)).
    //            Break points of the integration interval, where local
    //            difficulties of the integrand may occur (e.g. singularities
    //            or discontinuities) are provided by the user.
    //***LIBRARY   SLATEC (QUADPACK)
    //***CATEGORY  H2A2A1
    //***TYPE      DOUBLE PRECISION (QAGPE-S, DQAGPE-D)
    //***KEYWORDS  AUTOMATIC INTEGRATOR, EXTRAPOLATION, GENERAL-PURPOSE,
    //             GLOBALLY ADAPTIVE, QUADPACK, QUADRATURE,
    //             SINGULARITIES AT USER SPECIFIED POINTS
    //***AUTHOR  Piessens, Robert
    //             Applied Mathematics and Programming Division
    //             K. U. Leuven
    //           de Doncker, Elise
    //             Applied Mathematics and Programming Division
    //             K. U. Leuven
    //***DESCRIPTION
    //
    //        Computation of a definite integral
    //        Standard fortran subroutine
    //        Double precision version
    //
    //        PARAMETERS
    //         ON ENTRY
    //            F      - Double precision
    //                     Function subprogram defining the integrand
    //                     function F(X). **DJD** F needs to have
    //                   eval(x) method
    //                   Breaks() property -
    //                       Breaks.MinIndex
    //                       Breaks.MaxIndex
    //                       Breaks.element(I)
    //
    //            A      - Double precision
    //                     Lower limit of integration
    //
    //            B      - Double precision
    //                     Upper limit of integration
    //
    //            NPTS2  - Integer (removed, DJD)
    //                     Number equal to two more than the number of
    //                     user-supplied break points within the integration
    //                     range, NPTS2>=2.
    //                     If NPTS2<2, the routine will end with IER = 6.
    //
    //            POINTS - Double precision (removed, DJD)
    //                     Vector of dimension NPTS2, the first (NPTS2-2)
    //                     elements of which are the user provided break
    //                     POINTS. **DJD** Must be in ascending order.
    //
    //            EPSABS - Double precision
    //                     Absolute accuracy requested
    //            EPSREL - Double precision
    //                     Relative accuracy requested
    //                     If  EPSABS<=0
    //                     and EPSREL<MAX(50*REL.MACH.ACC.,0.5D-28),
    //                     the routine will end with IER = 6.
    //
    //            LIMIT  - Integer
    //                     Gives an upper bound on the number of subintervals
    //                     in the partition of (A,B), LIMIT>=NPTS2
    //                     If LIMIT<NPTS2, the routine will end with
    //                     IER = 6.
    //
    //         ON RETURN
    //            RESULT - Double precision
    //                     Approximation to the integral
    //
    //            ABSERR - Double precision
    //                     Estimate of the modulus of the absolute error,
    //                     which should equal or exceed ABS(I-RESULT)
    //
    //            NEVAL  - Integer
    //                     Number of integrand evaluations
    //
    //            IER    - Integer
    //                     IER = 0 Normal and reliable termination of the
    //                             routine. It is assumed that the requested
    //                             accuracy has been achieved.
    //                     IER>0 Abnormal termination of the routine.
    //                             The estimates for integral and error are
    //                             less reliable. It is assumed that the
    //                             requested accuracy has not been achieved.
    //            ERROR MESSAGES
    //                     IER = 1 Maximum number of subdivisions allowed
    //                             has been achieved. One can allow more
    //                             subdivisions by increasing the value of
    //                             LIMIT (and taking the according dimension
    //                             adjustments into account). However, if
    //                             this yields no improvement it is advised
    //                             to analyze the integrand in order to
    //                             determine the integration difficulties. If
    //                             the position of a local difficulty can be
    //                             determined (i.e. SINGULARITY,
    //                             DISCONTINUITY within the interval), it
    //                             should be supplied to the routine as an
    //                             element of the vector points. If necessary
    //                             an appropriate special-purpose integrator
    //                             must be used, which is designed for
    //                             handling the type of difficulty involved.
    //                         = 2 The occurrence of roundoff error is
    //                             detected, which prevents the requested
    //                             tolerance from being achieved.
    //                             The error may be under-estimated.
    //                         = 3 Extremely bad integrand behaviour occurs
    //                             At some points of the integration
    //                             interval.
    //                         = 4 The algorithm does not converge.
    //                             Roundoff error is detected in the
    //                             extrapolation table. It is presumed that
    //                             the requested tolerance cannot be
    //                             achieved, and that the returned result is
    //                             the best which can be obtained.
    //                         = 5 The integral is probably divergent, or
    //                             slowly convergent. It must be noted that
    //                             divergence can occur with any other value
    //                             of IER>0.
    //                         = 6 The input is invalid because
    //                             NPTS2<2 or
    //                             Break points are specified outside
    //                             the integration range or
    //                             (EPSABS<=0 and
    //                              EPSREL<MAX(50*REL.MACH.ACC.,0.5D-28))
    //                             or LIMIT<NPTS2.
    //                             RESULT, ABSERR, NEVAL, LAST, RLIST(1),
    //                             and ELIST(1) are set to zero. ALIST(1) and
    //                             BLIST(1) are set to A and B respectively.
    //
    //            ALIST  - Double precision (removed, DJD)
    //                     Vector of dimension at least LIMIT, the first
    //                      LAST  elements of which are the left end points
    //                     of the subintervals in the partition of the given
    //                     integration range (A,B)
    //
    //            BLIST  - Double precision (removed, DJD)
    //                     Vector of dimension at least LIMIT, the first
    //                      LAST  elements of which are the right end points
    //                     of the subintervals in the partition of the given
    //                     integration range (A,B)
    //
    //            RLIST  - Double precision (removed, DJD)
    //                     Vector of dimension at least LIMIT, the first
    //                      LAST  elements of which are the integral
    //                     approximations on the subintervals
    //
    //            ELIST  - Double precision (removed, DJD)
    //                     Vector of dimension at least LIMIT, the first
    //                      LAST  elements of which are the moduli of the
    //                     absolute error estimates on the subintervals
    //
    //            PTS    - Double precision (removed, DJD)
    //                     Vector of dimension at least NPTS2, containing the
    //                     integration limits and the break points of the
    //                     interval in ascending sequence.
    //
    //            LEVEL  - Integer (removed, DJD)
    //                     Vector of dimension at least LIMIT, containing the
    //                     subdivision levels of the subinterval, i.e. if
    //                     (AA,BB) is a subinterval of (P1,P2) where P1 as
    //                     well as P2 is a user-provided break point or
    //                     integration limit, then (AA,BB) has level L if
    //                     ABS(BB-AA) = ABS(P2-P1)*2**(-L).
    //
    //            NDIN   - Integer (removed, DJD)
    //                     Vector of dimension at least NPTS2, after first
    //                     integration over the intervals (PTS(I)),PTS(I+1),
    //                     I = 0,1, ..., NPTS2-2, the error estimates over
    //                     some of the intervals may have been increased
    //                     artificially, in order to put their subdivision
    //                     forward. If this happens for the subinterval
    //                     numbered K, NDIN(K) is put to 1, otherwise
    //                     NDIN(K) = 0.
    //
    //            IORD   - Integer (removed, DJD)
    //                     Vector of dimension at least LIMIT, the first K
    //                     elements of which are pointers to the
    //                     error estimates over the subintervals,
    //                     such that ELIST(IORD(1)), ..., ELIST(IORD(K))
    //                     form a decreasing sequence, with K = LAST
    //                     If LAST<=(LIMIT/2+2), and K = LIMIT+1-LAST
    //                     otherwise
    //
    //            LAST   - Integer (removed, DJD)
    //                     Number of subintervals actually produced in the
    //                     subdivisions process
   



    double ABSEPS, AREA, AREA1, AREA12;
    double AREA2, A1, A2, B1, B2;
    double CORREC, DEFABS, DEFAB1, DEFAB2;
    double DRES, EPMACH, ERLARG, ERLAST;
    double ERRBND, ERRMAX, ERROR1, ERRO12;
    double ERROR2, ERRSUM, ERTEST, OFLOW;
    double RESA, RESABS, RESEPS;
    double SIGN, UFLOW;
    double LIMLO, LIMHI, x;
    short i, ID, IERRO, IND1;
    short IND2, IROFF1, IROFF2;
    short IROFF3, j, JLOW, JUPBND;
    short k, KSGN, KTMIN, LAST;
    short LEVCUR, LEVMAX;
    short MAXERR, NINT, nPts;
    short NRES, NRMAX, NUMRL2, NPTS2;
    bool EXTRAP, NOEXT;

  

    double ALIST[MAXLIMIT], BLIST[MAXLIMIT], ELIST[MAXLIMIT];
    short IORD[MAXLIMIT], LEVEL[MAXLIMIT], NDIN[MAXPOINTS];
    double PTS[MAXPOINTS], RES3LA[4];
    double RLIST[MAXLIMIT], RLIST2[52];

    // Declare the function pointer
    //FUNCPTR vbFunc;

    // Point the function pointer at the passed-in address.
    //vbFunc = (FUNCPTR)f;

//
//            THE DIMENSION OF RLIST2 IS DETERMINED BY THE VALUE OF
//            LIMEXP IN SUBROUTINE EPSALG (RLIST2 SHOULD BE OF DIMENSION
//            (LIMEXP+2) AT LEAST).
//
//
//            LIST OF MAJOR VARIABLES
//            -----------------------
//
//           ALIST     - LIST OF LEFT END POINTS OF ALL SUBINTERVALS
//                       CONSIDERED UP TO NOW
//           BLIST     - LIST OF RIGHT END POINTS OF ALL SUBINTERVALS
//                       CONSIDERED UP TO NOW
//           RLIST(I)  - APPROXIMATION TO THE INTEGRAL OVER
//                       (ALIST(I),BLIST(I))
//           RLIST2    - ARRAY OF DIMENSION AT LEAST LIMEXP+2
//                       CONTAINING THE PART OF THE EPSILON TABLE WHICH
//                       IS STILL NEEDED FOR FURTHER COMPUTATIONS
//           ELIST(I)  - ERROR ESTIMATE APPLYING TO RLIST(I)
//           MAXERR    - POINTER TO THE INTERVAL WITH LARGEST ERROR
//                       ESTIMATE
//           ERRMAX    - ELIST(MAXERR)
//           ERLAST    - ERROR ON THE INTERVAL CURRENTLY SUBDIVIDED
//                       (BEFORE THAT SUBDIVISION HAS TAKEN PLACE)
//           AREA      - SUM OF THE INTEGRALS OVER THE SUBINTERVALS
//           ERRSUM    - SUM OF THE ERRORS OVER THE SUBINTERVALS
//           ERRBND    - REQUESTED ACCURACY MAX(EPSABS,EPSREL*
//                       ABS(RESULT))
//           *****1    - VARIABLE FOR THE LEFT SUBINTERVAL
//           *****2    - VARIABLE FOR THE RIGHT SUBINTERVAL
//           LAST      - INDEX FOR SUBDIVISION
//           NRES      - NUMBER OF CALLS TO THE EXTRAPOLATION ROUTINE
//           NUMRL2    - NUMBER OF ELEMENTS IN RLIST2. IF AN APPROPRIATE
//                       APPROXIMATION TO THE COMPOUNDED INTEGRAL HAS
//                       BEEN OBTAINED, IT IS PUT IN RLIST2(NUMRL2) AFTER
//                       NUMRL2 HAS BEEN INCREASED BY ONE.
//           ERLARG    - SUM OF THE ERRORS OVER THE INTERVALS LARGER
//                       THAN THE SMALLEST INTERVAL CONSIDERED UP TO NOW
//           EXTRAP    - LOGICAL VARIABLE DENOTING THAT THE ROUTINE
//                       IS ATTEMPTING TO PERFORM EXTRAPOLATION. I.E.
//                       BEFORE SUBDIVIDING THE SMALLEST INTERVAL WE
//                       TRY TO DECREASE THE VALUE OF ERLARG.
//           NOEXT     - LOGICAL VARIABLE DENOTING THAT EXTRAPOLATION IS
//                       NO LONGER ALLOWED (TRUE-VALUE)
//
//            MACHINE DEPENDENT CONSTANTS
//            ---------------------------
//
//           EPMACH IS THE LARGEST RELATIVE SPACING.
//           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
//           OFLOW IS THE LARGEST POSITIVE MAGNITUDE.
//
//***FIRST EXECUTABLE STATEMENT  DQAGPE
    //EPMACH = Slatec_D1mach(4)
    EPMACH = D1mach(4);
  
    * ier = 0;
    *neval = 0;
    LAST = 0;
    *result = (double)0;
    *ABSERR = (double)0;
    ALIST[1] = a;
    BLIST[1] = b;
    RLIST[1] = (double)0;
    ELIST[1] = (double)0;
    IORD[1] = 0;
    LEVEL[1] = 0;
    //
  
    if (a < b) {
        LIMLO = a;
        LIMHI = b;
    }
    else {
        LIMLO = b;
        LIMHI = b;
    }
    NPTS2 = 1;
    PTS[NPTS2] = LIMLO;

 
    for (i = MinIndex; i <= MaxIndex; i++)
    {
        x = element[i];
        if (x > LIMLO&& x < LIMHI) {
            if (x < PTS[NPTS2] || NPTS2 >= MAXPOINTS) {
                *ier = 6;
                goto l999;
            }
            NPTS2 = NPTS2 + 1;
            PTS[NPTS2] = x;
        }
    }


    NPTS2 = NPTS2 + 1;
    if (NPTS2 > MAXPOINTS) {
        *ier = 6;
        goto l999;
    }
    PTS[NPTS2] = LIMHI;
  
    nPts = NPTS2 - 2;
    if (NPTS2 < 2 || LIMIT <= nPts || LIMIT > MAXLIMIT || NPTS2 > MAXPOINTS
        || (EPSABS <= (double)0 && EPSREL < (double)50 * EPMACH))
        *ier = 6;
    if (*ier == 6)
        goto l999;
   

    SIGN = (double)1;
    if (a > b)
        SIGN = (double)(-1);
    NINT = nPts + 1;
    A1 = PTS[1];
    RESABS = (double)0;
    for (i = 1; i <= NINT; i++) {
        B1 = PTS[i + 1];
        dqk21(vbFunc, A1, B1, &AREA1, &ERROR1, &DEFABS, &RESA);
        *ABSERR = *ABSERR + ERROR1;
        *result = *result + AREA1;
        NDIN[i] = 0;
        if (ERROR1 == RESA && ERROR1 != (double)0)
            NDIN[i] = 1;
        RESABS = RESABS + DEFABS;
        LEVEL[i] = 0;
        ELIST[i] = ERROR1;
        ALIST[i] = A1;
        BLIST[i] = B1;
        RLIST[i] = AREA1;
        IORD[i] = i;
        A1 = B1;
    }
    ERRSUM = (double)0;
    for (i = 1; i <= NINT; i++) {
        if (NDIN[i] == 1)
            ELIST[i] = *ABSERR;
        ERRSUM = ERRSUM + ELIST[i];
    }
   

    LAST = NINT;
    *neval = 21 * NINT;
    DRES = fabs(*result);
    ERRBND = EPSREL * DRES;
    if (EPSABS > ERRBND)
        ERRBND = EPSABS;
    if (*ABSERR <= (double)100 * EPMACH * RESABS && *ABSERR > ERRBND)
        *ier = 2;
    if (NINT == 1)
        goto l80;
    for (i = 1; i <= nPts; i++) {
        JLOW = i + 1;
        IND1 = IORD[i];
        for (j = JLOW; j <= NINT; j++) {
            IND2 = IORD[j];
            if (ELIST[IND1] > ELIST[IND2])
                goto l60;
            IND1 = IND2;
            k = j;
        l60:
            ;
        }
        if (IND1 == IORD[i])
            goto l70;
        IORD[k] = IORD[i];
        IORD[i] = IND1;
    l70:
        ;
    }
    if (LIMIT < NPTS2)
        *ier = 1;
l80:
    if (*ier != 0 || *ABSERR <= ERRBND)
        goto l999;

  

    RLIST2[1] = *result;
    MAXERR = IORD[1];
    ERRMAX = ELIST[MAXERR];
    AREA = *result;
    NRMAX = 1;
    NRES = 0;
    NUMRL2 = 1;
    KTMIN = 0;
    EXTRAP = false;
    NOEXT = false;
    ERLARG = ERRSUM;
    ERTEST = ERRBND;
    LEVMAX = 1;
    IROFF1 = 0;
    IROFF2 = 0;
    IROFF3 = 0;
    IERRO = 0;
    UFLOW = D1mach(1);
    OFLOW = D1mach(2);
    *ABSERR = OFLOW;
    KSGN = -1;
    if (DRES >= ((double)1 - (double)50 * EPMACH) * RESABS)
        KSGN = 1;

 
    for (LAST = NPTS2; LAST <= LIMIT; LAST++) {
       

        LEVCUR = LEVEL[MAXERR] + 1;
        A1 = ALIST[MAXERR];
        B1 = 0.5 * (ALIST[MAXERR] + BLIST[MAXERR]);
        A2 = B1;
        B2 = BLIST[MAXERR];
        ERLAST = ERRMAX;
        dqk21(vbFunc, A1, B1, &AREA1, &ERROR1, &RESA, &DEFAB1);
        dqk21(vbFunc, A2, B2, &AREA2, &ERROR2, &RESA, &DEFAB2);

        * neval = *neval + 42;
        AREA12 = AREA1 + AREA2;
        ERRO12 = ERROR1 + ERROR2;
        ERRSUM = ERRSUM + ERRO12 - ERRMAX;
        AREA = AREA + AREA12 - RLIST[MAXERR];
        if (DEFAB1 == ERROR1 || DEFAB2 == ERROR2)
            goto l95;
        if (fabs(RLIST[MAXERR] - AREA12) > 0.00001 * fabs(AREA12)
            || ERRO12 < 0.99 * ERRMAX)
            goto l90;
        if (EXTRAP)
            IROFF2 = IROFF2 + 1;
        else
            IROFF1 = IROFF1 + 1;
    l90:
        if (LAST > 10 && ERRO12 > ERRMAX)
            IROFF3 = IROFF3 + 1;
    l95:
        LEVEL[MAXERR] = LEVCUR;
        LEVEL[LAST] = LEVCUR;
        RLIST[MAXERR] = AREA1;
        RLIST[LAST] = AREA2;
        ERRBND = EPSREL * fabs(AREA);
        if (EPSABS > ERRBND)
            ERRBND = EPSABS;


        if (IROFF1 + IROFF2 >= 10 || IROFF3 >= 20)
            *ier = 2;
        if (IROFF2 >= 5)
            IERRO = 3;
      
        if (LAST == LIMIT)
            *ier = 1;
     
        x = ((double)1 + (double)100 * EPMACH) * (fabs(A2) + (double)1000 * UFLOW);
        if (fabs(A1) < x && fabs(B2) < x)
            *ier = 4;
     

        if (ERROR2 > ERROR1)
            goto l100;
        ALIST[LAST] = A2;
        BLIST[MAXERR] = B1;
        BLIST[LAST] = B2;
        ELIST[MAXERR] = ERROR1;
        ELIST[LAST] = ERROR2;
        goto l110;
    l100:
        ALIST[MAXERR] = A2;
        ALIST[LAST] = A1;
        BLIST[LAST] = B1;
        RLIST[MAXERR] = AREA2;
        RLIST[LAST] = AREA1;
        ELIST[MAXERR] = ERROR2;
        ELIST[LAST] = ERROR1;


    l110:
        dqpsrt(LIMIT, LAST, &MAXERR, &ERRMAX, ELIST, IORD, &NRMAX);

  
        if (ERRSUM <= ERRBND)
            goto l190;
  
        if (*ier != 0)
            goto l170;
        if (NOEXT)
            goto l160;
        ERLARG = ERLARG - ERLAST;
        if (LEVCUR + 1 <= LEVMAX)
            ERLARG = ERLARG + ERRO12;
        if (EXTRAP)
            goto l120;

               
        if (LEVEL[MAXERR] + 1 <= LEVMAX)
            goto l160;
        EXTRAP = true;
        NRMAX = 2;
    l120:
        if (IERRO == 3 || ERLARG <= ERTEST)
            goto l140;


  
      

        ID = NRMAX;
        JUPBND = LAST;
        if (LAST > (2 + LIMIT / 2))
            JUPBND = LIMIT + 3 - LAST;
        for (k = ID; k <= JUPBND; k++) {
            MAXERR = IORD[NRMAX];
            ERRMAX = ELIST[MAXERR];
            // ***JUMP OUT OF DO-LOOP
            if (LEVEL[MAXERR] + 1 <= LEVMAX)
                goto l160;
            NRMAX = NRMAX + 1;
        }

   

    l140:
        NUMRL2 = NUMRL2 + 1;
        RLIST2[NUMRL2] = AREA;
        if (NUMRL2 <= 2)
            goto l155;

        dqelg(NUMRL2, RLIST2, &RESEPS, &ABSEPS, RES3LA, &NRES);
        KTMIN = KTMIN + 1;
        if (KTMIN > 5 && *ABSERR < 0.001 * ERRSUM)
            *ier = 5;
        if (ABSEPS >= *ABSERR)
            goto l150;
        KTMIN = 0;
        *ABSERR = ABSEPS;
        *result = RESEPS;
        CORREC = ERLARG;
        ERTEST = EPSREL * fabs(RESEPS);
        if (EPSABS > ERTEST)
            ERTEST = EPSABS;
        // ***JUMP OUT OF DO-LOOP
                //If (ABSERR < ERTEST) Then GoTo 170
        if (*ABSERR < ERTEST)
            goto l170;

 

    l150:
        if (NUMRL2 == 1)
            NOEXT = true;
        if (*ier >= 5)
            goto l170;
    l155:
        MAXERR = IORD[1];
        ERRMAX = ELIST[MAXERR];
        NRMAX = 1;
        EXTRAP = false;
        LEVMAX = LEVMAX + 1;
        ERLARG = ERRSUM;

    l160:
        ;
    }
   

l170:
    if (*ABSERR == OFLOW)
        goto l190;
    if ((*ier + IERRO) == 0)
        goto l180;
    if (IERRO == 3)
        *ABSERR = *ABSERR + CORREC;
    if (*ier == 0)
        *ier = 3;
    if (*result != (double)0 && AREA != (double)0)
        goto l175;
    if (*ABSERR > ERRSUM)
        goto l190;
    if (AREA == (double)0)
        goto l210;
    goto l180;
l175:
    if (*ABSERR / fabs(*result) > ERRSUM / fabs(AREA))
        goto l190;
   

l180:
    if (KSGN == (-1) && fabs(*result) < 0.01 * DEFABS && fabs(AREA) <= 0.01 * DEFABS)
        goto l210;
    if (0.01 > (*result / AREA) || (*result / AREA) > (double)100 || ERRSUM > fabs(AREA))
        *ier = 6;
    goto l210;

  

l190:
    *result = (double)0;
    for (k = 1; k <= LAST; k++)
        *result = *result + RLIST[k];
    *ABSERR = ERRSUM;
l210:
    if (*ier > 2)
        *ier = *ier - 1;
    *result = *result * SIGN;
l999:
    return;
}






