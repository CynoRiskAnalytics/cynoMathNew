#pragma once
#define NOMINMAX

#ifdef __cplusplus
#include <cstring>
#include <cstdlib>
#else
#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#endif

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN             // Exclude rarely-used stuff from Windows headers
#include <windows.h>
#include "OAIdl.h"
#else
#ifndef __stdcall
#define __stdcall
#endif

#ifndef _stdcall
#define _stdcall
#endif
#endif

#ifdef _WIN32
#ifdef CYNOMATHUTILITY_EXPORTS
#define CYNOMATHUTILITY_API __declspec(dllimport)
#else
#define CYNOMATHUTILITY_API __declspec(dllexport)
#endif
#else
#define CYNOMATHUTILITY_API __attribute__((visibility("default")))
#endif

#ifdef __cplusplus
extern "C" {
#endif


typedef double(__stdcall* FUNCPTR)(double);
typedef double(__stdcall* FUNCPTRARR)(double *);
double D1mach(short i);


CYNOMATHUTILITY_API int __stdcall cyno_Test(double x, FUNCPTR f);
CYNOMATHUTILITY_API int __stdcall cyno_Sort(int* index, double* x, int n);


CYNOMATHUTILITY_API double __stdcall cyno_StdNormInv(double p);
CYNOMATHUTILITY_API double __stdcall cyno_NormInv(double p, double xMean, double xStdev);
CYNOMATHUTILITY_API double __stdcall cyno_StdNormCDF(double x);
CYNOMATHUTILITY_API double __stdcall cyno_NormCDF(double x, double xMean, double xStdev);
CYNOMATHUTILITY_API double __stdcall cyno_StdNormDen(double x);
CYNOMATHUTILITY_API double __stdcall cyno_NormDen(double x, double xMean, double xStdev);
CYNOMATHUTILITY_API double __stdcall cyno_BetaDen(const double a, const double b, double x);
CYNOMATHUTILITY_API double _stdcall cyno_BetaCDF(const double a, const double b, double x);
CYNOMATHUTILITY_API double _stdcall cyno_BetaInv(const double a, const double b, double p);


CYNOMATHUTILITY_API void __stdcall cyno_Root(FUNCPTR f, double* x1, double* x2, double xGuess, double RelError, double ABSERROR, short* iFlag);
CYNOMATHUTILITY_API bool __stdcall cyno_lineSearch(long func, int numOfDim, double* xold, double* fold, double* g, double* p, double* x, double* f, double max_step);


CYNOMATHUTILITY_API void __stdcall cyno_IntegrationFin(FUNCPTR f, short MinIndex, short MaxIndex, double* element, double a, double b,
    double EPSABS, double EPSREL, short LIMIT,
    double* result, double* ABSERR, short* neval, short* ier);
CYNOMATHUTILITY_API void __stdcall cyno_IntegrationGHBase(double x[], double w[], short npoints);
CYNOMATHUTILITY_API void __stdcall cyno_IntegrationGH(double x[], double w[], short npoints);

CYNOMATHUTILITY_API void __stdcall cyno_Rand_Seed(long seed);
CYNOMATHUTILITY_API double __stdcall cyno_Rand(void);
CYNOMATHUTILITY_API void __stdcall cyno_RandMT_Seed(long seed);
CYNOMATHUTILITY_API double __stdcall cyno_MTRand(void);
CYNOMATHUTILITY_API double __stdcall cyno_NormRand(void);
CYNOMATHUTILITY_API double __stdcall cyno_NormRandInv(void);
CYNOMATHUTILITY_API bool __stdcall cyno_Rand1DArrayPlainAntithetic(double* outArray, long size);
CYNOMATHUTILITY_API bool _stdcall cyno_NormalRand1DArray(double* outArray, int size);
CYNOMATHUTILITY_API bool __stdcall cyno_NormalRand1DArrayAntithetic(double* outArray, int size, short randflag);

// Plain-array APIs for cross-platform callers.
CYNOMATHUTILITY_API bool __stdcall cyno_Rand1DArrayPlain(double* outArray, long size);
CYNOMATHUTILITY_API bool __stdcall cyno_Sobol1DArrayPlain(double* outArray, long size);
CYNOMATHUTILITY_API bool __stdcall cyno_SobolNormal1DArrayPlain(double* outArray, long size);
CYNOMATHUTILITY_API bool __stdcall cyno_Sobol2DArrayPlain(double* outArray, long dim, long samplesPerFactor);
CYNOMATHUTILITY_API bool __stdcall cyno_SobolNormal2DArrayPlain(double* outArray, long dim, long samplesPerFactor);
CYNOMATHUTILITY_API bool __stdcall cyno_CorrelatedNormalRand1DArrayPlain(const double* correlationMatrix, long dim, double* outArray, short randflag);
CYNOMATHUTILITY_API bool __stdcall cyno_CorrelatedNormalRand2DArrayPlain(const double* correlationMatrix, long dim, double* outArray, long samplesPerFactor, short randflag);
CYNOMATHUTILITY_API bool __stdcall cyno_CorrelatedNormalRand1DArrayPlainAntithetic(const double* correlationMatrix, long dim, double* outArray, double* antitheticOutArray, short randflag);
CYNOMATHUTILITY_API bool __stdcall cyno_CorrelatedNormalRand2DArrayPlainAntithetic(const double* correlationMatrix, long dim, double* outArray, long samplesPerFactor, short randflag);
CYNOMATHUTILITY_API bool __stdcall cyno_EquityPathPlain(double spot, const double* dt, const double* drift, const double* dividend, const double* sigma, long nSteps, double* outPath, short randflag);
CYNOMATHUTILITY_API bool __stdcall cyno_EquityPathAntitheticPlain(double spot, const double* dt, const double* drift, const double* dividend, const double* sigma, long nSteps, double* outPath, double* antitheticOutPath, short randflag);
CYNOMATHUTILITY_API bool __stdcall cyno_EquityPathSobolPlain(double spot, const double* dt, const double* drift, const double* dividend, const double* sigma, long nSteps, double* outPath);
CYNOMATHUTILITY_API bool __stdcall cyno_EquityPathsPlain(const double* spots, long nStocks, const double* correlationMatrix, const double* dt, const double* drift, const double* dividend, const double* sigma, long nSteps, double* outPaths, short randflag);
CYNOMATHUTILITY_API bool __stdcall cyno_EquityPathsAntitheticPlain(const double* spots, long nStocks, const double* correlationMatrix, const double* dt, const double* drift, const double* dividend, const double* sigma, long nSteps, double* outPaths, double* antitheticOutPaths, short randflag);
CYNOMATHUTILITY_API bool __stdcall cyno_EquityPathsSobolPlain(const double* spots, long nStocks, const double* correlationMatrix, const double* dt, const double* drift, const double* dividend, const double* sigma, long nSteps, double* outPaths);

CYNOMATHUTILITY_API bool __stdcall cyno_MatrixScalePlain(double a, const double* input, long rows, long cols, double* output);
CYNOMATHUTILITY_API bool __stdcall cyno_MatrixMultiplyPlain(const double* a, long aRows, long aCols, const double* b, long bRows, long bCols, double* output);
CYNOMATHUTILITY_API bool __stdcall cyno_VectorDotPlain(const double* a, const double* b, long size, double* output);
CYNOMATHUTILITY_API bool __stdcall cyno_InverseMatrixPlain(const double* input, long dim, double* output);
CYNOMATHUTILITY_API bool __stdcall cyno_MatrixTransposePlain(const double* input, long rows, long cols, double* output);
CYNOMATHUTILITY_API bool __stdcall cyno_MatrixPowerPlain(double* inOutMatrix, long dim, short power);
CYNOMATHUTILITY_API bool __stdcall cyno_CholeskyDecompPlain(double* inOutMatrix, long dim);
CYNOMATHUTILITY_API bool __stdcall cyno_EigenSystemPlain(const double* input, long dim, double* vr, double* vi, double* eigenVectors, bool computeVects, double symTol, bool sortVects);
CYNOMATHUTILITY_API bool __stdcall cyno_PDCorrMatrixPlain(double* inOutMatrix, long dim);
CYNOMATHUTILITY_API bool __stdcall cyno_ApproximateCorrMatrixPlain(double* inOutMatrix, long dim);
CYNOMATHUTILITY_API bool __stdcall cyno_LinearSystemSolverPlain(const double* jacobian, const double* rhs, long dim, double* output);


//double StdNormCDF(double x);

CYNOMATHUTILITY_API double __stdcall cyno_Call(double S, double K, double totalSigma);
CYNOMATHUTILITY_API double __stdcall cyno_Put(double S, double K, double totalSigma);
CYNOMATHUTILITY_API double __stdcall cyno_DigitalCall(double S, double K, double totalSigma);
CYNOMATHUTILITY_API double __stdcall cyno_DigitalPut(double S, double K, double totalSigma);
CYNOMATHUTILITY_API double __stdcall cyno_BS(double Tau, double TD, double S, double K, double rD, double rF, double sigma, const char* callPut);
CYNOMATHUTILITY_API double __stdcall cyno_BS_Call(double Tau, double TD, double S, double K, double rD, double rF, double sigma);
CYNOMATHUTILITY_API double __stdcall cyno_BS_Put(double Tau, double TD, double S, double K, double rD, double rF, double sigma);

CYNOMATHUTILITY_API double __stdcall cyno_singleBarrier(double S, double K, double B, double Tau, double TD, double rd, double rf, double sigma, const char* knockType, const char* callput);
CYNOMATHUTILITY_API double __stdcall cyno_oneTouch(double Tau, double TD, double S, double B, double rd, double rf, double sigma, const char* upDown, const char* nCurrency, const char* nPaidAt);
CYNOMATHUTILITY_API double __stdcall cyno_hitProbabilityOneBarrier(double Tau, double TD, double S, double H, double rD, double rF, double sigma, const char* nCurrency);
CYNOMATHUTILITY_API double __stdcall cyno_oneNoTouch(double Tau, double TD, double S, double B, double rd, double rf, double sigma, const char* upDown, const char* nCurrency, const char* nPaidAt);

CYNOMATHUTILITY_API double __stdcall cyno_doubleBarrier(double Tau, double TD, double S, double K, double rD, double rF, double sigma,
    double L, double H, int knockTypeIn, int callPutIn);
CYNOMATHUTILITY_API double __stdcall cyno_doubleNoTouch(double Tau, double TD, double S, double rD, double rF, double sigma, double L, double H, int rCurrIn);
CYNOMATHUTILITY_API double __stdcall cyno_doubleOneTouch(double Tau, double TD, double S, double rD, double rF, double sigma, double L, double H, int rCurrIn);
CYNOMATHUTILITY_API double __stdcall cyno_hitProbabilityDoubleBarrier(double Tau, double TD, double S, double rD, double rF, double sigma, double L, double H, int nCurrencyIn);

#ifdef __cplusplus
}
#endif
