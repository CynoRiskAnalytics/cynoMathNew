#pragma once

CYNOMATHUTILITY_API bool __stdcall cyno_EquityPathPlain(double spot, const double* dt, const double* drift, const double* dividend, const double* sigma, long nSteps, double* outPath, short randflag);
CYNOMATHUTILITY_API bool __stdcall cyno_EquityPathAntitheticPlain(double spot, const double* dt, const double* drift, const double* dividend, const double* sigma, long nSteps, double* outPath, double* antitheticOutPath, short randflag);
CYNOMATHUTILITY_API bool __stdcall cyno_EquityPathSobolPlain(double spot, const double* dt, const double* drift, const double* dividend, const double* sigma, long nSteps, double* outPath);
CYNOMATHUTILITY_API bool __stdcall cyno_EquityPathsPlain(const double* spots, long nStocks, const double* correlationMatrix, const double* dt, const double* drift, const double* dividend, const double* sigma, long nSteps, double* outPaths, short randflag);
CYNOMATHUTILITY_API bool __stdcall cyno_EquityPathsAntitheticPlain(const double* spots, long nStocks, const double* correlationMatrix, const double* dt, const double* drift, const double* dividend, const double* sigma, long nSteps, double* outPaths, double* antitheticOutPaths, short randflag);
CYNOMATHUTILITY_API bool __stdcall cyno_EquityPathsSobolPlain(const double* spots, long nStocks, const double* correlationMatrix, const double* dt, const double* drift, const double* dividend, const double* sigma, long nSteps, double* outPaths);
