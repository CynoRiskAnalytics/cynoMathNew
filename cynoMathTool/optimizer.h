/* optimizer.h
Feb 9, 2020 - created by Jesse Wang

*/


#include "pch.h"
#include "framework.h"
#include "math.h"


static const double LSALF = 1e-4;
static const double LSTOLX = 1e-7;

CYNOMATHUTILITY_API bool __stdcall cyno_lineSearch(long func, int numOfDim, double* xold, double* fold, double* g, double* p, double* x, double* f, double max_step);
