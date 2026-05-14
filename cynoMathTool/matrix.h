#include "pch.h"
#include "framework.h"
#include "math.h"


bool cyno_ConvertFromSafeArray1D(SAFEARRAY** a, long dim, double* pa);

bool cyno_ConvertFromSafeArray2D(SAFEARRAY** a, long dim1, long dim2, double** pa);
#include "pch.h"
#include "framework.h"
#include "math.h"


static const double RADIX = 2;
static const double LUTINY = 1E-20;

bool cyno_ConvertToSafeArray1D(SAFEARRAY** a, long dim, double* pa);

bool cyno_ConvertToSafeArray2D(SAFEARRAY** a, long dim1, long dim2, double** pa);

void _stdcall cyno_MatrixMultiScalar1D(double a, long dim, double* b, double* c);

void _stdcall cyno_MatrixMultiScalar2D(double a, long dim1, long dim2, double** b, double** c);

bool _stdcall cyno_MatrixMulti2D1D(long adim1, long adim2, double** a, long bdim, double* b, double* c);

bool _stdcall cyno_MatrixMulti2D2D(long adim1, long adim2, double** a, long bdim1, long bdim2, double** b, double** c);

bool _stdcall cyno_MatrixMulti1D1D(long adim, double* a, long bdim, double* b, double* c);

void _stdcall cyno_MatrixTrans(long dim1, long dim2, double** a, double** b);

