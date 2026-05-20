#pragma once

#include "framework.h"
#include <math.h>

static const double RADIX = 2.0;
static const double LUTINY = 1E-20;

void _stdcall cyno_MatrixMultiScalar1D(double a, long dim, double* b, double* c);
void _stdcall cyno_MatrixMultiScalar2D(double a, long dim1, long dim2, double** b, double** c);
bool _stdcall cyno_MatrixMulti2D1D(long adim1, long adim2, double** a, long bdim, double* b, double* c);
bool _stdcall cyno_MatrixMulti2D2D(long adim1, long adim2, double** a, long bdim1, long bdim2, double** b, double** c);
bool _stdcall cyno_MatrixMulti1D1D(long adim, double* a, long bdim, double* b, double* c);
void _stdcall cyno_MatrixTrans(long dim1, long dim2, double** a, double** b);

