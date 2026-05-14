
#include "pch.h"
#include "framework.h"
#include "math.h"
#include "matrix.h"
#include "MTRandomNumber.h"
#include "time.h"
#include <stdlib.h>


//Feb 5, 2020 - New random number generators by Jesse Wang

typedef unsigned long long Ullong;
//typedef unsigned __int64 Ullong;


typedef struct Rand2 {
	Ullong v;
	Ullong w;
} Rand2;

static Rand2* p = NULL;
static short g_normRandHasSpare = 0;
static double g_normRandSpare = 0.0;

static Ullong rand2_int64(Rand2* state) {
	state->v ^= state->v >> 17;
	state->v ^= state->v << 31;
	state->v ^= state->v >> 8;
	state->w = 4294957665U * (state->w & 0xffffffff) + (state->w >> 32);
	return state->v ^ state->w;
}

static void rand2_seed(Rand2* state, Ullong seed) {
	state->v = 4101842887655102017ULL;
	state->w = 1;
	state->v ^= seed;
	state->w = rand2_int64(state);
	state->v = rand2_int64(state);
}

static double rand2_doub(Rand2* state) {
	return 5.42101086242752217E-20 * rand2_int64(state);
}

static double* alloc_double_buffer(long size) {
	return (double*)malloc((size_t)size * sizeof(double));
}

static double** alloc_row_view(double* data, long rows, long cols) {
	double** row_view;
	long i;

	row_view = (double**)malloc((size_t)rows * sizeof(double*));
	if (row_view == NULL) {
		return NULL;
	}

	for (i = 0; i < rows; i++) {
		row_view[i] = data + (i * cols);
	}

	return row_view;
}

static void fill_random_vector(double* outArray, long size, short randflag) {
	long i;

	for (i = 0; i < size; i++) {
		outArray[i] = (randflag == 1) ? cyno_NormRand() : cyno_NormRandInv();
	}
}

static void fill_random_vector_antithetic(double* outArray, long size, short randflag) {
	long i;

	for (i = 0; i < size; i += 2) {
		double sample = (randflag == 1) ? cyno_NormRand() : cyno_NormRandInv();
		outArray[i] = sample;
		if (i + 1 < size) {
			outArray[i + 1] = -sample;
		}
	}
}

static void fill_uniform_vector_antithetic(double* outArray, long size) {
	long i;

	for (i = 0; i < size; i += 2) {
		double sample = cyno_Rand();
		outArray[i] = sample;
		if (i + 1 < size) {
			outArray[i + 1] = 1.0 - sample;
		}
	}
}



//Purpose
//		Set random number generator cyno_Rand() seed.
//Parameters
//		seed	input	seed value
//Return Values
//		none
//Remarks
// Numeric Receipe 3 ranq2
CYNOMATHUTILITY_API void __stdcall cyno_Rand_Seed(long seed) {

	free(p);
	p = (Rand2*)malloc(sizeof(Rand2));
	g_normRandHasSpare = 0;
	g_normRandSpare = 0.0;
	if (p != NULL) {
		rand2_seed(p, (Ullong)seed);
	}
}

//Purpose
//		Generate random number between 0 and 1
//Parameters
//		none
//Return Values
//		double
//Remarks
//		ranq2 from Numerica Receipe 3.
//		if seed is not initialized explicitly, it will be initialized using system time.
//Numeric Receipe 3 ranq2

CYNOMATHUTILITY_API double __stdcall cyno_Rand(void) {

	if (p == NULL) {
		p = (Rand2*)malloc(sizeof(Rand2));
		if (p == NULL) {
			return 0.0;
		}
		rand2_seed(p, (Ullong)time(0));
	}
	return rand2_doub(p);
}


//Purpose
//		Generate random number between 0 and 1 by MT.
//Parameters
//		none
//Return Values
//		double
//Remarks
//		mt19937
// Feb 5, 2020 - newly generated using the latest 64 bit MT 

CYNOMATHUTILITY_API void __stdcall cyno_RandMT_Seed(long seed) {

	init_genrand64(seed);

}


CYNOMATHUTILITY_API double __stdcall cyno_MTRand(void) {
	return genrand64_real1();
}

//Purpose
//		Generates random normally distributed deviate (with zero mean, and unit variance) using
//		Box-Muller transform
//Parameters
//		none
//Return Values
//		double


CYNOMATHUTILITY_API double __stdcall cyno_NormRand(void) {
	double fac, rsq, v1, v2;


	if (g_normRandHasSpare == 0) {

		do {
			v1 = (double)2 * cyno_Rand() - (double)1;
			v2 = (double)2 * cyno_Rand() - (double)1;
			rsq = v1 * v1 + v2 * v2;
		} while (rsq >= (double)1 || rsq == (double)0);
		fac = sqrt((double)-2 * log(rsq) / rsq);
		g_normRandSpare = v1 * fac;
		g_normRandHasSpare = 1;

		return v2 * fac;
	}
	else {
		g_normRandHasSpare = 0;
		return g_normRandSpare;
	}
}

//Purpose
//		Generates normal random number by inverse function
//Parameters
//		none
//Return Values
//		double
//Remarks

CYNOMATHUTILITY_API double __stdcall cyno_NormRandInv(void) {
	return cyno_StdNormInv(cyno_Rand());
}


//Purpose
//		Generates 1D correlated normal random number by either BM or normal inverse function
//Parameters
//		psaMatrix	input	correlation matrix
//		psaOutArr	output	1D array that stores the generated random numbers
//		randflag	input	Which method is used to generate random number. 1(default)=Box-Muller other value=normal inverse function
//Return Values
//		bool
//Remarks
//		1) The output array psaOutArr must be one dimension and be sized properly
//			before call this function.
//		2) C++ bool corresponds to Byte in VBA
//			true - 1 or any non zero
//			false - 0
//			When VBA references this function, return type should be Byte
//			should the need arise, A wrapper function should be created to convert Byte to boolean

CYNOMATHUTILITY_API bool _stdcall cyno_CorrelatedNormalRand1DArray(SAFEARRAY** psaMatrix, SAFEARRAY** psaOutArr, short randflag) {
	long i;
	long size, dim1, dim2;
	bool ret;
	double* randomVector;
	double* correlatedRandomVector;
	double* pMatrixData;
	double** pMatrix;

	if ((*psaOutArr)->cDims == 1 && (*psaMatrix)->cDims == 2) {
		size = (*psaOutArr)->rgsabound[0].cElements;
		randomVector = alloc_double_buffer(size);
		correlatedRandomVector = alloc_double_buffer(size);
		if (randomVector == NULL || correlatedRandomVector == NULL) {
			free(randomVector);
			free(correlatedRandomVector);
			return false;
		}
		if (randflag == 1) {
			for (i = 0; i < size; i++) randomVector[i] = cyno_NormRand();
		}
		else {
			for (i = 0; i < size; i++) randomVector[i] = cyno_NormRandInv();
		}

		dim1 = (*psaMatrix)->rgsabound[1].cElements;
		dim2 = (*psaMatrix)->rgsabound[0].cElements;

		pMatrixData = alloc_double_buffer(dim1 * dim2);
		pMatrix = (pMatrixData == NULL) ? NULL : alloc_row_view(pMatrixData, dim1, dim2);
		if (pMatrix == NULL) {
			free(randomVector);
			free(correlatedRandomVector);
			free(pMatrixData);
			return false;
		}

		ret = cyno_ConvertFromSafeArray2D(psaMatrix, dim1, dim2, pMatrix);
		ret = cyno_MatrixMulti2D1D(dim1, dim2, pMatrix, size, randomVector, correlatedRandomVector);
		ret = cyno_ConvertToSafeArray1D(psaOutArr, size, correlatedRandomVector);

		free(pMatrix);
		free(pMatrixData);
		free(correlatedRandomVector);
		free(randomVector);

		return ret;
	}

	return false;
}

//Purpose
//		Generates 1D random number
//Parameters
//		psaOutArr	output	1D array that stores the generated random numbers
//Return Values
//		bool
//Remarks
//		The output array psaOutArr must be one dimension and be sized properly
//		before call this function.

CYNOMATHUTILITY_API bool _stdcall cyno_Rand1DArray(SAFEARRAY** psaOutArr) {
	long i;
	long size;
	bool ret;
	double* randomVector;

	if ((*psaOutArr)->cDims == 1) {
		size = (*psaOutArr)->rgsabound[0].cElements;
		randomVector = alloc_double_buffer(size);
		if (randomVector == NULL) {
			return false;
		}
		for (i = 0; i < size; i++) {
			randomVector[i] = cyno_Rand();
		}

		ret = cyno_ConvertToSafeArray1D(psaOutArr, size, randomVector);
		free(randomVector);

		return ret;
	}
	return false;
}

//Purpose
//		Generates 1D normal random number by BM -can be used to simulate one path
//Parameters
//		psaOutArr	output	1D array that stores the generated normal random numbers
//Return Values
//		bool
//Remarks
//		The output array psaOutArr must be one dimension and be sized properly
//		before call this function.
//CYNOMATHUTILITY_API bool _stdcall cyno_NormalRand1DArray(SAFEARRAY** psaOutArr) {
//	long i;
//	long size;
//	double* randomVector;
//	bool ret;
//
//	if ((*psaOutArr)->cDims == 1) {
//		size = (*psaOutArr)->rgsabound[0].cElements;
//		randomVector = new double[size];
//		for (i = 0; i < size; i++) {
//			randomVector[i] = cyno_NormRand();
//		}
//
//		ret = cyno_ConvertToSafeArray1D(psaOutArr, size, randomVector);
//
//		delete[] randomVector;
//
//		return ret;
//	}
//	return false;
//}

//CYNOMATHUTILITY_API bool _stdcall cyno_NormalRand1DArray(SAFEARRAY* psaOutArr) {
//	if (psaOutArr == nullptr || psaOutArr->cDims != 1) return false;
//
//	double* pData = nullptr;
//	SafeArrayAccessData(psaOutArr, (void**)&pData);
//
//	if (!pData) return false;
//
//	long size = psaOutArr->rgsabound[0].cElements;
//	for (long i = 0; i < size; ++i)
//		pData[i] = cyno_NormRand(); // your custom random generator
//
//	SafeArrayUnaccessData(psaOutArr);
//	return true;
//}

CYNOMATHUTILITY_API bool _stdcall cyno_NormalRand1DArray(double* outArray, int size) {
	if (!outArray || size <= 0 )
		return false;

	for (int i = 0; i < size; i++) outArray[i] = cyno_NormRand();
		
	return true;
}

CYNOMATHUTILITY_API bool __stdcall cyno_NormalRand1DArrayAntithetic(double* outArray, int size, short randflag) {
	if (!outArray || size <= 0) {
		return false;
	}

	fill_random_vector_antithetic(outArray, size, randflag);
	return true;
}


//Purpose
//		Generates 2D correlated normal random number by either BM or normal inverse function
//		This function has a special usage for Razor simulation
//Parameters
//		psaMatrix	in	correlation matrix
//		psaOutArr	out	2D array that stores the generated normal random numbers
//		randflag	in	Which method is used to generate random number. 1(default)=Box-Muller other value=normal inverse function
//Return Values
//		bool
//Remarks
//		The output array psaOutArr must be two dimension and be sized properly
//		before call this function. First dimension of psaOutArr is number of risk factors, second dimension is number of
//		time nodes

CYNOMATHUTILITY_API bool _stdcall cyno_CorrelatedNormalRand2DArray(SAFEARRAY** psaMatrix, SAFEARRAY** psaOutArr, short randflag) {
	long i, j;
	long size1, size2, dim1, dim2;
	bool ret;
	double* randomVectorData;
	double* correlatedRandomVectorData;
	double** randomVector;
	double** correlatedRandomVector;
	double* pMatrixData;
	double** pMatrix;

	if ((*psaOutArr)->cDims == 2 && (*psaMatrix)->cDims == 2) {
		size1 = (*psaOutArr)->rgsabound[1].cElements;
		size2 = (*psaOutArr)->rgsabound[0].cElements;

		dim1 = (*psaMatrix)->rgsabound[1].cElements;
		dim2 = (*psaMatrix)->rgsabound[0].cElements;

		//Correlation Matrix must be a square.
		if (dim1 != dim2)
			return false;

		randomVectorData = alloc_double_buffer(size1 * size2);
		correlatedRandomVectorData = alloc_double_buffer(size1 * size2);
		randomVector = (randomVectorData == NULL) ? NULL : alloc_row_view(randomVectorData, size1, size2);
		correlatedRandomVector = (correlatedRandomVectorData == NULL) ? NULL : alloc_row_view(correlatedRandomVectorData, size1, size2);
		if (randomVectorData == NULL || correlatedRandomVectorData == NULL || randomVector == NULL || correlatedRandomVector == NULL) {
			free(randomVector);
			free(correlatedRandomVector);
			free(randomVectorData);
			free(correlatedRandomVectorData);
			return false;
		}

		if (randflag == 1) {
			for (i = 0; i < size1; i++)
				for (j = 0; j < size2; j++)
					randomVector[i][j] = cyno_NormRand();
		}
		else {
			for (i = 0; i < size1; i++)
				for (j = 0; j < size2; j++)
					randomVector[i][j] = cyno_NormRandInv();
		}


		pMatrixData = alloc_double_buffer(dim1 * dim2);
		pMatrix = (pMatrixData == NULL) ? NULL : alloc_row_view(pMatrixData, dim1, dim2);
		if (pMatrix == NULL) {
			free(randomVector);
			free(correlatedRandomVector);
			free(randomVectorData);
			free(correlatedRandomVectorData);
			free(pMatrixData);
			return false;
		}

		ret = cyno_ConvertFromSafeArray2D(psaMatrix, dim1, dim2, pMatrix);

		if (dim1 == size1)
			ret = cyno_MatrixMulti2D2D(dim1, dim2, pMatrix, size1, size2, randomVector, correlatedRandomVector);
		else
			ret = cyno_MatrixMulti2D2D(size1, size2, randomVector, dim1, dim2, pMatrix, correlatedRandomVector);

		ret = cyno_ConvertToSafeArray2D(psaOutArr, size1, size2, correlatedRandomVector);

		free(pMatrix);
		free(pMatrixData);
		free(randomVector);
		free(correlatedRandomVector);
		free(randomVectorData);
		free(correlatedRandomVectorData);

		return ret;
	}
	return false;
}

CYNOMATHUTILITY_API bool __stdcall cyno_Rand1DArrayPlain(double* outArray, long size) {
	long i;

	if (outArray == NULL || size <= 0) {
		return false;
	}

	for (i = 0; i < size; i++) {
		outArray[i] = cyno_Rand();
	}

	return true;
}

CYNOMATHUTILITY_API bool __stdcall cyno_Rand1DArrayPlainAntithetic(double* outArray, long size) {
	if (outArray == NULL || size <= 0) {
		return false;
	}

	fill_uniform_vector_antithetic(outArray, size);
	return true;
}


CYNOMATHUTILITY_API bool __stdcall cyno_CorrelatedNormalRand1DArrayPlain(const double* correlationMatrix, long dim, double* outArray, short randflag) {
	double* randomVector;
	double* matrixData;
	double** matrixRows;
	bool ret;

	if (correlationMatrix == NULL || outArray == NULL || dim <= 0) {
		return false;
	}

	randomVector = alloc_double_buffer(dim);
	matrixData = alloc_double_buffer(dim * dim);
	matrixRows = (matrixData == NULL) ? NULL : alloc_row_view(matrixData, dim, dim);
	if (randomVector == NULL || matrixData == NULL || matrixRows == NULL) {
		free(randomVector);
		free(matrixData);
		free(matrixRows);
		return false;
	}

	memcpy(matrixData, correlationMatrix, (size_t)(dim * dim) * sizeof(double));
	fill_random_vector(randomVector, dim, randflag);
	ret = cyno_MatrixMulti2D1D(dim, dim, matrixRows, dim, randomVector, outArray);

	free(matrixRows);
	free(matrixData);
	free(randomVector);
	return ret;
}

CYNOMATHUTILITY_API bool __stdcall cyno_CorrelatedNormalRand1DArrayPlainAntithetic(const double* correlationMatrix, long dim, double* outArray, double* antitheticOutArray, short randflag) {
	double* randomVector;
	double* antitheticVector;
	double* matrixData;
	double** matrixRows;
	bool ret;

	if (correlationMatrix == NULL || outArray == NULL || antitheticOutArray == NULL || dim <= 0) {
		return false;
	}

	randomVector = alloc_double_buffer(dim);
	antitheticVector = alloc_double_buffer(dim);
	matrixData = alloc_double_buffer(dim * dim);
	matrixRows = (matrixData == NULL) ? NULL : alloc_row_view(matrixData, dim, dim);
	if (randomVector == NULL || antitheticVector == NULL || matrixData == NULL || matrixRows == NULL) {
		free(randomVector);
		free(antitheticVector);
		free(matrixData);
		free(matrixRows);
		return false;
	}

	memcpy(matrixData, correlationMatrix, (size_t)(dim * dim) * sizeof(double));
	fill_random_vector(randomVector, dim, randflag);
	for (long i = 0; i < dim; i++) {
		antitheticVector[i] = -randomVector[i];
	}

	ret = cyno_MatrixMulti2D1D(dim, dim, matrixRows, dim, randomVector, outArray);
	if (ret) {
		ret = cyno_MatrixMulti2D1D(dim, dim, matrixRows, dim, antitheticVector, antitheticOutArray);
	}

	free(matrixRows);
	free(matrixData);
	free(antitheticVector);
	free(randomVector);
	return ret;
}

CYNOMATHUTILITY_API bool __stdcall cyno_CorrelatedNormalRand2DArrayPlain(const double* correlationMatrix, long dim, double* outArray, long samplesPerFactor, short randflag) {
	double* matrixData;
	double* randomData;
	double** matrixRows;
	double** randomRows;
	double** outputRows;
	bool ret;
	long i;

	if (correlationMatrix == NULL || outArray == NULL || dim <= 0 || samplesPerFactor <= 0) {
		return false;
	}

	matrixData = alloc_double_buffer(dim * dim);
	randomData = alloc_double_buffer(dim * samplesPerFactor);
	matrixRows = (matrixData == NULL) ? NULL : alloc_row_view(matrixData, dim, dim);
	randomRows = (randomData == NULL) ? NULL : alloc_row_view(randomData, dim, samplesPerFactor);
	outputRows = alloc_row_view(outArray, dim, samplesPerFactor);
	if (matrixData == NULL || randomData == NULL || matrixRows == NULL || randomRows == NULL || outputRows == NULL) {
		free(matrixData);
		free(randomData);
		free(matrixRows);
		free(randomRows);
		free(outputRows);
		return false;
	}

	memcpy(matrixData, correlationMatrix, (size_t)(dim * dim) * sizeof(double));
	for (i = 0; i < dim * samplesPerFactor; i++) {
		randomData[i] = (randflag == 1) ? cyno_NormRand() : cyno_NormRandInv();
	}

	ret = cyno_MatrixMulti2D2D(dim, dim, matrixRows, dim, samplesPerFactor, randomRows, outputRows);

	free(outputRows);
	free(randomRows);
	free(matrixRows);
	free(randomData);
	free(matrixData);
	return ret;
}

CYNOMATHUTILITY_API bool __stdcall cyno_CorrelatedNormalRand2DArrayPlainAntithetic(const double* correlationMatrix, long dim, double* outArray, long samplesPerFactor, short randflag) {
	double* matrixData;
	double* randomData;
	double** matrixRows;
	double** randomRows;
	double** outputRows;
	bool ret;
	long factor;
	long sample;

	if (correlationMatrix == NULL || outArray == NULL || dim <= 0 || samplesPerFactor <= 0) {
		return false;
	}

	matrixData = alloc_double_buffer(dim * dim);
	randomData = alloc_double_buffer(dim * samplesPerFactor);
	matrixRows = (matrixData == NULL) ? NULL : alloc_row_view(matrixData, dim, dim);
	randomRows = (randomData == NULL) ? NULL : alloc_row_view(randomData, dim, samplesPerFactor);
	outputRows = alloc_row_view(outArray, dim, samplesPerFactor);
	if (matrixData == NULL || randomData == NULL || matrixRows == NULL || randomRows == NULL || outputRows == NULL) {
		free(matrixData);
		free(randomData);
		free(matrixRows);
		free(randomRows);
		free(outputRows);
		return false;
	}

	memcpy(matrixData, correlationMatrix, (size_t)(dim * dim) * sizeof(double));
	for (factor = 0; factor < dim; factor++) {
		fill_random_vector_antithetic(randomRows[factor], samplesPerFactor, randflag);
	}

	ret = cyno_MatrixMulti2D2D(dim, dim, matrixRows, dim, samplesPerFactor, randomRows, outputRows);

	free(outputRows);
	free(randomRows);
	free(matrixRows);
	free(randomData);
	free(matrixData);
	return ret;
}

