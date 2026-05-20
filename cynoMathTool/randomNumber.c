
#include "pch.h"
#include "framework.h"
#include "math.h"
#include "matrix.h"
#include "MTRandomNumber.h"
#include "time.h"
#include <stdlib.h>


//Feb 5, 2020 - New random number generators by Jesse Wang
//May 19, 2026 -Jesse Wang
//Major updates and performance enhancement. Now this is much faster than both c# and java when random number is larger than 100,000

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

static void apply_lower_triangular_1d(const double* factorMatrix, long dim, const double* independentNormals, double* outArray) {
	long row;
	long col;

	for (row = 0; row < dim; row++) {
		double sum = 0.0;
		long rowOffset = row * dim;
		for (col = 0; col <= row; col++) {
			sum += factorMatrix[rowOffset + col] * independentNormals[col];
		}
		outArray[row] = sum;
	}
}

static void apply_lower_triangular_2d_sample(const double* factorMatrix, long dim, const double* independentNormals, double* outArray, long samplesPerFactor, long sampleIndex) {
	long row;
	long col;

	for (row = 0; row < dim; row++) {
		double sum = 0.0;
		long rowOffset = row * dim;
		for (col = 0; col <= row; col++) {
			sum += factorMatrix[rowOffset + col] * independentNormals[col];
		}
		outArray[row * samplesPerFactor + sampleIndex] = sum;
	}
}

static void negate_correlated_sample(double* outArray, long dim, long samplesPerFactor, long sourceIndex, long targetIndex) {
	long row;

	for (row = 0; row < dim; row++) {
		outArray[row * samplesPerFactor + targetIndex] = -outArray[row * samplesPerFactor + sourceIndex];
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

	if (correlationMatrix == NULL || outArray == NULL || dim <= 0) {
		return false;
	}

	randomVector = alloc_double_buffer(dim);
	if (randomVector == NULL) {
		free(randomVector);
		return false;
	}

	fill_random_vector(randomVector, dim, randflag);
	apply_lower_triangular_1d(correlationMatrix, dim, randomVector, outArray);
	free(randomVector);
	return true;
}

CYNOMATHUTILITY_API bool __stdcall cyno_CorrelatedNormalRand1DArrayPlainAntithetic(const double* correlationMatrix, long dim, double* outArray, double* antitheticOutArray, short randflag) {
	double* randomVector;

	if (correlationMatrix == NULL || outArray == NULL || antitheticOutArray == NULL || dim <= 0) {
		return false;
	}

	randomVector = alloc_double_buffer(dim);
	if (randomVector == NULL) {
		free(randomVector);
		return false;
	}

	fill_random_vector(randomVector, dim, randflag);
	apply_lower_triangular_1d(correlationMatrix, dim, randomVector, outArray);
	for (long i = 0; i < dim; i++) {
		antitheticOutArray[i] = -outArray[i];
	}

	free(randomVector);
	return true;
}

CYNOMATHUTILITY_API bool __stdcall cyno_CorrelatedNormalRand2DArrayPlain(const double* correlationMatrix, long dim, double* outArray, long samplesPerFactor, short randflag) {
	double* randomVector;
	long sample;
	long row;
	long col;

	if (correlationMatrix == NULL || outArray == NULL || dim <= 0 || samplesPerFactor <= 0) {
		return false;
	}

	randomVector = alloc_double_buffer(dim);
	if (randomVector == NULL) {
		return false;
	}

	for (sample = 0; sample < samplesPerFactor; sample++) {
		for (col = 0; col < dim; col++) {
			randomVector[col] = (randflag == 1) ? cyno_NormRand() : cyno_NormRandInv();
		}
		for (row = 0; row < dim; row++) {
			double sum = 0.0;
			long rowOffset = row * dim;
			for (col = 0; col <= row; col++) {
				sum += correlationMatrix[rowOffset + col] * randomVector[col];
			}
			outArray[row * samplesPerFactor + sample] = sum;
		}
	}

	free(randomVector);
	return true;
}

CYNOMATHUTILITY_API bool __stdcall cyno_CorrelatedNormalRand2DArrayPlainAntithetic(const double* correlationMatrix, long dim, double* outArray, long samplesPerFactor, short randflag) {
	double* randomVector;
	long sample;

	if (correlationMatrix == NULL || outArray == NULL || dim <= 0 || samplesPerFactor <= 0) {
		return false;
	}

	randomVector = alloc_double_buffer(dim);
	if (randomVector == NULL) {
		return false;
	}

	for (sample = 0; sample < samplesPerFactor; sample += 2) {
		fill_random_vector(randomVector, dim, randflag);
		apply_lower_triangular_2d_sample(correlationMatrix, dim, randomVector, outArray, samplesPerFactor, sample);
		if (sample + 1 < samplesPerFactor) {
			negate_correlated_sample(outArray, dim, samplesPerFactor, sample, sample + 1);
		}
	}

	free(randomVector);
	return true;
}

