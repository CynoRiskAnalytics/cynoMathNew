#include "pch.h"
#include "framework.h"
#include <float.h>
#include <stdint.h>

//May 14, 2026 -Jesse Wang
//This is purely Codex generated, need to review it later


#define CYNO_SOBOL_MAX_BITS 32
#define CYNO_SOBOL_EXACT_DIM 8
#define CYNO_SOBOL_MAX_DIM 32

typedef struct SobolParams {
    int degree;
    unsigned int a;
    unsigned int m[5];
} SobolParams;

static const SobolParams SOBOL_PARAMS[CYNO_SOBOL_EXACT_DIM] = {
    {0, 0, {1, 0, 0, 0, 0}},
    {1, 0, {1, 0, 0, 0, 0}},
    {2, 1, {1, 3, 0, 0, 0}},
    {3, 1, {1, 3, 1, 0, 0}},
    {3, 2, {1, 1, 1, 0, 0}},
    {4, 1, {1, 1, 3, 3, 0}},
    {4, 4, {1, 3, 5, 13, 0}},
    {5, 2, {1, 1, 5, 5, 17}}
};

static const int HALTON_PRIMES[CYNO_SOBOL_MAX_DIM - CYNO_SOBOL_EXACT_DIM] = {
    23, 29, 31, 37, 41, 43, 47, 53,
    59, 61, 67, 71, 73, 79, 83, 89,
    97, 101, 103, 107, 109, 113, 127, 131
};

static int sobol_rightmost_zero_count(unsigned int value) {
    int c = 1;

    while (value & 1U) {
        value >>= 1U;
        c++;
    }

    return c;
}

static void sobol_build_direction_numbers(int dimIndex, uint32_t direction[CYNO_SOBOL_MAX_BITS + 1]) {
    int i;
    int j;
    int s = SOBOL_PARAMS[dimIndex].degree;
    unsigned int a = SOBOL_PARAMS[dimIndex].a;

    if (dimIndex == 0) {
        for (i = 1; i <= CYNO_SOBOL_MAX_BITS; i++) {
            direction[i] = 1U << (32 - i);
        }
        return;
    }

    for (i = 1; i <= s; i++) {
        direction[i] = SOBOL_PARAMS[dimIndex].m[i - 1] << (32 - i);
    }

    for (i = s + 1; i <= CYNO_SOBOL_MAX_BITS; i++) {
        uint32_t value = direction[i - s] ^ (direction[i - s] >> s);
        for (j = 1; j <= s - 1; j++) {
            if ((a >> (s - 1 - j)) & 1U) {
                value ^= direction[i - j];
            }
        }
        direction[i] = value;
    }
}

static double sobol_unit_from_uint32(uint32_t value) {
    return ((double)value) / 4294967296.0;
}

static double halton_value(long index, int base) {
    double result = 0.0;
    double factor = 1.0 / (double)base;

    while (index > 0) {
        result += factor * (double)(index % base);
        index /= base;
        factor /= (double)base;
    }

    return result;
}

static double sobol_clamp_open_unit(double value) {
    if (value <= 0.0) {
        return DBL_EPSILON;
    }
    if (value >= 1.0) {
        return 1.0 - DBL_EPSILON;
    }
    return value;
}

static bool sobol_fill_uniform(double* outArray, long dim, long samplesPerFactor) {
    uint32_t direction[CYNO_SOBOL_EXACT_DIM][CYNO_SOBOL_MAX_BITS + 1];
    uint32_t x[CYNO_SOBOL_EXACT_DIM];
    long sample;
    long d;

    if (outArray == NULL || dim <= 0 || dim > CYNO_SOBOL_MAX_DIM || samplesPerFactor <= 0) {
        return false;
    }

    for (d = 0; d < dim && d < CYNO_SOBOL_EXACT_DIM; d++) {
        sobol_build_direction_numbers((int)d, direction[d]);
        x[d] = 0U;
    }

    for (sample = 0; sample < samplesPerFactor; sample++) {
        if (sample == 0) {
            for (d = 0; d < dim; d++) {
                outArray[d * samplesPerFactor + sample] = 0.0;
            }
            continue;
        }

        {
            int c = sobol_rightmost_zero_count((unsigned int)(sample - 1));
            for (d = 0; d < dim; d++) {
                if (d < CYNO_SOBOL_EXACT_DIM) {
                    x[d] ^= direction[d][c];
                    outArray[d * samplesPerFactor + sample] = sobol_unit_from_uint32(x[d]);
                }
                else {
                    outArray[d * samplesPerFactor + sample] = halton_value(sample, HALTON_PRIMES[d - CYNO_SOBOL_EXACT_DIM]);
                }
            }
        }
    }

    return true;
}

CYNOMATHUTILITY_API bool __stdcall cyno_Sobol1DArrayPlain(double* outArray, long size) {
    return sobol_fill_uniform(outArray, 1, size);
}

CYNOMATHUTILITY_API bool __stdcall cyno_SobolNormal1DArrayPlain(double* outArray, long size) {
    long i;

    if (!cyno_Sobol1DArrayPlain(outArray, size)) {
        return false;
    }

    for (i = 0; i < size; i++) {
        outArray[i] = cyno_StdNormInv(sobol_clamp_open_unit(outArray[i]));
    }

    return true;
}

CYNOMATHUTILITY_API bool __stdcall cyno_Sobol2DArrayPlain(double* outArray, long dim, long samplesPerFactor) {
    return sobol_fill_uniform(outArray, dim, samplesPerFactor);
}

CYNOMATHUTILITY_API bool __stdcall cyno_SobolNormal2DArrayPlain(double* outArray, long dim, long samplesPerFactor) {
    long i;
    long total;

    if (!cyno_Sobol2DArrayPlain(outArray, dim, samplesPerFactor)) {
        return false;
    }

    total = dim * samplesPerFactor;
    for (i = 0; i < total; i++) {
        outArray[i] = cyno_StdNormInv(sobol_clamp_open_unit(outArray[i]));
    }

    return true;
}
