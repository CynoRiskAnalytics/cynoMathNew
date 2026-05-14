#include "pch.h"
#include "framework.h"
#include "equityPath.h"
#include "math.h"
#include <stdlib.h>

static double* equity_alloc_double_buffer(long size) {
    return (double*)malloc((size_t)size * sizeof(double));
}

static bool validate_equity_path_inputs(const double* dt, const double* drift, const double* dividend, const double* sigma, long nSteps, double* outPath) {
    long i;

    if (dt == NULL || drift == NULL || dividend == NULL || sigma == NULL || outPath == NULL || nSteps <= 0) {
        return false;
    }

    for (i = 0; i < nSteps; i++) {
        if (dt[i] < 0.0 || sigma[i] < 0.0) {
            return false;
        }
    }

    return true;
}

static void evolve_equity_path(double spot, const double* dt, const double* drift, const double* dividend, const double* sigma, const double* normals, long nSteps, double* outPath) {
    long i;

    outPath[0] = spot;
    for (i = 0; i < nSteps; i++) {
        double variance = sigma[i] * sigma[i];
        double step = (drift[i] - dividend[i] - 0.5 * variance) * dt[i]
            + sigma[i] * sqrt(dt[i]) * normals[i];
        outPath[i + 1] = outPath[i] * exp(step);
    }
}

CYNOMATHUTILITY_API bool __stdcall cyno_EquityPathPlain(double spot, const double* dt, const double* drift, const double* dividend, const double* sigma, long nSteps, double* outPath, short randflag) {
    double* normals;
    long i;

    if (!validate_equity_path_inputs(dt, drift, dividend, sigma, nSteps, outPath)) {
        return false;
    }

    normals = equity_alloc_double_buffer(nSteps);
    if (normals == NULL) {
        return false;
    }

    if (randflag == 1) {
        if (!cyno_NormalRand1DArray(normals, (int)nSteps)) {
            free(normals);
            return false;
        }
    }
    else {
        for (i = 0; i < nSteps; i++) {
            normals[i] = cyno_NormRandInv();
        }
    }

    evolve_equity_path(spot, dt, drift, dividend, sigma, normals, nSteps, outPath);

    free(normals);
    return true;
}

CYNOMATHUTILITY_API bool __stdcall cyno_EquityPathAntitheticPlain(double spot, const double* dt, const double* drift, const double* dividend, const double* sigma, long nSteps, double* outPath, double* antitheticOutPath, short randflag) {
    double* normals;
    double* antitheticNormals;
    long i;

    if (!validate_equity_path_inputs(dt, drift, dividend, sigma, nSteps, outPath) || antitheticOutPath == NULL) {
        return false;
    }

    normals = equity_alloc_double_buffer(nSteps);
    antitheticNormals = equity_alloc_double_buffer(nSteps);
    if (normals == NULL || antitheticNormals == NULL) {
        free(normals);
        free(antitheticNormals);
        return false;
    }

    if (randflag == 1) {
        if (!cyno_NormalRand1DArray(normals, (int)nSteps)) {
            free(normals);
            free(antitheticNormals);
            return false;
        }
    }
    else {
        for (i = 0; i < nSteps; i++) {
            normals[i] = cyno_NormRandInv();
        }
    }

    for (i = 0; i < nSteps; i++) {
        antitheticNormals[i] = -normals[i];
    }

    evolve_equity_path(spot, dt, drift, dividend, sigma, normals, nSteps, outPath);
    evolve_equity_path(spot, dt, drift, dividend, sigma, antitheticNormals, nSteps, antitheticOutPath);

    free(antitheticNormals);
    free(normals);
    return true;
}

CYNOMATHUTILITY_API bool __stdcall cyno_EquityPathSobolPlain(double spot, const double* dt, const double* drift, const double* dividend, const double* sigma, long nSteps, double* outPath) {
    double* normals;
    bool ok;

    if (!validate_equity_path_inputs(dt, drift, dividend, sigma, nSteps, outPath)) {
        return false;
    }

    normals = equity_alloc_double_buffer(nSteps);
    if (normals == NULL) {
        return false;
    }

    ok = cyno_SobolNormal1DArrayPlain(normals, nSteps);
    if (ok) {
        evolve_equity_path(spot, dt, drift, dividend, sigma, normals, nSteps, outPath);
    }

    free(normals);
    return ok;
}
