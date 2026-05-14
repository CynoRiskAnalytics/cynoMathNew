#include "pch.h"
#include "framework.h"
#include "equityPath.h"
#include "math.h"
#include <stdlib.h>

static double* equity_alloc_double_buffer(long size) {
    return (double*)malloc((size_t)size * sizeof(double));
}

static double path_input_at(const double* values, long stockIndex, long stepIndex, long nSteps) {
    return values[stockIndex * nSteps + stepIndex];
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

static bool validate_multi_equity_path_inputs(const double* spots, long nStocks, const double* correlationMatrix, const double* dt, const double* drift, const double* dividend, const double* sigma, long nSteps, double* outPaths) {
    long i;

    if (spots == NULL || correlationMatrix == NULL || outPaths == NULL || nStocks <= 0) {
        return false;
    }

    if (!validate_equity_path_inputs(dt, drift, dividend, sigma, nSteps, outPaths)) {
        return false;
    }

    for (i = 0; i < nStocks; i++) {
        if (spots[i] <= 0.0) {
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

static void evolve_multi_equity_paths(const double* spots, long nStocks, const double* dt, const double* drift, const double* dividend, const double* sigma, const double* correlatedNormals, long nSteps, double* outPaths) {
    long stock;
    long step;

    for (stock = 0; stock < nStocks; stock++) {
        outPaths[stock * (nSteps + 1)] = spots[stock];
    }

    for (step = 0; step < nSteps; step++) {
        for (stock = 0; stock < nStocks; stock++) {
            long pathOffset = stock * (nSteps + 1);
            double sigmaStep = path_input_at(sigma, stock, step, nSteps);
            double variance = sigmaStep * sigmaStep;
            double z = correlatedNormals[step * nStocks + stock];
            double driftStep = path_input_at(drift, stock, step, nSteps);
            double dividendStep = path_input_at(dividend, stock, step, nSteps);
            double stepTerm = (driftStep - dividendStep - 0.5 * variance) * dt[step]
                + sigmaStep * sqrt(dt[step]) * z;
            outPaths[pathOffset + step + 1] = outPaths[pathOffset + step] * exp(stepTerm);
        }
    }
}

static void apply_lower_triangular_correlation(const double* choleskyMatrix, long nStocks, const double* independentNormals, double* correlatedNormals) {
    long row;
    long col;

    for (row = 0; row < nStocks; row++) {
        double sum = 0.0;
        for (col = 0; col <= row; col++) {
            sum += choleskyMatrix[row * nStocks + col] * independentNormals[col];
        }
        correlatedNormals[row] = sum;
    }
}

static bool build_cholesky_matrix(const double* correlationMatrix, long nStocks, double** choleskyOut) {
    double* matrixCopy;

    matrixCopy = equity_alloc_double_buffer(nStocks * nStocks);
    if (matrixCopy == NULL) {
        return false;
    }

    memcpy(matrixCopy, correlationMatrix, (size_t)(nStocks * nStocks) * sizeof(double));
    if (!cyno_CholeskyDecompPlain(matrixCopy, nStocks)) {
        free(matrixCopy);
        return false;
    }

    *choleskyOut = matrixCopy;
    return true;
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

CYNOMATHUTILITY_API bool __stdcall cyno_EquityPathsPlain(const double* spots, long nStocks, const double* correlationMatrix, const double* dt, const double* drift, const double* dividend, const double* sigma, long nSteps, double* outPaths, short randflag) {
    double* choleskyMatrix;
    double* independentNormals;
    double* correlatedNormals;
    long step;
    long stock;
    bool ok = false;

    if (!validate_multi_equity_path_inputs(spots, nStocks, correlationMatrix, dt, drift, dividend, sigma, nSteps, outPaths)) {
        return false;
    }

    if (!build_cholesky_matrix(correlationMatrix, nStocks, &choleskyMatrix)) {
        return false;
    }

    independentNormals = equity_alloc_double_buffer(nStocks);
    correlatedNormals = equity_alloc_double_buffer(nSteps * nStocks);
    if (independentNormals == NULL || correlatedNormals == NULL) {
        free(independentNormals);
        free(correlatedNormals);
        free(choleskyMatrix);
        return false;
    }

    for (step = 0; step < nSteps; step++) {
        if (randflag == 1) {
            if (!cyno_NormalRand1DArray(independentNormals, (int)nStocks)) {
                goto cleanup;
            }
        }
        else {
            for (stock = 0; stock < nStocks; stock++) {
                independentNormals[stock] = cyno_NormRandInv();
            }
        }

        apply_lower_triangular_correlation(choleskyMatrix, nStocks, independentNormals, correlatedNormals + step * nStocks);
    }

    evolve_multi_equity_paths(spots, nStocks, dt, drift, dividend, sigma, correlatedNormals, nSteps, outPaths);
    ok = true;

cleanup:
    free(correlatedNormals);
    free(independentNormals);
    free(choleskyMatrix);
    return ok;
}

CYNOMATHUTILITY_API bool __stdcall cyno_EquityPathsAntitheticPlain(const double* spots, long nStocks, const double* correlationMatrix, const double* dt, const double* drift, const double* dividend, const double* sigma, long nSteps, double* outPaths, double* antitheticOutPaths, short randflag) {
    double* choleskyMatrix;
    double* independentNormals;
    double* correlatedNormals;
    double* correlatedAntitheticNormals;
    long step;
    long stock;
    bool ok = false;

    if (!validate_multi_equity_path_inputs(spots, nStocks, correlationMatrix, dt, drift, dividend, sigma, nSteps, outPaths) || antitheticOutPaths == NULL) {
        return false;
    }

    if (!build_cholesky_matrix(correlationMatrix, nStocks, &choleskyMatrix)) {
        return false;
    }

    independentNormals = equity_alloc_double_buffer(nStocks);
    correlatedNormals = equity_alloc_double_buffer(nSteps * nStocks);
    correlatedAntitheticNormals = equity_alloc_double_buffer(nSteps * nStocks);
    if (independentNormals == NULL || correlatedNormals == NULL || correlatedAntitheticNormals == NULL) {
        free(independentNormals);
        free(correlatedNormals);
        free(correlatedAntitheticNormals);
        free(choleskyMatrix);
        return false;
    }

    for (step = 0; step < nSteps; step++) {
        if (randflag == 1) {
            if (!cyno_NormalRand1DArray(independentNormals, (int)nStocks)) {
                goto cleanup_multi_antithetic;
            }
        }
        else {
            for (stock = 0; stock < nStocks; stock++) {
                independentNormals[stock] = cyno_NormRandInv();
            }
        }

        apply_lower_triangular_correlation(choleskyMatrix, nStocks, independentNormals, correlatedNormals + step * nStocks);
        for (stock = 0; stock < nStocks; stock++) {
            correlatedAntitheticNormals[step * nStocks + stock] = -correlatedNormals[step * nStocks + stock];
        }
    }

    evolve_multi_equity_paths(spots, nStocks, dt, drift, dividend, sigma, correlatedNormals, nSteps, outPaths);
    evolve_multi_equity_paths(spots, nStocks, dt, drift, dividend, sigma, correlatedAntitheticNormals, nSteps, antitheticOutPaths);
    ok = true;

cleanup_multi_antithetic:
    free(correlatedAntitheticNormals);
    free(correlatedNormals);
    free(independentNormals);
    free(choleskyMatrix);
    return ok;
}

CYNOMATHUTILITY_API bool __stdcall cyno_EquityPathsSobolPlain(const double* spots, long nStocks, const double* correlationMatrix, const double* dt, const double* drift, const double* dividend, const double* sigma, long nSteps, double* outPaths) {
    double* choleskyMatrix;
    double* sobolNormals;
    double* stepNormals;
    double* correlatedNormals;
    long step;
    long stock;
    bool ok = false;

    if (!validate_multi_equity_path_inputs(spots, nStocks, correlationMatrix, dt, drift, dividend, sigma, nSteps, outPaths)) {
        return false;
    }

    if (!build_cholesky_matrix(correlationMatrix, nStocks, &choleskyMatrix)) {
        return false;
    }

    sobolNormals = equity_alloc_double_buffer(nStocks * nSteps);
    stepNormals = equity_alloc_double_buffer(nStocks);
    correlatedNormals = equity_alloc_double_buffer(nStocks * nSteps);
    if (sobolNormals == NULL || stepNormals == NULL || correlatedNormals == NULL) {
        free(sobolNormals);
        free(stepNormals);
        free(correlatedNormals);
        free(choleskyMatrix);
        return false;
    }

    if (!cyno_SobolNormal2DArrayPlain(sobolNormals, nStocks, nSteps)) {
        goto cleanup_multi_sobol;
    }

    for (step = 0; step < nSteps; step++) {
        for (stock = 0; stock < nStocks; stock++) {
            stepNormals[stock] = sobolNormals[stock * nSteps + step];
        }
        apply_lower_triangular_correlation(choleskyMatrix, nStocks, stepNormals, correlatedNormals + step * nStocks);
    }

    evolve_multi_equity_paths(spots, nStocks, dt, drift, dividend, sigma, correlatedNormals, nSteps, outPaths);
    ok = true;

cleanup_multi_sobol:
    free(stepNormals);
    free(sobolNormals);
    free(correlatedNormals);
    free(choleskyMatrix);
    return ok;
}
