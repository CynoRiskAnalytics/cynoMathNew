#include <jni.h>
#include <stdbool.h>

#include "../../cynoMathTool/framework.h"

JNIEXPORT jdouble JNICALL Java_CynoMathExample_cynoBS(
    JNIEnv* env,
    jclass clazz,
    jdouble tau,
    jdouble td,
    jdouble s,
    jdouble k,
    jdouble rd,
    jdouble rf,
    jdouble sigma,
    jstring callPut) {
    const char* callPutChars;
    double result;

    (void)clazz;

    if (callPut == NULL) {
        return 0.0;
    }

    callPutChars = (*env)->GetStringUTFChars(env, callPut, NULL);
    if (callPutChars == NULL) {
        return 0.0;
    }

    result = cyno_BS(tau, td, s, k, rd, rf, sigma, callPutChars);
    (*env)->ReleaseStringUTFChars(env, callPut, callPutChars);
    return result;
}

JNIEXPORT jboolean JNICALL Java_CynoMathExample_cynoRand1DArrayPlain(
    JNIEnv* env,
    jclass clazz,
    jdoubleArray outArray) {
    jdouble* values;
    jsize size;
    bool ok;

    (void)clazz;

    if (outArray == NULL) {
        return JNI_FALSE;
    }

    size = (*env)->GetArrayLength(env, outArray);
    values = (*env)->GetDoubleArrayElements(env, outArray, NULL);
    if (values == NULL) {
        return JNI_FALSE;
    }

    ok = cyno_Rand1DArrayPlain(values, (long)size);
    (*env)->ReleaseDoubleArrayElements(env, outArray, values, 0);
    return ok ? JNI_TRUE : JNI_FALSE;
}

JNIEXPORT jboolean JNICALL Java_CynoMathExample_cynoSobol1DArrayPlain(
    JNIEnv* env,
    jclass clazz,
    jdoubleArray outArray) {
    jdouble* values;
    jsize size;
    bool ok;

    (void)clazz;

    if (outArray == NULL) {
        return JNI_FALSE;
    }

    size = (*env)->GetArrayLength(env, outArray);
    values = (*env)->GetDoubleArrayElements(env, outArray, NULL);
    if (values == NULL) {
        return JNI_FALSE;
    }

    ok = cyno_Sobol1DArrayPlain(values, (long)size);
    (*env)->ReleaseDoubleArrayElements(env, outArray, values, 0);
    return ok ? JNI_TRUE : JNI_FALSE;
}

JNIEXPORT jboolean JNICALL Java_CynoMathExample_cynoNormalRand1DArrayAntithetic(
    JNIEnv* env,
    jclass clazz,
    jdoubleArray outArray,
    jshort randflag) {
    jdouble* values;
    jsize size;
    bool ok;

    (void)clazz;

    if (outArray == NULL) {
        return JNI_FALSE;
    }

    size = (*env)->GetArrayLength(env, outArray);
    values = (*env)->GetDoubleArrayElements(env, outArray, NULL);
    if (values == NULL) {
        return JNI_FALSE;
    }

    ok = cyno_NormalRand1DArrayAntithetic(values, (int)size, randflag);
    (*env)->ReleaseDoubleArrayElements(env, outArray, values, 0);
    return ok ? JNI_TRUE : JNI_FALSE;
}

JNIEXPORT jboolean JNICALL Java_CynoMathExample_cynoMatrixMultiplyPlain(
    JNIEnv* env,
    jclass clazz,
    jdoubleArray a,
    jint aRows,
    jint aCols,
    jdoubleArray b,
    jint bRows,
    jint bCols,
    jdoubleArray output) {
    jdouble* aValues;
    jdouble* bValues;
    jdouble* outValues;
    bool ok;

    (void)clazz;

    if (a == NULL || b == NULL || output == NULL) {
        return JNI_FALSE;
    }

    aValues = (*env)->GetDoubleArrayElements(env, a, NULL);
    bValues = (*env)->GetDoubleArrayElements(env, b, NULL);
    outValues = (*env)->GetDoubleArrayElements(env, output, NULL);
    if (aValues == NULL || bValues == NULL || outValues == NULL) {
        if (aValues != NULL) {
            (*env)->ReleaseDoubleArrayElements(env, a, aValues, JNI_ABORT);
        }
        if (bValues != NULL) {
            (*env)->ReleaseDoubleArrayElements(env, b, bValues, JNI_ABORT);
        }
        if (outValues != NULL) {
            (*env)->ReleaseDoubleArrayElements(env, output, outValues, 0);
        }
        return JNI_FALSE;
    }

    ok = cyno_MatrixMultiplyPlain(aValues, aRows, aCols, bValues, bRows, bCols, outValues);

    (*env)->ReleaseDoubleArrayElements(env, a, aValues, JNI_ABORT);
    (*env)->ReleaseDoubleArrayElements(env, b, bValues, JNI_ABORT);
    (*env)->ReleaseDoubleArrayElements(env, output, outValues, 0);
    return ok ? JNI_TRUE : JNI_FALSE;
}
