using System;
using System.Runtime.InteropServices;

namespace cynoFinTool.Numerics
{
    [UnmanagedFunctionPointer(CallingConvention.StdCall)]
    public delegate double myDel(double Value);
    static public class cynoMathAPI
    {
#if WINDOWS
        private const string NativeLibraryName = "cynoMathTools.dll";
#else
        private const string NativeLibraryName = "libcynomath.so";
#endif

        //1. Normal functions
        [DllImport(NativeLibraryName)]
        public static extern double cyno_StdNormCDF(double x);

        [DllImport(NativeLibraryName)]
        public static extern double cyno_NormCDF(double x, double mean, double sig);

        [DllImport(NativeLibraryName)]
        public static extern double cyno_StdNormDen(double x);

        [DllImport(NativeLibraryName)]
        public static extern double cyno_NormDen(double x, double mean, double sig);

        [DllImport(NativeLibraryName)]
        public static extern double cyno_StdNormInv(double x);

        [DllImport(NativeLibraryName)]
        public static extern double cyno_NormInv(double x, double mean, double sig);

        [DllImport(NativeLibraryName)]
        public static extern double cyno_BetaDen(double a, double b, double x);

        [DllImport(NativeLibraryName)]
        public static extern double cyno_BetaCDF(double a, double b, double x);

        [DllImport(NativeLibraryName)]
        public static extern double cyno_BetaInv(double a, double b, double p);

        //3. Random number
        [DllImport(NativeLibraryName)]
        public extern static void cyno_Rand_Seed(int ByVal);


        [DllImport(NativeLibraryName)]
        public extern static int cyno_Test(double byval, myDel f);

        [DllImport(NativeLibraryName)]
        public extern static double cyno_Rand();

        [DllImport(NativeLibraryName)]
        public extern static double cyno_MTRand();

        [DllImport(NativeLibraryName)]
        public extern static void cyno_RandMT_Seed(int seed);

        [DllImport(NativeLibraryName)]
        public static extern double cyno_NormRand();

        [DllImport(NativeLibraryName)]
        public static extern bool cyno_Rand1DArrayPlainAntithetic([Out] double[] outArray, long size);

        /// <summary>
        /// generate normal random numbers in a 1D array
        /// </summary>
        /// <param name="rand_array">The output array with random numbers</param>
        /// <returns>true if the call was successful</returns>
        [DllImport(NativeLibraryName)]
        public static extern bool cyno_NormalRand1DArray([Out] double[] rand_array, int size);

        [DllImport(NativeLibraryName)]
        public static extern bool cyno_NormalRand1DArrayAntithetic([Out] double[] randArray, int size, short randflag);

        [DllImport(NativeLibraryName)]
        public static extern bool cyno_Rand1DArrayPlain([Out] double[] outArray, long size);

        [DllImport(NativeLibraryName)]
        public static extern bool cyno_Sobol1DArrayPlain([Out] double[] outArray, long size);

        [DllImport(NativeLibraryName)]
        public static extern bool cyno_SobolNormal1DArrayPlain([Out] double[] outArray, long size);

        [DllImport(NativeLibraryName)]
        public static extern bool cyno_Sobol2DArrayPlain([Out] double[] outArray, long dim, long samplesPerFactor);

        [DllImport(NativeLibraryName)]
        public static extern bool cyno_SobolNormal2DArrayPlain([Out] double[] outArray, long dim, long samplesPerFactor);

        [DllImport(NativeLibraryName)]
        public static extern bool cyno_CorrelatedNormalRand1DArrayPlain([In] double[] correlationMatrix, long dim, [Out] double[] outArray, short randflag);

        [DllImport(NativeLibraryName)]
        public static extern bool cyno_CorrelatedNormalRand2DArrayPlain([In] double[] correlationMatrix, long dim, [Out] double[] outArray, long samplesPerFactor, short randflag);

        [DllImport(NativeLibraryName)]
        public static extern bool cyno_CorrelatedNormalRand1DArrayPlainAntithetic([In] double[] correlationMatrix, long dim, [Out] double[] outArray, [Out] double[] antitheticOutArray, short randflag);

        [DllImport(NativeLibraryName)]
        public static extern bool cyno_CorrelatedNormalRand2DArrayPlainAntithetic([In] double[] correlationMatrix, long dim, [Out] double[] outArray, long samplesPerFactor, short randflag);

        [DllImport(NativeLibraryName)]
        public extern static double cyno_NormRandInv();

        [DllImport(NativeLibraryName)]
        public static extern bool cyno_MatrixScalePlain(double a, [In] double[] input, long rows, long cols, [Out] double[] output);

        [DllImport(NativeLibraryName)]
        public static extern bool cyno_MatrixMultiplyPlain([In] double[] a, long aRows, long aCols, [In] double[] b, long bRows, long bCols, [Out] double[] output);

        [DllImport(NativeLibraryName)]
        public static extern bool cyno_VectorDotPlain([In] double[] a, [In] double[] b, long size, out double output);

        [DllImport(NativeLibraryName)]
        public static extern bool cyno_InverseMatrixPlain([In] double[] input, long dim, [Out] double[] output);

        [DllImport(NativeLibraryName)]
        public static extern bool cyno_MatrixTransposePlain([In] double[] input, long rows, long cols, [Out] double[] output);

        [DllImport(NativeLibraryName)]
        public static extern bool cyno_MatrixPowerPlain([In, Out] double[] inOutMatrix, long dim, short power);

        [DllImport(NativeLibraryName)]
        public static extern bool cyno_CholeskyDecompPlain([In, Out] double[] inOutMatrix, long dim);

        [DllImport(NativeLibraryName)]
        public static extern bool cyno_EigenSystemPlain([In] double[] input, long dim, [Out] double[] vr, [Out] double[] vi, [Out] double[] eigenVectors, bool computeVects, double symTol, bool sortVects);

        [DllImport(NativeLibraryName)]
        public static extern bool cyno_PDCorrMatrixPlain([In, Out] double[] inOutMatrix, long dim);

        [DllImport(NativeLibraryName)]
        public static extern bool cyno_ApproximateCorrMatrixPlain([In, Out] double[] inOutMatrix, long dim);

        [DllImport(NativeLibraryName)]
        public static extern bool cyno_LinearSystemSolverPlain([In] double[] jacobian, [In] double[] rhs, long dim, [Out] double[] output);

        //4. Sorting
        [DllImport(NativeLibraryName)]
        public extern static int cyno_Sort(int[] index, double[] x, int arraySize);

        //5. Root searching
        [DllImport(NativeLibraryName)]
        public extern static void cyno_Root(myDel f, ref double x1, ref double x2, double xguess, double RelError, double AbsError, ref short iFlag);

        //6. Integration(Finite Domain)
        [DllImport(NativeLibraryName)]
        public extern static void cyno_IntegrationFin(myDel f, short minIndex, short maxIndex, ref double element, double a, double b, double EPSABS, double EPSREL, short LIMIT, out double result, out double ABSERR, out short neval, out short ier);

        //7. Integration(Infinite Domain - Gaussian Hermite)
        [DllImport(NativeLibraryName)]
        public extern static void cyno_IntegrationGHBase([Out] double[] x, [Out] double[] w, short npoints);

        [DllImport(NativeLibraryName)]
        public extern static void cyno_IntegrationGH([Out] double[] x, [Out] double[] w, short npoints);

        //8. Integration( Half Infinite Domain - Gaussian-Laguerre)
        [DllImport(NativeLibraryName)]
        public extern static void cyno_IntegrationGL([Out] double[] x, [Out] double[] w, double inputAlpha, short npoints);

        [DllImport(NativeLibraryName)]
        public static extern double cyno_Call(double S, double K, double totalSigma);

        [DllImport(NativeLibraryName)]
        public static extern double cyno_Put(double S, double K, double totalSigma);

        [DllImport(NativeLibraryName)]
        public static extern double cyno_DigitalCall(double S, double K, double totalSigma);

        [DllImport(NativeLibraryName)]
        public static extern double cyno_DigitalPut(double S, double K, double totalSigma);

        [DllImport(NativeLibraryName)]
        public static extern double cyno_BS(double Tau, double TD, double S, double K, double rD, double rF, double sigma, string callPut);

        [DllImport(NativeLibraryName)]
        public static extern double cyno_BS_Call(double Tau, double TD, double S, double K, double rD, double rF, double sigma);

        [DllImport(NativeLibraryName)]
        public static extern double cyno_BS_Put(double Tau, double TD, double S, double K, double rD, double rF, double sigma);

        [DllImport(NativeLibraryName)]
        public static extern double cyno_singleBarrier(double S, double K, double B, double Tau, double TD, double rd, double rf, double sigma, string knockType, string callput);

        [DllImport(NativeLibraryName)]
        public static extern double cyno_oneTouch(double Tau, double TD, double S, double B, double rd, double rf, double sigma, string upDown, string nCurrency, string nPaidAt);

        [DllImport(NativeLibraryName)]
        public static extern double cyno_hitProbabilityOneBarrier(double Tau, double TD, double S, double H, double rD, double rF, double sigma, string nCurrency);

        [DllImport(NativeLibraryName)]
        public static extern double cyno_oneNoTouch(double Tau, double TD, double S, double B, double rd, double rf, double sigma, string upDown, string nCurrency, string nPaidAt);

        [DllImport(NativeLibraryName)]
        public static extern double cyno_doubleBarrier(double Tau, double TD, double S, double K, double rD, double rF, double sigma, double L, double H, int knockTypeIn, int callPutIn);

        [DllImport(NativeLibraryName)]
        public static extern double cyno_doubleNoTouch(double Tau, double TD, double S, double rD, double rF, double sigma, double L, double H, int rCurrIn);

        [DllImport(NativeLibraryName)]
        public static extern double cyno_doubleOneTouch(double Tau, double TD, double S, double rD, double rF, double sigma, double L, double H, int rCurrIn);

        [DllImport(NativeLibraryName)]
        public static extern double cyno_hitProbabilityDoubleBarrier(double Tau, double TD, double S, double rD, double rF, double sigma, double L, double H, int nCurrencyIn);

    }
}
 
