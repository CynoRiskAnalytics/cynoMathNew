using System.Runtime.InteropServices;

internal static class CynoMathApi
{
#if WINDOWS
    private const string NativeLibraryName = "cynoMathTools.dll";
#else
    private const string NativeLibraryName = "libcynomath.so";
#endif

    [DllImport(NativeLibraryName)]
    internal static extern double cyno_BS(double tau, double td, double s, double k, double rd, double rf, double sigma, string callPut);

    [DllImport(NativeLibraryName)]
    internal static extern bool cyno_Rand1DArrayPlain([Out] double[] output, long size);

    [DllImport(NativeLibraryName)]
    internal static extern bool cyno_Sobol1DArrayPlain([Out] double[] output, long size);

    [DllImport(NativeLibraryName)]
    internal static extern bool cyno_NormalRand1DArrayAntithetic([Out] double[] output, int size, short randflag);

    [DllImport(NativeLibraryName)]
    internal static extern bool cyno_MatrixMultiplyPlain(
        [In] double[] a, long aRows, long aCols,
        [In] double[] b, long bRows, long bCols,
        [Out] double[] output);
}

internal static class Program
{
    private static void Main()
    {
        double callPrice = CynoMathApi.cyno_BS(1.0, 1.0, 100.0, 95.0, 0.03, 0.01, 0.20, "CALL");
        Console.WriteLine($"Black-Scholes call price: {callPrice:F10}");

        double[] uniforms = new double[5];
        bool randomOk = CynoMathApi.cyno_Rand1DArrayPlain(uniforms, uniforms.Length);
        Console.WriteLine($"Random fill ok: {randomOk}");
        for (int i = 0; i < uniforms.Length; i++)
        {
            Console.WriteLine($"uniform[{i}] = {uniforms[i]:F10}");
        }

        double[] sobol = new double[8];
        bool sobolOk = CynoMathApi.cyno_Sobol1DArrayPlain(sobol, sobol.Length);
        Console.WriteLine($"Sobol fill ok: {sobolOk}");
        for (int i = 0; i < sobol.Length; i++)
        {
            Console.WriteLine($"sobol[{i}] = {sobol[i]:F10}");
        }

        double[] antitheticNormals = new double[6];
        bool antitheticOk = CynoMathApi.cyno_NormalRand1DArrayAntithetic(antitheticNormals, antitheticNormals.Length, 1);
        Console.WriteLine($"Antithetic normal fill ok: {antitheticOk}");
        for (int i = 0; i < antitheticNormals.Length; i += 2)
        {
            Console.WriteLine($"pair {i / 2}: {antitheticNormals[i]:F10}, {antitheticNormals[i + 1]:F10}");
        }

        double[] left =
        {
            1.0, 2.0,
            3.0, 4.0
        };
        double[] right =
        {
            5.0, 6.0,
            7.0, 8.0
        };
        double[] product = new double[4];

        bool matrixOk = CynoMathApi.cyno_MatrixMultiplyPlain(left, 2, 2, right, 2, 2, product);
        Console.WriteLine($"Matrix multiply ok: {matrixOk}");
        Console.WriteLine($"[{product[0]:F1} {product[1]:F1}]");
        Console.WriteLine($"[{product[2]:F1} {product[3]:F1}]");
    }
}
