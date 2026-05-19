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
    internal static extern bool cyno_Sobol2DArrayPlain([Out] double[] output, long dim, long samplesPerFactor);

    [DllImport(NativeLibraryName)]
    internal static extern bool cyno_CorrelatedNormalRand2DArrayPlain([In] double[] correlationMatrix, long dim, [Out] double[] output, long samplesPerFactor, short randflag);

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

        double[,] corrMatrix2D =
        {
            { 1.0, 0.4 },
            { 0.4, 1.0 }
        };
        double[] corr = FlattenMatrix(corrMatrix2D);

        double[] sobol2DFlat = new double[2 * 4];
        bool sobol2DOk = CynoMathApi.cyno_Sobol2DArrayPlain(sobol2DFlat, 2, 4);
        Console.WriteLine($"Sobol 2D fill ok: {sobol2DOk}");
        PrintMatrix("Sobol 2D uniforms", ToMatrix(sobol2DFlat, 2, 4));

        double[] correlated2DFlat = new double[2 * 4];
        bool correlated2DOk = CynoMathApi.cyno_CorrelatedNormalRand2DArrayPlain(corr, 2, correlated2DFlat, 4, 1);
        Console.WriteLine($"Correlated normal 2D fill ok: {correlated2DOk}");
        PrintMatrix("Correlation matrix loaded into C", corrMatrix2D);
        PrintMatrix("Correlated normal 2D", ToMatrix(correlated2DFlat, 2, 4));

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

    private static double[] FlattenMatrix(double[,] matrix)
    {
        int rows = matrix.GetLength(0);
        int cols = matrix.GetLength(1);
        double[] flat = new double[rows * cols];

        for (int row = 0; row < rows; row++)
        {
            for (int col = 0; col < cols; col++)
            {
                flat[row * cols + col] = matrix[row, col];
            }
        }

        return flat;
    }

    private static double[,] ToMatrix(double[] flat, int rows, int cols)
    {
        double[,] matrix = new double[rows, cols];

        for (int row = 0; row < rows; row++)
        {
            for (int col = 0; col < cols; col++)
            {
                matrix[row, col] = flat[row * cols + col];
            }
        }

        return matrix;
    }

    private static void PrintMatrix(string title, double[,] matrix)
    {
        Console.WriteLine(title);
        int rows = matrix.GetLength(0);
        int cols = matrix.GetLength(1);

        for (int row = 0; row < rows; row++)
        {
            for (int col = 0; col < cols; col++)
            {
                Console.Write($"{matrix[row, col]:F10} ");
            }

            Console.WriteLine();
        }
    }
}
