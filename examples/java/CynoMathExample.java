import com.sun.jna.Library;
import com.sun.jna.Native;

public final class CynoMathExample {
    private interface CynoMathLibrary extends Library {
        CynoMathLibrary INSTANCE = Native.load("cynomath", CynoMathLibrary.class);

        double cyno_BS(double tau, double td, double s, double k, double rd, double rf, double sigma, String callPut);
        boolean cyno_Rand1DArrayPlain(double[] outArray, long size);
        boolean cyno_Sobol2DArrayPlain(double[] outArray, long dim, long samplesPerFactor);
        boolean cyno_CorrelatedNormalRand2DArrayPlain(double[] correlationMatrix, long dim, double[] outArray, long samplesPerFactor, short randflag);
        boolean cyno_Sobol1DArrayPlain(double[] outArray, long size);
        boolean cyno_NormalRand1DArrayAntithetic(double[] outArray, int size, short randflag);
        boolean cyno_MatrixMultiplyPlain(
            double[] a, long aRows, long aCols,
            double[] b, long bRows, long bCols,
            double[] output);
    }

    public static void main(String[] args) {
        CynoMathLibrary lib = CynoMathLibrary.INSTANCE;

        double callPrice = lib.cyno_BS(1.0, 1.0, 100.0, 95.0, 0.03, 0.01, 0.20, "CALL");
        System.out.println("Black-Scholes call price: " + callPrice);

        double[] uniforms = new double[5];
        boolean randomOk = lib.cyno_Rand1DArrayPlain(uniforms, uniforms.length);
        System.out.println("Random fill ok: " + randomOk);
        for (int i = 0; i < uniforms.length; i++) {
            System.out.printf("uniform[%d] = %.10f%n", i, uniforms[i]);
        }

        double[][] corrMatrix2D = {
            { 1.0, 0.4 },
            { 0.4, 1.0 }
        };
        double[] corr = flattenMatrix(corrMatrix2D);

        double[] sobol2D = new double[2 * 4];
        boolean sobol2DOk = lib.cyno_Sobol2DArrayPlain(sobol2D, 2, 4);
        System.out.println("Sobol 2D fill ok: " + sobol2DOk);
        printMatrix("Sobol 2D uniforms", toMatrix(sobol2D, 2, 4));

        double[] correlated2D = new double[2 * 4];
        boolean correlated2DOk = lib.cyno_CorrelatedNormalRand2DArrayPlain(corr, 2, correlated2D, 4, (short)1);
        System.out.println("Correlated normal 2D fill ok: " + correlated2DOk);
        printMatrix("Correlation matrix loaded into C", corrMatrix2D);
        printMatrix("Correlated normal 2D", toMatrix(correlated2D, 2, 4));

        double[] sobol = new double[8];
        boolean sobolOk = lib.cyno_Sobol1DArrayPlain(sobol, sobol.length);
        System.out.println("Sobol fill ok: " + sobolOk);
        for (int i = 0; i < sobol.length; i++) {
            System.out.printf("sobol[%d] = %.10f%n", i, sobol[i]);
        }

        double[] antitheticNormals = new double[6];
        boolean antitheticOk = lib.cyno_NormalRand1DArrayAntithetic(antitheticNormals, antitheticNormals.length, (short)1);
        System.out.println("Antithetic normal fill ok: " + antitheticOk);
        for (int i = 0; i < antitheticNormals.length; i += 2) {
            System.out.printf("pair %d = %.10f, %.10f%n", i / 2, antitheticNormals[i], antitheticNormals[i + 1]);
        }

        double[] left = {
            1.0, 2.0,
            3.0, 4.0
        };
        double[] right = {
            5.0, 6.0,
            7.0, 8.0
        };
        double[] product = new double[4];

        boolean matrixOk = lib.cyno_MatrixMultiplyPlain(left, 2, 2, right, 2, 2, product);
        System.out.println("Matrix multiply ok: " + matrixOk);
        System.out.printf("[%.1f %.1f]%n", product[0], product[1]);
        System.out.printf("[%.1f %.1f]%n", product[2], product[3]);
    }

    private static double[] flattenMatrix(double[][] matrix) {
        int rows = matrix.length;
        int cols = matrix[0].length;
        double[] flat = new double[rows * cols];

        for (int row = 0; row < rows; row++) {
            for (int col = 0; col < cols; col++) {
                flat[row * cols + col] = matrix[row][col];
            }
        }

        return flat;
    }

    private static double[][] toMatrix(double[] flat, int rows, int cols) {
        double[][] matrix = new double[rows][cols];

        for (int row = 0; row < rows; row++) {
            for (int col = 0; col < cols; col++) {
                matrix[row][col] = flat[row * cols + col];
            }
        }

        return matrix;
    }

    private static void printMatrix(String title, double[][] matrix) {
        System.out.println(title);
        for (double[] row : matrix) {
            for (double value : row) {
                System.out.printf("%.10f ", value);
            }
            System.out.println();
        }
    }
}
