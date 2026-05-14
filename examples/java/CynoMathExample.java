import com.sun.jna.Library;
import com.sun.jna.Native;

public final class CynoMathExample {
    private interface CynoMathLibrary extends Library {
        CynoMathLibrary INSTANCE = Native.load("cynomath", CynoMathLibrary.class);

        double cyno_BS(double tau, double td, double s, double k, double rd, double rf, double sigma, String callPut);
        boolean cyno_Rand1DArrayPlain(double[] outArray, long size);
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
}
