import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

public final class CynoMathExample {
    static {
        loadNativeLibraries();
    }

    private CynoMathExample() {
    }

    private static void loadNativeLibraries() {
        try {
            Path repoRoot = detectRepoRoot();
            System.load(repoRoot.resolve("build/libcynomath.so").toString());
            System.load(repoRoot.resolve("examples/java/build/libcynomathjni.so").toString());
        } catch (Exception ex) {
            throw new RuntimeException("Failed to load JNI libraries", ex);
        }
    }

    private static Path detectRepoRoot() throws Exception {
        Path codeSource = Paths.get(CynoMathExample.class.getProtectionDomain().getCodeSource().getLocation().toURI());
        Path candidate = codeSource.toAbsolutePath().resolve("../../../../").normalize();
        if (Files.exists(candidate.resolve("build/libcynomath.so"))) {
            return candidate;
        }
        return Paths.get("").toAbsolutePath();
    }

    private static native double cynoBS(
        double tau,
        double td,
        double s,
        double k,
        double rd,
        double rf,
        double sigma,
        String callPut);

    private static native boolean cynoRand1DArrayPlain(double[] outArray);

    private static native boolean cynoSobol1DArrayPlain(double[] outArray);

    private static native boolean cynoNormalRand1DArrayAntithetic(double[] outArray, short randflag);

    private static native boolean cynoMatrixMultiplyPlain(
        double[] a,
        int aRows,
        int aCols,
        double[] b,
        int bRows,
        int bCols,
        double[] output);

    public static void main(String[] args) {
        double callPrice = cynoBS(1.0, 1.0, 100.0, 95.0, 0.03, 0.01, 0.20, "CALL");
        System.out.println("Black-Scholes call price: " + callPrice);

        double[] uniforms = new double[5];
        boolean randomOk = cynoRand1DArrayPlain(uniforms);
        System.out.println("Random fill ok: " + randomOk);
        for (int i = 0; i < uniforms.length; i++) {
            System.out.printf("uniform[%d] = %.10f%n", i, uniforms[i]);
        }

        double[] sobol = new double[8];
        boolean sobolOk = cynoSobol1DArrayPlain(sobol);
        System.out.println("Sobol fill ok: " + sobolOk);
        for (int i = 0; i < sobol.length; i++) {
            System.out.printf("sobol[%d] = %.10f%n", i, sobol[i]);
        }

        double[] antitheticNormals = new double[6];
        boolean antitheticOk = cynoNormalRand1DArrayAntithetic(antitheticNormals, (short)1);
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

        boolean matrixOk = cynoMatrixMultiplyPlain(left, 2, 2, right, 2, 2, product);
        System.out.println("Matrix multiply ok: " + matrixOk);
        System.out.printf("[%.1f %.1f]%n", product[0], product[1]);
        System.out.printf("[%.1f %.1f]%n", product[2], product[3]);
    }
}
