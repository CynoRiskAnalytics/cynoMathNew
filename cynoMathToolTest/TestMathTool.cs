using Microsoft.VisualStudio.TestTools.UnitTesting;
using cynoFinTool.Numerics;

namespace cynoMathToolTest
{
    [TestClass]
    public class TestMathTool
    {
        private const double TightTolerance = 1e-9;
        private const double MediumTolerance = 1e-6;

        private static void AssertAlmostEqual(double expected, double actual, double tolerance = TightTolerance)
        {
            Assert.IsTrue(
                Math.Abs(expected - actual) <= tolerance,
                $"Expected {expected}, actual {actual}, tolerance {tolerance}");
        }

        private static double SquareMinusTwo(double x) => x * x - 2.0;

        private static double Square(double x) => x * x;

        private static double[] BuildExpectedEquityPath(double spot, double[] dt, double[] drift, double[] dividend, double[] sigma, double[] normals)
        {
            double[] path = new double[normals.Length + 1];
            path[0] = spot;
            for (int i = 0; i < normals.Length; i++)
            {
                double variance = sigma[i] * sigma[i];
                double step = (drift[i] - dividend[i] - 0.5 * variance) * dt[i]
                    + sigma[i] * Math.Sqrt(dt[i]) * normals[i];
                path[i + 1] = path[i] * Math.Exp(step);
            }
            return path;
        }

        [TestMethod]
        public void NormalDistributionFunctionsMatchKnownValues()
        {
            AssertAlmostEqual(0.5, cynoMathAPI.cyno_StdNormCDF(0.0));
            AssertAlmostEqual(0.3989422804014327, cynoMathAPI.cyno_StdNormDen(0.0), 1e-12);
            AssertAlmostEqual(0.0, cynoMathAPI.cyno_StdNormInv(0.5), 1e-12);

            AssertAlmostEqual(0.5, cynoMathAPI.cyno_NormCDF(2.0, 2.0, 3.0));
            AssertAlmostEqual(2.0, cynoMathAPI.cyno_NormInv(0.5, 2.0, 3.0), 1e-12);
            AssertAlmostEqual(
                cynoMathAPI.cyno_StdNormDen(1.25) / 3.0,
                cynoMathAPI.cyno_NormDen(2.0 + 3.0 * 1.25, 2.0, 3.0),
                1e-12);
        }

        [TestMethod]
        public void BetaFunctionsAreSelfConsistent()
        {
            const double a = 2.5;
            const double b = 3.5;
            const double x = 0.4;

            double cdf = cynoMathAPI.cyno_BetaCDF(a, b, x);
            double inv = cynoMathAPI.cyno_BetaInv(a, b, cdf);
            double density = cynoMathAPI.cyno_BetaDen(a, b, x);

            Assert.IsTrue(density > 0.0);
            Assert.IsTrue(cdf > 0.0 && cdf < 1.0);
            AssertAlmostEqual(x, inv, MediumTolerance);
        }

        [TestMethod]
        public void TestCallbackAndRootFindingWork()
        {
            int callbackValue = cynoMathAPI.cyno_Test(3.0, x => x * x);
            Assert.AreEqual(9, callbackValue);

            double x1 = 0.0;
            double x2 = 2.0;
            short flag = 0;

            cynoMathAPI.cyno_Root(SquareMinusTwo, ref x1, ref x2, 1.0, 1e-12, 1e-12, ref flag);

            Assert.AreEqual(1, flag);
            AssertAlmostEqual(Math.Sqrt(2.0), x1, 1e-9);
        }

        [TestMethod]
        public void SortingReturnsAscendingIndexOrder()
        {
            double[] values = { 3.0, 1.0, 2.0, -4.0 };
            int[] index = new int[values.Length];

            cynoMathAPI.cyno_Sort(index, values, values.Length);

            double[] ordered = index.Select(i => values[i]).ToArray();
            CollectionAssert.AreEqual(new[] { -4.0, 1.0, 2.0, 3.0 }, ordered);
        }

        [TestMethod]
        public void FiniteAndQuadratureIntegrationsProduceExpectedMoments()
        {
            double element = 0.0;
            cynoMathAPI.cyno_IntegrationFin(Square, 0, 0, ref element, 0.0, 1.0, 1e-12, 1e-12, 100, out double result, out double absErr, out short neval, out short ier);

            AssertAlmostEqual(1.0 / 3.0, result, 1e-10);
            Assert.AreEqual(0, ier);
            Assert.IsTrue(absErr >= 0.0);
            Assert.IsTrue(neval > 0);

            double[] ghX = new double[5];
            double[] ghW = new double[5];
            cynoMathAPI.cyno_IntegrationGHBase(ghX, ghW, 4);
            AssertAlmostEqual(Math.Sqrt(Math.PI), ghW.Skip(1).Sum(), 1e-10);

            double[] ghXFast = new double[5];
            double[] ghWFast = new double[5];
            cynoMathAPI.cyno_IntegrationGH(ghXFast, ghWFast, 4);
            AssertAlmostEqual(Math.Sqrt(Math.PI), ghWFast.Skip(1).Sum(), 1e-10);

            double[] glX = new double[4];
            double[] glW = new double[4];
            cynoMathAPI.cyno_IntegrationGL(glX, glW, 0.0, 4);
            AssertAlmostEqual(1.0, glW.Sum(), 1e-10);
        }

        [TestMethod]
        public void RandomGeneratorsAreSeedableAndFillArrays()
        {
            cynoMathAPI.cyno_Rand_Seed(12345);
            double first = cynoMathAPI.cyno_Rand();
            cynoMathAPI.cyno_Rand_Seed(12345);
            double second = cynoMathAPI.cyno_Rand();
            AssertAlmostEqual(first, second, 1e-15);

            cynoMathAPI.cyno_RandMT_Seed(98765);
            double mtFirst = cynoMathAPI.cyno_MTRand();
            cynoMathAPI.cyno_RandMT_Seed(98765);
            double mtSecond = cynoMathAPI.cyno_MTRand();
            AssertAlmostEqual(mtFirst, mtSecond, 1e-15);

            double[] normals = new double[128];
            bool ok = cynoMathAPI.cyno_NormalRand1DArray(normals, normals.Length);
            Assert.IsTrue(ok);
            Assert.IsTrue(normals.All(x => !double.IsNaN(x) && !double.IsInfinity(x)));

            double[] uniforms = new double[32];
            Assert.IsTrue(cynoMathAPI.cyno_Rand1DArrayPlain(uniforms, uniforms.Length));
            Assert.IsTrue(uniforms.All(x => x >= 0.0 && x <= 1.0));

            double[] antitheticUniforms = new double[6];
            Assert.IsTrue(cynoMathAPI.cyno_Rand1DArrayPlainAntithetic(antitheticUniforms, antitheticUniforms.Length));
            AssertAlmostEqual(1.0, antitheticUniforms[0] + antitheticUniforms[1], 1e-12);
            AssertAlmostEqual(1.0, antitheticUniforms[2] + antitheticUniforms[3], 1e-12);
            AssertAlmostEqual(1.0, antitheticUniforms[4] + antitheticUniforms[5], 1e-12);

            double[] sobol = new double[8];
            Assert.IsTrue(cynoMathAPI.cyno_Sobol1DArrayPlain(sobol, sobol.Length));
            CollectionAssert.AreEqual(
                new[] { 0.0, 0.5, 0.75, 0.25, 0.375, 0.875, 0.625, 0.125 },
                sobol);

            double[] sobolNormals = new double[8];
            Assert.IsTrue(cynoMathAPI.cyno_SobolNormal1DArrayPlain(sobolNormals, sobolNormals.Length));
            Assert.IsTrue(double.IsNegativeInfinity(sobolNormals[0]) == false);
            Assert.IsTrue(sobolNormals.All(x => !double.IsNaN(x) && !double.IsInfinity(x)));

            double[] sobol2D = new double[2 * 4];
            Assert.IsTrue(cynoMathAPI.cyno_Sobol2DArrayPlain(sobol2D, 2, 4));
            CollectionAssert.AreEqual(new[] { 0.0, 0.5, 0.75, 0.25 }, sobol2D.Take(4).ToArray());
            Assert.IsTrue(sobol2D.Skip(4).All(x => x >= 0.0 && x < 1.0));

            double[] sobolHighDim = new double[12 * 4];
            Assert.IsTrue(cynoMathAPI.cyno_Sobol2DArrayPlain(sobolHighDim, 12, 4));
            Assert.IsTrue(sobolHighDim.All(x => x >= 0.0 && x < 1.0));
            Assert.IsTrue(sobolHighDim.Skip(8 * 4).Distinct().Count() > 4);

            double[] sobolNormal2D = new double[2 * 4];
            Assert.IsTrue(cynoMathAPI.cyno_SobolNormal2DArrayPlain(sobolNormal2D, 2, 4));
            Assert.IsTrue(sobolNormal2D.All(x => !double.IsNaN(x) && !double.IsInfinity(x)));

            double[] antitheticNormals = new double[6];
            Assert.IsTrue(cynoMathAPI.cyno_NormalRand1DArrayAntithetic(antitheticNormals, antitheticNormals.Length, 1));
            AssertAlmostEqual(0.0, antitheticNormals[0] + antitheticNormals[1], 1e-12);
            AssertAlmostEqual(0.0, antitheticNormals[2] + antitheticNormals[3], 1e-12);
            AssertAlmostEqual(0.0, antitheticNormals[4] + antitheticNormals[5], 1e-12);

            double[] corr = {
                1.0, 0.0,
                0.0, 1.0
            };
            double[] correlated1D = new double[2];
            Assert.IsTrue(cynoMathAPI.cyno_CorrelatedNormalRand1DArrayPlain(corr, 2, correlated1D, 1));
            Assert.IsTrue(correlated1D.All(x => !double.IsNaN(x) && !double.IsInfinity(x)));

            double[] correlated2D = new double[2 * 4];
            Assert.IsTrue(cynoMathAPI.cyno_CorrelatedNormalRand2DArrayPlain(corr, 2, correlated2D, 4, 1));
            Assert.IsTrue(correlated2D.All(x => !double.IsNaN(x) && !double.IsInfinity(x)));

            double[] correlatedAntithetic = new double[2];
            double[] correlatedAntitheticMirror = new double[2];
            Assert.IsTrue(cynoMathAPI.cyno_CorrelatedNormalRand1DArrayPlainAntithetic(corr, 2, correlatedAntithetic, correlatedAntitheticMirror, 1));
            AssertAlmostEqual(0.0, correlatedAntithetic[0] + correlatedAntitheticMirror[0], 1e-12);
            AssertAlmostEqual(0.0, correlatedAntithetic[1] + correlatedAntitheticMirror[1], 1e-12);

            double[] correlated2DAntithetic = new double[2 * 4];
            Assert.IsTrue(cynoMathAPI.cyno_CorrelatedNormalRand2DArrayPlainAntithetic(corr, 2, correlated2DAntithetic, 4, 1));
            AssertAlmostEqual(0.0, correlated2DAntithetic[0] + correlated2DAntithetic[1], 1e-12);
            AssertAlmostEqual(0.0, correlated2DAntithetic[2] + correlated2DAntithetic[3], 1e-12);
            AssertAlmostEqual(0.0, correlated2DAntithetic[4] + correlated2DAntithetic[5], 1e-12);
            AssertAlmostEqual(0.0, correlated2DAntithetic[6] + correlated2DAntithetic[7], 1e-12);

            for (int i = 0; i < 64; i++)
            {
                double sample = cynoMathAPI.cyno_NormRand();
                double inverseSample = cynoMathAPI.cyno_NormRandInv();
                Assert.IsFalse(double.IsNaN(sample));
                Assert.IsFalse(double.IsInfinity(sample));
                Assert.IsTrue(inverseSample > -10.0 && inverseSample < 10.0);
            }
        }

        [TestMethod]
        public void EquityPathGeneratorsUseStepSpecificInputs()
        {
            const double spot = 100.0;
            double[] dt = { 0.25, 0.5, 0.25 };
            double[] drift = { 0.03, 0.04, 0.05 };
            double[] dividend = { 0.01, 0.015, 0.02 };
            double[] sigma = { 0.20, 0.25, 0.30 };

            cynoMathAPI.cyno_Rand_Seed(24680);
            double[] path = new double[dt.Length + 1];
            Assert.IsTrue(cynoMathAPI.cyno_EquityPathPlain(spot, dt, drift, dividend, sigma, dt.Length, path, 1));

            cynoMathAPI.cyno_Rand_Seed(24680);
            double[] normals = new double[dt.Length];
            Assert.IsTrue(cynoMathAPI.cyno_NormalRand1DArray(normals, normals.Length));
            double[] expected = BuildExpectedEquityPath(spot, dt, drift, dividend, sigma, normals);
            CollectionAssert.AreEqual(expected, path);

            cynoMathAPI.cyno_Rand_Seed(13579);
            double[] antitheticPath = new double[dt.Length + 1];
            double[] mirroredPath = new double[dt.Length + 1];
            Assert.IsTrue(cynoMathAPI.cyno_EquityPathAntitheticPlain(spot, dt, drift, dividend, sigma, dt.Length, antitheticPath, mirroredPath, 1));

            cynoMathAPI.cyno_Rand_Seed(13579);
            double[] antitheticNormals = new double[dt.Length];
            Assert.IsTrue(cynoMathAPI.cyno_NormalRand1DArray(normals, normals.Length));
            for (int i = 0; i < dt.Length; i++)
            {
                antitheticNormals[i] = -normals[i];
            }
            double[] expectedAntithetic = BuildExpectedEquityPath(spot, dt, drift, dividend, sigma, normals);
            double[] expectedMirrored = BuildExpectedEquityPath(spot, dt, drift, dividend, sigma, antitheticNormals);
            CollectionAssert.AreEqual(expectedAntithetic, antitheticPath);
            CollectionAssert.AreEqual(expectedMirrored, mirroredPath);

            double[] sobolPath = new double[dt.Length + 1];
            Assert.IsTrue(cynoMathAPI.cyno_EquityPathSobolPlain(spot, dt, drift, dividend, sigma, dt.Length, sobolPath));
            Assert.IsTrue(sobolPath.All(x => x > 0.0 && !double.IsNaN(x) && !double.IsInfinity(x)));
        }

        [TestMethod]
        public void VanillaBlackScholesFunctionsSatisfyParity()
        {
            const double s = 100.0;
            const double k = 95.0;
            const double tau = 1.0;
            const double td = 1.0;
            const double rd = 0.03;
            const double rf = 0.01;
            const double sigma = 0.2;
            double totalSigma = sigma * Math.Sqrt(tau);

            double call = cynoMathAPI.cyno_Call(s, k, totalSigma);
            double put = cynoMathAPI.cyno_Put(s, k, totalSigma);
            AssertAlmostEqual(k - s + call, put, 1e-12);

            double digitalCall = cynoMathAPI.cyno_DigitalCall(s, k, totalSigma);
            double digitalPut = cynoMathAPI.cyno_DigitalPut(s, k, totalSigma);
            AssertAlmostEqual(1.0, digitalCall + digitalPut, 1e-12);

            double bsCall = cynoMathAPI.cyno_BS(tau, td, s, k, rd, rf, sigma, "CALL");
            double bsPut = cynoMathAPI.cyno_BS(tau, td, s, k, rd, rf, sigma, "PUT");
            AssertAlmostEqual(bsCall, cynoMathAPI.cyno_BS_Call(tau, td, s, k, rd, rf, sigma), 1e-12);
            AssertAlmostEqual(bsPut, cynoMathAPI.cyno_BS_Put(tau, td, s, k, rd, rf, sigma), 1e-12);
            AssertAlmostEqual(s * Math.Exp(-rf * td) - k * Math.Exp(-rd * td), bsCall - bsPut, 1e-9);
        }

        [TestMethod]
        public void SingleBarrierAndTouchFunctionsSatisfyIdentities()
        {
            const double s = 100.0;
            const double k = 100.0;
            const double barrier = 120.0;
            const double tau = 1.0;
            const double td = 1.0;
            const double rd = 0.02;
            const double rf = 0.01;
            const double sigma = 0.25;

            double vanillaCall = cynoMathAPI.cyno_BS(tau, td, s, k, rd, rf, sigma, "CALL");
            double upInCall = cynoMathAPI.cyno_singleBarrier(s, k, barrier, tau, td, rd, rf, sigma, "UPANDIN", "CALL");
            double upOutCall = cynoMathAPI.cyno_singleBarrier(s, k, barrier, tau, td, rd, rf, sigma, "UPANDOUT", "CALL");
            AssertAlmostEqual(vanillaCall, upInCall + upOutCall, 1e-6);

            double oneTouchDom = cynoMathAPI.cyno_oneTouch(tau, td, s, barrier, rd, rf, sigma, "UP", "DOM", "ATMATURITY");
            double oneNoTouchDom = cynoMathAPI.cyno_oneNoTouch(tau, td, s, barrier, rd, rf, sigma, "UP", "DOM", "ATMATURITY");
            AssertAlmostEqual(Math.Exp(-rd * td), oneTouchDom + oneNoTouchDom, 1e-6);

            double hitProbability = cynoMathAPI.cyno_hitProbabilityOneBarrier(tau, td, s, barrier, rd, rf, sigma, "DOM");
            Assert.IsTrue(hitProbability >= 0.0 && hitProbability <= 1.0 + 1e-6);
        }

        [TestMethod]
        public void DoubleBarrierFunctionsSatisfyParityRelations()
        {
            const double s = 100.0;
            const double k = 100.0;
            const double low = 80.0;
            const double high = 120.0;
            const double tau = 1.0;
            const double td = 1.0;
            const double rd = 0.02;
            const double rf = 0.01;
            const double sigma = 0.2;

            double vanillaCall = cynoMathAPI.cyno_BS(tau, td, s, k, rd, rf, sigma, "CALL");
            double knockOutCall = cynoMathAPI.cyno_doubleBarrier(tau, td, s, k, rd, rf, sigma, low, high, 0, 0);
            double knockInCall = cynoMathAPI.cyno_doubleBarrier(tau, td, s, k, rd, rf, sigma, low, high, 1, 0);
            AssertAlmostEqual(vanillaCall, knockOutCall + knockInCall, 1e-6);

            double noTouchDom = cynoMathAPI.cyno_doubleNoTouch(tau, td, s, rd, rf, sigma, low, high, 0);
            double oneTouchDom = cynoMathAPI.cyno_doubleOneTouch(tau, td, s, rd, rf, sigma, low, high, 0);
            AssertAlmostEqual(Math.Exp(-rd * td), noTouchDom + oneTouchDom, 1e-6);

            double probability = cynoMathAPI.cyno_hitProbabilityDoubleBarrier(tau, td, s, rd, rf, sigma, low, high, 0);
            Assert.IsTrue(probability >= 0.0 && probability <= 1.0 + 1e-6);
        }

        [TestMethod]
        public void PlainArrayMatrixApisProduceExpectedResults()
        {
            double[] matrix = {
                1.0, 2.0,
                3.0, 4.0
            };
            double[] scaled = new double[4];
            Assert.IsTrue(cynoMathAPI.cyno_MatrixScalePlain(2.0, matrix, 2, 2, scaled));
            CollectionAssert.AreEqual(new[] { 2.0, 4.0, 6.0, 8.0 }, scaled);

            double[] product = new double[4];
            Assert.IsTrue(cynoMathAPI.cyno_MatrixMultiplyPlain(matrix, 2, 2, matrix, 2, 2, product));
            CollectionAssert.AreEqual(new[] { 7.0, 10.0, 15.0, 22.0 }, product);

            Assert.IsTrue(cynoMathAPI.cyno_VectorDotPlain(new[] { 1.0, 2.0, 3.0 }, new[] { 4.0, 5.0, 6.0 }, 3, out double dot));
            AssertAlmostEqual(32.0, dot);

            double[] transpose = new double[4];
            Assert.IsTrue(cynoMathAPI.cyno_MatrixTransposePlain(matrix, 2, 2, transpose));
            CollectionAssert.AreEqual(new[] { 1.0, 3.0, 2.0, 4.0 }, transpose);

            double[] power = { 1.0, 2.0, 3.0, 4.0 };
            Assert.IsTrue(cynoMathAPI.cyno_MatrixPowerPlain(power, 2, 2));
            CollectionAssert.AreEqual(new[] { 7.0, 10.0, 15.0, 22.0 }, power);

            double[] inverse = new double[4];
            Assert.IsTrue(cynoMathAPI.cyno_InverseMatrixPlain(matrix, 2, inverse));
            AssertAlmostEqual(-2.0, inverse[0], 1e-9);
            AssertAlmostEqual(1.0, inverse[1], 1e-9);
            AssertAlmostEqual(1.5, inverse[2], 1e-9);
            AssertAlmostEqual(-0.5, inverse[3], 1e-9);

            double[] chol = {
                4.0, 2.0,
                2.0, 3.0
            };
            Assert.IsTrue(cynoMathAPI.cyno_CholeskyDecompPlain(chol, 2));
            AssertAlmostEqual(2.0, chol[0], 1e-9);
            AssertAlmostEqual(0.0, chol[1], 1e-9);
            AssertAlmostEqual(1.0, chol[2], 1e-9);
            AssertAlmostEqual(Math.Sqrt(2.0), chol[3], 1e-9);

            double[] eigVals = new double[2];
            double[] eigImag = new double[2];
            double[] eigVecs = new double[4];
            Assert.IsTrue(cynoMathAPI.cyno_EigenSystemPlain(matrix, 2, eigVals, eigImag, eigVecs, true, 0.0, false));
            Assert.IsTrue(eigVals.All(x => !double.IsNaN(x) && !double.IsInfinity(x)));
            Assert.IsTrue(eigImag.All(x => Math.Abs(x) < 1e-9));

            double[] corr = {
                1.0, 0.2,
                0.2, 1.0
            };
            double[] corrPd = (double[])corr.Clone();
            Assert.IsTrue(cynoMathAPI.cyno_PDCorrMatrixPlain(corrPd, 2));
            AssertAlmostEqual(corrPd[1], corrPd[2], 1e-9);

            double[] corrApprox = {
                1.0, 1.2,
                1.2, 1.0
            };
            Assert.IsTrue(cynoMathAPI.cyno_ApproximateCorrMatrixPlain(corrApprox, 2));
            AssertAlmostEqual(corrApprox[1], corrApprox[2], 1e-6);

            double[] rhs = { 5.0, 11.0 };
            double[] solution = new double[2];
            Assert.IsTrue(cynoMathAPI.cyno_LinearSystemSolverPlain(matrix, rhs, 2, solution));
            AssertAlmostEqual(1.0, solution[0], 1e-9);
            AssertAlmostEqual(2.0, solution[1], 1e-9);
        }
    }
}
