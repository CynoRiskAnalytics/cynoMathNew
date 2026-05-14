// Translated from C# to C++
#include "pch.h"
#include "framework.h"
#include <math.h>
#include <string.h>



//namespace cynoFinTool {

    static double cyno_max_double(double a, double b) {
        return (a > b) ? a : b;
    }

    static double cyno_min_double(double a, double b) {
        return (a < b) ? a : b;
    }

    static int cyno_round_to_int(double value) {
        return (int)((value >= 0.0) ? floor(value + 0.5) : ceil(value - 0.5));
    }

    double doubleStdNormCDF(double x) {
        return 0.5 * erfc(-x / sqrt(2.0));
        //comment
    }

    const double PI = 3.14159265358979323846;

    double phiDKO(double z1, double z2, double beta, double gamma, double nu) {

        if (gamma == 0.0) {
			gamma = 1e-12; // Avoid division by zero
		}

        double term1 = 0.5 * nu * (nu * gamma * gamma - 2.0 * beta);
        double arg1 = gamma * nu - (beta + z1) / gamma;
        double arg2 = gamma * nu - (beta + z2) / gamma;
        return exp(term1) * (doubleStdNormCDF(arg1) - doubleStdNormCDF(arg2));
    }

    double FDKO(double tau, double mu, double sigma, double xL, double xH, double nu, double z1, double z2, int nMax) {
        double alpha1 = cyno_max_double(z1, xL);
        double alpha2 = cyno_min_double(z2, xH);

        double res = 0.0;
        for (int i = -nMax; i < nMax; ++i) {
            double u = 2.0 * i * (xH - xL);
            res += exp(-mu * u / (sigma * sigma)) * phiDKO(alpha1, alpha2, -mu * tau + u, sigma * sqrt(tau), nu)
                - exp(mu * (2.0 * xH - u) / (sigma * sigma)) * phiDKO(alpha1, alpha2, -mu * tau - 2.0 * xH + u, sigma * sqrt(tau), nu);
        }
        return res;
    }

    double BS(double Tau, double TD, double S, double K, double rD, double rF, double sigma, int callPutIn)
    {
        //call=0 put=1

        //if (Tau <= 0.0 || sigma <= 0.0) return std::max((callPut == "CALL") ? S - K : K - S, 0.0);
        if (Tau <= 0.0 || sigma <= 0.0) return cyno_max_double((callPutIn ==0) ? S - K : K - S, 0.0);


        double sqrtTau = sqrt(Tau);
        double d1 = (log(S / K) + (rD - rF) * TD + 0.5 * sigma * sigma * Tau) / (sigma * sqrtTau);
        double d2 = d1 - sigma * sqrtTau;

        if (callPutIn == 0) {
            return S * exp(-rF * TD) * doubleStdNormCDF(d1) - K * exp(-rD * TD) * doubleStdNormCDF(d2);
        }
        else {
            return K * exp(-rD * TD) * doubleStdNormCDF(-d2) - S * exp(-rF * TD) * doubleStdNormCDF(-d1);
        }
    }

    CYNOMATHUTILITY_API double __stdcall cyno_doubleBarrier(double Tau, double TD, double S, double K, double rD, double rF, double sigma,
        double L, double H, int knockTypeIn, int callPutIn) {

        //if (S >= H) throw std::invalid_argument("Higher barrier has to be higher than spot");
        //if (S <= L) throw std::invalid_argument("Lower barrier has to be lower than spot");
        if (S / L > 5) L = S / 5.0;
        if (H / S > 5) H = 5.0 * S;


       // std::string knockType = knockTypeIn;
       // std::string callPut = callPutIn;

        //std::transform(knockType.begin(), knockType.end(), knockType.begin(), ::toupper);
        //std::transform(callPut.begin(), callPut.end(), callPut.begin(), ::toupper);

		//knockout=0, knockin=1
		//call=0, put=1

        double muh = (rD - rF) * TD / Tau - 0.5 * sigma * sigma;
        double mup = sqrt(muh * muh + 2.0 * rD * sigma * sigma);

        double x0 = log(K / S);
        double xL = log(L / S);
        double xH = log(H / S);
        double eps = 1e-8;

        int nMax = cyno_round_to_int(sqrt(cyno_max_double(25.0, (-2 * log(eps) / Tau - pow(muh / sigma, 2)) * log(H / L) / PI / sigma)));

        double res = 0.0;

       // return BS(Tau, TD, S, K, rD, rF, sigma, 0);

      //  res= cyno_doubleBarrier(Tau, TD, S, K, rD, rF, sigma, L, H, 0, callPutIn);

        //return res;

        if (knockTypeIn == 0) {
            if (callPutIn == 0) {
                res = exp(-rD * TD) * (
                    S * FDKO(Tau, muh, sigma, xL, xH, 1, x0, xH, nMax) -
                    K * FDKO(Tau, muh, sigma, xL, xH, 0, x0, xH, nMax));
                res = cyno_max_double(0.0, res);
            }
            else {
                res = exp(-rD * TD) * (
                    K * FDKO(Tau, muh, sigma, xL, xH, 0, xL, x0, nMax) -
                    S * FDKO(Tau, muh, sigma, xL, xH, 1, xL, x0, nMax));
                res = cyno_max_double(0.0, res);		
            }
        }
        else if (knockTypeIn == 1) {
            if (callPutIn == 0) {
                res = BS(Tau, TD, S, K, rD, rF, sigma, 0) -
                    cyno_doubleBarrier(Tau, TD, S, K, rD, rF, sigma, L, H,0, callPutIn);
		
            }
            else {
                res = BS(Tau, TD, S, K, rD, rF, sigma, 1) -
                    cyno_doubleBarrier(Tau, TD, S, K, rD, rF, sigma, L, H, 0, callPutIn);
			
            }
        }
        return res;
    }






    CYNOMATHUTILITY_API double __stdcall cyno_doubleNoTouch(double Tau, double TD, double S, double rD, double rF, double sigma, double L, double H, int rCurrIn) 
    {
		//DOM=0, FOR=1

        if (rCurrIn != 0 && rCurrIn != 1) {
            return 0.0;
        }

        //std::transform(rCurr.begin(), rCurr.end(), rCurr.begin(), ::toupper);

       // if (S >= H || S <= L) throw std::invalid_argument("Spot has to be between lower barrier and higher barrier");

       // return 0;

        double muh = (rD - rF) * TD / Tau - 0.5 * sigma * sigma;
        double mup = sqrt(muh * muh + 2.0 * rD * sigma * sigma);
        double xL = log(L / S);
        double xH = log(H / S);
        double eps = 1e-8;

        int nMax = (int)(sqrt(cyno_max_double(25.0, (-2 * log(eps) / Tau - pow(muh / sigma, 2)) * log(H / L) / PI / sigma)));

        if (S / L > 5.0) L = S / 5.0;
        if (H / S > 5.0) H = 5.0 * S;

        double res = exp(-rD * TD) * FDKO(Tau, muh, sigma, xL, xH, 0, xL, xH, nMax);

        if (rCurrIn == 1) {
            res = L * cyno_doubleNoTouch(Tau, TD, S, rD, rF, sigma, L, H, 0) +
                cyno_doubleBarrier(Tau, TD, S, L, rD, rF, sigma, L, H, 0, 0);
        }

        return res;
    }

    double cyno_doubleOneTouch(double Tau, double TD, double S, double rD, double rF, double sigma, double L, double H, int rCurrIn) {
        if (rCurrIn != 0 && rCurrIn != 1) {
            return 0.0;
        }
        //std::transform(rCurr.begin(), rCurr.end(), rCurr.begin(), ::toupper);

        if (rCurrIn == 1) {
            return exp(-rF * TD) * S - cyno_doubleNoTouch(Tau, TD, S, rD, rF, sigma, L, H, rCurrIn);
        }
        else {
            return exp(-rD * TD) - cyno_doubleNoTouch(Tau, TD, S, rD, rF, sigma, L, H, rCurrIn);
        }
    }

    double cyno_hitProbabilityDoubleBarrier(double Tau, double TD, double S, double rD, double rF, double sigma, double L, double H, int nCurrencyIn) {
        if (nCurrencyIn != 0 && nCurrencyIn != 1) {
            return 0.0;
        }
        //std::transform(nCurrency.begin(), nCurrency.end(), nCurrency.begin(), ::toupper);

        if (S < H && S > L) {
            double res = exp(rD * TD) * cyno_doubleOneTouch(Tau, TD, S, rD, rF, sigma, L, H, nCurrencyIn);
            if (nCurrencyIn == 1) res /= S;
            return res;
        }
        else {
            return 0.0;
        }
    }


//} // namespace cynoFinTool
