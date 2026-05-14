#include "pch.h"
#include "framework.h"
#include "BlackScholes.h"
#include <math.h>
#include <string.h>



//namespace cynoFinTool {

    static double bs_max_double(double a, double b) {
        return (a > b) ? a : b;
    }

    // Placeholder for standard normal CDF function
    double StdNormCDF(double x) {
        return 0.5 * erfc(-x / sqrt(2.0));
    }
    

    double cyno_Call(double S, double K, double totalSigma) {
        
        if (totalSigma == 0.0 || K == 0.0 || S / K <= 0.0) {
            return bs_max_double(S - K, 0.0);
        }
        double dPlus = (log(S / K) + 0.5 * totalSigma * totalSigma) / totalSigma;
        double dMinus = (log(S / K) - 0.5 * totalSigma * totalSigma) / totalSigma;
        return S * StdNormCDF(dPlus) - K * StdNormCDF(dMinus);
    }

    CYNOMATHUTILITY_API double __stdcall cyno_Put(double S, double K, double totalSigma) {
        return K - S + cyno_Call(S, K, totalSigma);
    }

    CYNOMATHUTILITY_API double __stdcall cyno_DigitalCall(double S, double K, double totalSigma) {
        if (totalSigma == 0.0 || K == 0.0 || S / K <= 0.0) {
            return (S - K > 0.0) ? 1.0 : 0.0;
        }
        double dMinus = (log(S / K) - 0.5 * totalSigma * totalSigma) / totalSigma;
        return StdNormCDF(dMinus);
    }

    CYNOMATHUTILITY_API double __stdcall cyno_DigitalPut(double S, double K, double totalSigma) {
        if (totalSigma == 0.0 || K == 0.0 || S / K <= 0.0) {
            return (S - K < 0.0) ? 1.0 : 0.0;
        }
        double dMinus = (log(S / K) - 0.5 * totalSigma * totalSigma) / totalSigma;
        return StdNormCDF(-dMinus);
    }

    CYNOMATHUTILITY_API double __stdcall cyno_BS(double Tau, double TD, double S, double K, double rD, double rF, double sigma, const char* callPut)
    {
      //  std::transform(callPut.begin(), callPut.end(), callPut.begin(), ::toupper);

        if (Tau <= 0.0 || sigma <= 0.0) return bs_max_double((strcmp(callPut, "CALL") == 0) ? S - K : K - S, 0.0);

        double sqrtTau = sqrt(Tau);
        double d1 = (log(S / K) + (rD - rF) * TD + 0.5 * sigma * sigma * Tau) / (sigma * sqrtTau);
        double d2 = d1 - sigma * sqrtTau;

        if (strcmp(callPut, "CALL") == 0) {
            return S * exp(-rF * TD) * StdNormCDF(d1) - K * exp(-rD * TD) * StdNormCDF(d2);
        }
        else {
            return K * exp(-rD * TD) * StdNormCDF(-d2) - S * exp(-rF * TD) * StdNormCDF(-d1);
        }
    }
    CYNOMATHUTILITY_API double __stdcall cyno_BS_Call(double Tau, double TD, double S, double K, double rD, double rF, double sigma) {
		//char* thisCallPut;
  //      thisCallPut == "CALL";
  //      if (Tau <= 0.0 || sigma <= 0.0) return std::max(S - K, 0.0);
        
		// Call the main Black-Scholes function
       return cyno_BS(Tau, TD, S, K, rD, rF, sigma, "CALL");
    }

    CYNOMATHUTILITY_API double __stdcall cyno_BS_Put(double Tau, double TD, double S, double K, double rD, double rF, double sigma) {

       // char* thisCallPut;
       // thisCallPut = &"PUT";
      //  if (Tau <= 0.0 || sigma <= 0.0) return std::max(K - S, 0.0);
        
		// Call the main Black-Scholes function
       return cyno_BS(Tau, TD, S, K, rD, rF, sigma, "PUT");
    }



//} // namespace BlackScholes
