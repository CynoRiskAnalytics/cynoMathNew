// single_barrier_option.cpp
#include "pch.h"
#include "framework.h"
#include <math.h>
#include <string.h>

//namespace cynoFinTool {

    static double single_max_double(double a, double b) {
        return (a > b) ? a : b;
    }

    double singleStdNormCDF(double x) {
        return 0.5 * erfc(-x / sqrt(2.0));
    }


  //  std::tuple<double, double, double, double> CalcIntermediateValues(double S, double K, double B, double tau, double TD, double rd, double rf, double sigma, const std::string& knockType, const std::string& callput) {
        
        void CalcIntermediateValues(
            double S, double K, double B,
            double tau, double TD,
            double rd, double rf, double sigma,
            const char* knockType,
            const char* callput,
            double* A1, double* A2,
            double* A3, double* A4){ 
        
        
        const char* callPutUpper = callput;
       // std::transform(callPutUpper.begin(), callPutUpper.end(), callPutUpper.begin(), ::toupper);
        const char* knockTypeUpper = knockType;
       // std::transform(knockTypeUpper.begin(), knockTypeUpper.end(), knockTypeUpper.begin(), ::toupper);

        double phi;
        if (strcmp(callPutUpper, "CALL") == 0) { 
            phi = 1.0; }
        else {
            phi = -1.0; }
        double eta;
        if (strcmp(knockTypeUpper, "UPANDIN") == 0 || strcmp(knockTypeUpper, "UPANDOUT") == 0) {
            eta = -1.0; }
        else {
            eta = 1.0; } 

        double df_d = exp(-rd * TD);
        double df_f = exp(-rf * TD);

        double df_d_end = exp(-rd * tau);
        double df_f_end = exp(-rf * tau);

        double x = single_max_double(S, 1e-12);
        double k = single_max_double(K, 1e-12);
        double b = single_max_double(B, 1e-12);

        double theta1 = (rd - rf) / sigma + sigma / 2;
        double theta2 = (rd - rf) / sigma - sigma / 2;

        double sqrt_tau = sqrt(tau);

        double d1 = (log(x / k) + sigma * theta1 * tau) / (sigma * sqrt_tau);
        double d2 = (log(x / k) + sigma * theta2 * tau) / (sigma * sqrt_tau);

        double x1 = (log(x / b) + sigma * theta1 * tau) / (sigma * sqrt_tau);
        double x2 = (log(x / b) + sigma * theta2 * tau) / (sigma * sqrt_tau);

        double z1 = (log(b * b / x / k) + sigma * theta1 * tau) / (sigma * sqrt_tau);
        double z2 = (log(b * b / x / k) + sigma * theta2 * tau) / (sigma * sqrt_tau);

        double y1 = (log(b / x) + sigma * theta1 * tau) / (sigma * sqrt_tau);
        double y2 = (log(b / x) + sigma * theta2 * tau) / (sigma * sqrt_tau);

        *A1 = phi * x * df_f * singleStdNormCDF(phi * d1) - phi * k * df_d * singleStdNormCDF(phi * d2);
        *A2 = phi * x * df_f * singleStdNormCDF(phi * x1) - phi * k * df_d * singleStdNormCDF(phi * x2);

        double pow_term = pow(b / x, 2 * theta2 / sigma);
        *A3 = phi * pow_term * (x * df_f * pow(b / x, 2.0) * singleStdNormCDF(eta * z1) - k * df_d * singleStdNormCDF(eta * z2));
        *A4 = phi * pow_term * (x * df_f * pow(b / x, 2.0) * singleStdNormCDF(eta * y1) - k * df_d * singleStdNormCDF(eta * y2));

        //return { A1, A2, A3, A4 };
    }

    double cyno_singleBarrier(double S, double K, double B, double Tau, double TD, double rd, double rf, double sigma, const char* knockType, const char* callput) {
        
        //
        
        const char* knock = knockType;
       // std::transform(knock.begin(), knock.end(), knock.begin(), ::toupper);
        const char* cp = callput;
        //std::transform(cp.begin(), cp.end(), cp.begin(), ::toupper);

        //if (S >= B && knock.find("UP") != std::string::npos)
        //    throw std::invalid_argument("Barrier must be above spot for up barriers.");
        //if (S <= B && knock.find("DOWN") != std::string::npos)
        //    throw std::invalid_argument("Barrier must be below spot for down barriers.");

       // auto [A1, A2, A3, A4] = CalcIntermediateValues(S, K, B, Tau, TD, rd, rf, sigma, knock, cp);
		double A1, A2, A3, A4;

        CalcIntermediateValues(S, K, B, Tau, TD, rd, rf, sigma, knock, cp, &A1, &A2, &A3, &A4);
        
        double res = 0.0;

        if (strcmp(knock, "UPANDIN")==0) {
            if (strcmp(cp, "CALL")==0) res = (K > B) ? A1 : (A2 - A3 + A4);
            else res = (K > B) ? (A1 - A2 + A4) : A3;
        }
        else if (strcmp(knock, "DOWNANDIN") == 0) {
            if (strcmp(cp, "CALL") == 0) res = (K > B) ? A3 : (A1 - A2 + A4);
            else res = (K > B) ? (A2 - A3 + A4) : A1;
        }
        else if (strcmp(knock, "UPANDOUT") == 0) {
            if (strcmp(cp, "CALL") == 0) res = (K > B) ? 0.0 : (A1 - A2 + A3 - A4);
            else res = (K > B) ? (A2 - A4) : (A1 - A3);
        }
        else if (strcmp(knock, "DOWNANDOUT") == 0) {
            if (strcmp(cp, "CALL") == 0) res = (K > B) ? (A1 - A3) : (A2 - A4);
            else res = (K > B) ? (A1 - A2 + A3 - A4) : 0.0;
        }
        else {
          //  throw std::invalid_argument("Unsupported knock type: " + knock);
        }

        return res;
    }

    double cyno_oneTouch(double Tau, double TD, double S, double B, double rd, double rf, double sigma, const char* upDown, const char* nCurrency, const char* nPaidAt) {
        
        const char* up = upDown;
        const char* pay = nPaidAt;
        const char* curr = nCurrency;
        //std::transform(up.begin(), up.end(), up.begin(), ::toupper);
        //std::transform(pay.begin(), pay.end(), pay.begin(), ::toupper);
        //std::transform(curr.begin(), curr.end(), curr.begin(), ::toupper);

        //if (S >= B && up == "UP") throw std::invalid_argument("Barrier must be above spot for up barrier");
        //if (S <= B && up == "DOWN") throw std::invalid_argument("Barrier must be below spot for down barrier");

        double epsp;
        if (strcmp(up, "UP") ==0) { epsp = -1; } else { epsp = 1.0; }
        double mu = rd - rf;
        double muh = mu * TD / Tau - 0.5 * sigma * sigma;
        double mup = sqrt(muh * muh + 2.0 * rd * sigma * sigma);

        double z, res = 0.0;

        if (strcmp(pay, "ATHIT") == 0) {
            z = log(B / S) / (sigma * sqrt(Tau)) + mup / sigma * sqrt(Tau);
            res = singleStdNormCDF(epsp * z) * pow(B / S, (muh + mup) / (sigma * sigma)) +
                singleStdNormCDF(epsp * z - 2.0 * epsp * mup / sigma * sqrt(Tau)) * pow(B / S, (muh - mup) / (sigma * sigma));
        }
        else {
            z = log(B / S) / (sigma * sqrt(Tau)) + muh / sigma * sqrt(Tau);
            res = exp(-rd * TD) * (singleStdNormCDF(epsp * z) * pow(B / S, 2.0 * muh / (sigma * sigma)) +
                singleStdNormCDF(epsp * z - 2.0 * epsp * muh / sigma * sqrt(Tau)));
        }

        if (strcmp(curr, "FOR") == 0) {
            if (strcmp(pay, "ATMATURITY") == 0) {
                if (strcmp(up, "DOWN") == 0) {
                    res = B * cyno_oneTouch(Tau, TD, S, B, rd, rf, sigma, up, "DOM", pay) +
                        cyno_singleBarrier(S, B, B, Tau, TD, rd, rf, sigma, "DOWNANDIN", "CALL") -
                        cyno_singleBarrier(S, B, B, Tau, TD, rd, rf, sigma, "DOWNANDIN", "PUT");
                }
                else {
                    res = B * cyno_oneTouch(Tau, TD, S, B, rd, rf, sigma, up, "DOM", pay) +
                        cyno_singleBarrier(S, B, B, Tau, TD, rd, rf, sigma, "UPANDIN", "CALL") -
                        cyno_singleBarrier(S, B, B, Tau, TD, rd, rf, sigma, "UPANDIN", "PUT");
                }
            }
            else {
                res = B * res;
            }
        }

        return res;
    }

    double cyno_hitProbabilityOneBarrier(double Tau, double TD, double S, double H, double rD, double rF, double sigma, const char* nCurrency) {
        double res = 0.0;
        if (S < H) {
            res = exp(rD * TD) * cyno_oneTouch(Tau, TD, S, H, rD, rF, sigma, "UP", nCurrency, "ATMATURITY");
            if (strcmp(nCurrency, "FOR") == 0) res /= S;
        }
        else {
            res = exp(rD * TD) * cyno_oneTouch(Tau, TD, S, H, rD, rF, sigma, "DOWN", nCurrency, "ATMATURITY");
            if (strcmp(nCurrency, "FOR") == 0) res /= S;
        }
        return res;
    }

    double cyno_oneNoTouch(double Tau, double TD, double S, double B, double rd, double rf, double sigma, const char* upDown, const char* nCurrency, const char* nPaidAt) {
      
        
        if (strcmp(nPaidAt, "ATHIT") == 0) return 0.0;
        double res = 0.0;

        if (strcmp(nCurrency, "FOR") == 0) {
            res = exp(-rf * TD) * S - cyno_oneTouch(Tau, TD, S, B, rd, rf, sigma, upDown, nCurrency, nPaidAt);
        }
        else {
            res = exp(-rd * TD) - cyno_oneTouch(Tau, TD, S, B, rd, rf, sigma, upDown, nCurrency, nPaidAt);
            //res =  cyno_oneTouch(Tau, TD, S, B, rd, rf, sigma, upDown, nCurrency, nPaidAt);
        }

        return res;
    }

//} // namespace cynoFinTool
