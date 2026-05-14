
#pragma once
//namespace cynoFinTool {
//	// Computes the price of a European call option (no discounting, totalSigma = sigma * sqrt(T))
//	double cyno_Call(double S, double K, double totalSigma);
//	// Computes the price of a European put option via put-call parity (no discounting)
//	double cyno_Put(double S, double K, double totalSigma);
//	// Digital (cash-or-nothing) call option
//	double cyno_DigitalCall(double S, double K, double totalSigma);
//	// Digital (cash-or-nothing) put option
//	double cyno_DigitalPut(double S, double K, double totalSigma);
//} // namespace cynoFinTool
//
//#ifdef CYNOMATHUTILITY_EXPORTS
//#define CYNOMATHUTILITY_API __declspec(dllexport)
//#else
//#define CYNOMATHUTILITY_API __declspec(dllimport)
//#endif
//
//// Computes the price of a European call option (no discounting, totalSigma = sigma * sqrt(T))
//CYNOMATHUTILITY_API double __stdcall cyno_Call(double S, double K, double totalSigma);
//
//// Computes the price of a European put option via put-call parity (no discounting)
//CYNOMATHUTILITY_API double __stdcall cyno_Put(double S, double K, double totalSigma);
//
//// Digital (cash-or-nothing) call option
//CYNOMATHUTILITY_API double __stdcall cyno_DigitalCall(double S, double K, double totalSigma);
//
//// Digital (cash-or-nothing) put option
//CYNOMATHUTILITY_API double __stdcall cyno_DigitalPut(double S, double K, double totalSigma);
//
//// Full Black-Scholes pricing formula with discounting
//// Tau: time to maturity, TD: discounting time
//// rD: domestic rate, rF: foreign rate
//// sigma: annualized volatility
//// callPut: "CALL" or "PUT"
//CYNOMATHUTILITY_API double __stdcall cyno_BS(double Tau, double TD, double S, double K, double rD, double rF, double sigma, const char* callPut);
//
//#endif // BLACK_SCHOLES_H
