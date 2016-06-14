namespace myfunc{
	static const double LOG2E = log(2.0);

	// still under construction
	//	fast function to calculate expnential value of x
	static double fast_exp(const double x);
	//	returns the exponential value of b when b is between 0 and ln 2
	static double PolynominalApproximationOfExpBetween0andLn2(const double b);
	//	returns the value of 2^n
	static double nthPowerOf2(const int n);
};

#include "fast_functions.inl"