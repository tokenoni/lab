namespace myfunc{
		// still under construction
	inline double fast_exp(const double x){
		double a = LOG2E * x;
		a -= (a < 0);
		int n = (int)a;
		double b = x- n*LOG2E;

	//	 exponential value of b when b is between 0 and ln 2
		double y;
		y = 1.185268231308989403584147407056378360798378534739e-2;
		y *= b;
		y += 3.87412011356070379615759057344100690905653320886699e-2;
		y *= b;
		y += 0.16775408658617866431779970932853611481292418818223;
		y *= b;
		y += 0.49981934577169208735732248650232562589934399402426;
		y *= b;
		y += 1.00001092396453942157124178508842412412025643386873;
		y *= b;
		y += 0.99999989311082729779536722205742989232069120354073;

	//	returns the value of 2^n
		typedef union {
			double d;
			unsigned short s[4];
		} ieee754;

		ieee754 u;
		u.d = 0.0;
		u.s[3] = (unsigned short)(((n + 1023) << 4) & 0x7FF0);
	//	return y*u.d;
	return 0.0;
	}

};