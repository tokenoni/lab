#pragma once
#include<complex>
#include<math.h>

namespace MyFunc{
	static double sign(double x, double y){
		if(y < 0.0) return -1.0*fabs(x);
		else return fabs(x);
	}

	static double gamma(double x){
		double pi= 3.14159265358979324e+00;
		double pv= 6.09750757539068576e+00;
		double pr= 5.63606561897560650e-03;
		double p0= 1.22425977326879918e-01;
		double p1= 8.51370813165034183e-01;
		double p2= 2.25023047535618168e+00;
		double p3= 2.09629553538949977e+00;
		double p4= 5.02197227033920907e-01;
		double p5= 1.12405826571654074e-02;
		double q1= 1.00000000000065532e+00;
		double q2= 1.99999999962010231e+00;
		double q3= 3.00000004672652415e+00;
		double q4= 3.99999663000075089e+00;
		double q5= 5.00035898848319255e+00;

		double w = x;
		double y;
		if(x < 0.0)	{
			w=1.0-x;
			y=((((((p5/(w+q5)+p4)/(w+q4)+p3)/(w+q3)+p2)/(w+q2)+p1)/(w+q1)+p0)/w+pr)*exp((w-0.5)*log(w+pv)-w);
			return y=pi/(y*sin(pi*x));
		}else{
			y=((((((p5/(w+q5)+p4)/(w+q4)+p3)/(w+q3)+p2)/(w+q2)+p1)/(w+q1)+p0)/w+pr)*exp((w-0.5)*log(w+pv)-w);
			return y;
		}
	}

	static double ln_gamma(double x){
		double pi= 3.14159265358979324e+00;
		double pv= 6.09750757539068576e+00;
		double pr= 5.63606561897560650e-03;
		double p0= 1.22425977326879918e-01;
		double p1= 8.51370813165034183e-01;
		double p2= 2.25023047535618168e+00;
		double p3= 2.09629553538949977e+00;
		double p4= 5.02197227033920907e-01;
		double p5= 1.12405826571654074e-02;
		double q1= 1.00000000000065532e+00;
		double q2= 1.99999999962010231e+00;
		double q3= 3.00000004672652415e+00;
		double q4= 3.99999663000075089e+00;
		double q5= 5.00035898848319255e+00;

		double w = x;
		double y;
		if(x < 0.0){
			w = 1.0 - x;
			y=log((((((p5/(w+q5)+p4)/(w+q4)+p3)/(w+q3)+p2)/(w+q2)+p1)/(w+q1)+p0)/w+pr)+(w-0.5)*log(w+pv)-w;
			return log(pi/sin(pi*x))-y;
		}else{
			y=log((((((p5/(w+q5)+p4)/(w+q4)+p3)/(w+q3)+p2)/(w+q2)+p1)/(w+q1)+p0)/w+pr)+(w-0.5)*log(w+pv)-w;
			return y;
		}
	}

	static std::complex<double> gamma(std::complex<double> z){
		double pi= 3.14159265358979324e+00;
		double pv= 7.31790632447016203e+00;
		double pr= 1.66327323256597418e-03;
		double p0= 5.07600435957593046e-02;
		double p1= 5.35012283109043333e-01;
		double p2= 2.40511969085241713e+00;
		double p3= 4.62624880633891119e+00;
		double p4= 3.36415438064135324e+00;
		double p5= 6.71454520397933746e-01;
		double p6= 1.49097237136722798e-02;
		double q1= 9.99999999999975753e-01;
		double q2= 2.00000000000603851e+00;
		double q3= 2.99999999944915534e+00;
		double q4= 4.00000003016801681e+00;
		double q5= 4.99999857982434025e+00;
		double q6= 6.00009857740312429e+00;

		std::complex<double> w = z;
		std::complex<double> y;
		if(z.real() < 0.0){
			w = 1.0 - z;
			y=(((((((p6/(w+q6)+p5)/(w+q5)+p4)/(w+q4)+p3)/(w+q3)+p2)/(w+q2)+p1)/(w+q1)+p0)/w+pr)*exp((w-0.5)*log(w+pv)-w);
			return pi/(y*sin(pi*z));
		}else{
			y=(((((((p6/(w+q6)+p5)/(w+q5)+p4)/(w+q4)+p3)/(w+q3)+p2)/(w+q2)+p1)/(w+q1)+p0)/w+pr)*exp((w-0.5)*log(w+pv)-w);
			return y;
		}
	}

	static std::complex<double> ln_gamma(std::complex<double> z){
		double pi= 3.14159265358979324e+00;
		double pv= 7.31790632447016203e+00;
		double pr= 1.66327323256597418e-03;
		double p0= 5.07600435957593046e-02;
		double p1= 5.35012283109043333e-01;
		double p2= 2.40511969085241713e+00;
		double p3= 4.62624880633891119e+00;
		double p4= 3.36415438064135324e+00;
		double p5= 6.71454520397933746e-01;
		double p6= 1.49097237136722798e-02;
		double q1= 9.99999999999975753e-01;
		double q2= 2.00000000000603851e+00;
		double q3= 2.99999999944915534e+00;
		double q4= 4.00000003016801681e+00;
		double q5= 4.99999857982434025e+00;
		double q6= 6.00009857740312429e+00;

		std::complex<double> w = z;
		std::complex<double> y;

		if(z.real() < 0.0){
			y=log(((((((p6/(w+q6)+p5)/(w+q5)+p4)/(w+q4)+p3)/(w+q3)+p2)/(w+q2)+p1)/(w+q1)+p0)/w+pr)+(w-0.5)*log(w+pv)-w;
			y=log(pi/sin(pi*z))-y;
		}else{
			y=log(((((((p6/(w+q6)+p5)/(w+q5)+p4)/(w+q4)+p3)/(w+q3)+p2)/(w+q2)+p1)/(w+q1)+p0)/w+pr)+(w-0.5)*log(w+pv)-w;
		}
		double r = 0.5-y.imag()/(2.0*pi);
		double s = (int)r-1.0;
		std::complex<double> val(0.0,(2*pi)*((int)(r-s)+s));
		return val;
	}

	static double erfc(double x){
		double pv= 1.26974899965115684e+01;
		double ph= 6.10399733098688199e+00;
		double p0= 2.96316885199227378e-01;
		double p1= 1.81581125134637070e-01;
		double p2= 6.81866451424939493e-02;
		double p3= 1.56907543161966709e-02;
		double p4= 2.21290116681517573e-03;
		double p5= 1.91395813098742864e-04;
		double p6= 9.71013284010551623e-06;
		double p7= 1.66642447174307753e-07;
		double q0= 6.12158644495538758e-02;
		double q1= 5.50942780056002085e-01;
		double q2= 1.53039662058770397e+00;
		double q3= 2.99957952311300634e+00;
		double q4= 4.95867777128246701e+00;
		double q5= 7.41471251099335407e+00;
		double q6= 1.04765104356545238e+01;
		double q7= 1.48455557345597957e+01;
	
		double y = x*x;
		y = exp(-y)*x*(p7/(y+q7)+p6/(y+q6)+p5/(y+q5)+p4/(y+q4)+p3/(y+q3)+p2/(y+q2)+p1/(y+q1)+p0/(y+q0));
		if(x < ph) y=y+2/(exp(pv*x)+1);
		return y;
	}

	static double erf(double x){
		double p0= 1.12837916709551257e+00;
		double p1=-3.76126389031833602e-01;
		double p2= 1.12837916706621301e-01;
		double p3=-2.68661698447642378e-02;
		double p4= 5.22387877685618101e-03;
		double p5=-8.49202435186918470e-04;

		double y = fabs(x);
		if(y > 0.125) return sign(1.0-erfc(y), x);
		else{
			y=x*x;
			return (((((p5*y+p4)*y+p3)*y+p2)*y+p1)*y+p0)*x;
		}
	}
	static std::complex<double> erfc(std::complex<double> x){
		double pv = 1.27813464856668857e+01;
		double ph = 6.64067324283344283e+00;
		double p0 = 2.94608570191793668e-01;
		double p1 = 1.81694307871527086e-01;
		double p2 = 6.91087778921425355e-02;
		double p3 = 1.62114197106901582e-02;
		double p4 = 2.34533471539159422e-03;
		double p5 = 2.09259199579049675e-04;
		double p6 = 1.15149016557480535e-05;
		double p7 = 3.90779571296927748e-07;
		double p8 = 8.17898509247247602e-09;
		double p9 = 1.05575446466983499e-10;
		double p10= 8.40470321453263734e-13;
		double p11= 4.12646136715431977e-15;
		double p12= 1.24947948599560084e-17;
		double q0 = 6.04152433382652546e-02;
		double q1 = 5.43737190044387291e-01;
		double q2 = 1.51038108345663136e+00;
		double q3 = 2.96034692357499747e+00;
		double q4 = 4.89363471039948562e+00;
		double q5 = 7.31024444393009580e+00;
		double q6 = 1.02101761241668280e+01;
		double q7 = 1.35934297511096823e+01;
		double q8 = 1.74600053247586586e+01;
		double q9 = 2.18099028451137569e+01;
		double q10= 2.66431223121749773e+01;
		double q11= 3.19596637259423197e+01;
		double q12= 3.77595270864157841e+01;
		double r0 = 1.56478036351085356e-01;
		double r1 = 2.45771407110492625e-01;
		double r2 = 1.19035163906534275e-01;
		double r3 = 3.55561834455977740e-02;
		double r4 = 6.55014550718381002e-03;
		double r5 = 7.44188068433574137e-04;
		double r6 = 5.21447276257559040e-05;
		double r7 = 2.25337799750608244e-06;
		double r8 = 6.00556181041662576e-08;
		double r9 = 9.87118243564461826e-10;
		double r10= 1.00064645539515792e-11;
		double r11= 6.25587539334288736e-14;
		double r12= 2.41207864479170276e-16;
		double s1 = 2.41660973353061018e-01;
		double s2 = 9.66643893412244073e-01;
		double s3 = 2.17494876017754917e+00;
		double s4 = 3.86657557364897629e+00;
		double s5 = 6.04152433382652546e+00;
		double s6 = 8.69979504071019666e+00;
		double s7 = 1.18413876942999899e+01;
		double s8 = 1.54663022945959052e+01;
		double s9 = 1.95745388415979425e+01;
		double s10= 2.41660973353061018e+01;
		double s11= 2.92409777757203832e+01;
		double s12= 3.47991801628407866e+01;

		std::complex<double> y=x*x;

		if((fabs(x.real()) + fabs(x.imag())) < ph){
			std::complex<double> z=exp(pv*x);
			if(z.real() >= 0.0)
				y=exp(-y)*x*(p12/(y+q12)+p11/(y+q11)+p10/(y+q10)+p9/(y+q9)+p8/(y+q8)+p7/(y+q7)+p6/(y+q6)+p5/(y+q5)+p4/(y+q4)+p3/(y+q3)+p2/(y+q2)+p1/(y+q1)+p0/(y+q0))+2.0/(1.0+z);
			else
				y=exp(-y)*x*(r12/(y+s12)+r11/(y+s11)+r10/(y+s10)+r9/(y+s9)+r8/(y+s8)+r7/(y+s7)+r6/(y+s6)+r5/(y+s5)+r4/(y+s4)+r3/(y+s3)+r2/(y+s2)+r1/(y+s1)+r0/y)+2.0/(1.0-z);
		}else{
			y=exp(-y)*x*(p12/(y+q12)+p11/(y+q11)+p10/(y+q10)+p9/(y+q9)+p8/(y+q8)+p7/(y+q7)+p6/(y+q6)+p5/(y+q5)+p4/(y+q4)+p3/(y+q3)+p2/(y+q2)+p1/(y+q1)+p0/(y+q0));
			if(x.real() < 0.0) y=y+2.0;
		}
		return y;
	}
	
	static std::complex<double> erf(std::complex<double> x){
		double p0= 1.12837916709551257e+00;
		double p1=-3.76126389031837525e-01;
		double p2= 1.12837916709551257e-01;
		double p3=-2.68661706451312518e-02;
		double p4= 5.22397762544218784e-03;
		double p5=-8.54832702345085283e-04;
		double p6= 1.20553329817896643e-04;

		std::complex<double> y;
		if((fabs(x.real())+fabs(x.imag())) > 0.125){
			if(x.real() >= 0.0)	y = 1.0 - erfc( x);
			else                y = erfc(-x) - 1.0;
			return y;
		}else{
			y=x*x;
			return ((((((p6*y+p5)*y+p4)*y+p3)*y+p2)*y+p1)*y+p0)*x;
		}
	}
	

};
