namespace mygsl{
	//	
	//	member functions of SphericalHarmonics expansion
	// 
	inline SphericalHarmonicsExpansion::SphericalHarmonicsExpansion(const size_t rank_)
			: rank_max(4), rank((rank_>=rank_max)?rank_max:rank_),
		  	  coef_l0(    1.0/ 2.0 * sqrt(  1.0/( 1.0*myconst::pi))),
			  //	l=1	 	 			 
			  coef_l1_m0( 1.0/ 2.0 * sqrt(  3.0/( 1.0*myconst::pi))),
			  coef_l1_m1(-1.0/ 2.0 * sqrt(  3.0/( 2.0*myconst::pi))),
			  //	l=2	 	 			 
			  coef_l2_m0( 1.0/ 4.0 * sqrt(  5.0/( 1.0*myconst::pi))),
			  coef_l2_m1(-1.0/ 2.0 * sqrt( 15.0/( 3.0*myconst::pi))),
			  coef_l2_m2( 1.0/ 2.0 * sqrt( 15.0/( 2.0*myconst::pi))),
			  //	l=3	 	 			 
			  coef_l3_m0( 1.0/ 4.0 * sqrt(  7.0/( 1.0*myconst::pi))),
			  coef_l3_m1(-1.0/ 8.0 * sqrt( 21.0/( 1.0*myconst::pi))),
			  coef_l3_m2( 1.0/ 4.0 * sqrt(105.0/( 2.0*myconst::pi))),
			  coef_l3_m3(-1.0/ 8.0 * sqrt( 35.0/( 1.0*myconst::pi))),
			  //	l=4	 				 
			  coef_l4_m0( 3.0/16.0 * sqrt(  1.0/( 1.0*myconst::pi))),
			  coef_l4_m1(-3.0/ 8.0 * sqrt(  5.0/( 1.0*myconst::pi))),
			  coef_l4_m2( 3.0/ 8.0 * sqrt(  5.0/( 2.0*myconst::pi))),
			  coef_l4_m3(-3.0/ 8.0 * sqrt( 35.0/( 1.0*myconst::pi))),
			  coef_l4_m4( 3.0/16.0 * sqrt( 35.0/( 2.0*myconst::pi))),
			  size(rankToSize(rank_))
	{	
		sh_values = new complex [size];	
		sh_values[0] = complex(coef_l0, 0.0);
	};

	inline SphericalHarmonicsExpansion::~SphericalHarmonicsExpansion(void)
	{		delete[] sh_values;	}

	inline void SphericalHarmonicsExpansion::run(const double theta, const double phi){
		if(rank == 0) return;
		double cos_theta = cos(theta);
		double sin_theta = sin(theta);
		//	j=1
//		complex exp_iphi = exp(complex(0.0, phi));
		complex exp_iphi = exp_euler(phi);
		//	m=0
		sh_values[1] = complex(coef_l1_m0 * cos_theta, 0.0);
		//	m=1, -1
		sh_values[2] = exp_iphi * coef_l1_m1 * sin_theta;
		sh_values[3] = -sh_values[2].conj();
		if(rank==1)	return;
		
		//	j=2
//		complex exp_2iphi = exp(complex(0.0, 2.0*phi));
		complex exp_2iphi = exp_euler(2.0*phi);
		double cos_theta2 = cos_theta*cos_theta;
		double sin_theta2 = sin_theta*sin_theta;
		//	m=0
		sh_values[4] = complex(coef_l2_m0 * (3.0*cos_theta2-1.0), 0.0);
		//	m=1, -1
		sh_values[5] = exp_iphi  * coef_l2_m1 * sin_theta*cos_theta;
		sh_values[6] = -sh_values[5].conj();
		//	m=2, -2
		sh_values[7] = exp_2iphi * coef_l2_m2 * sin_theta2;
		sh_values[8] =  sh_values[7].conj();
		if(rank==2)	return;

		// j=3
//		complex exp_3iphi = exp(complex(0.0, 3.0*phi));
		complex exp_3iphi = exp_euler(3.0*phi);
		double cos_theta3 = cos_theta2*cos_theta;
		double sin_theta3 = sin_theta2*sin_theta;
		//	m=0
		sh_values[9] = complex(coef_l3_m0 * (5.0*cos_theta3- 3.0*cos_theta), 0.0);
		//	m=1, -1
		sh_values[10]= exp_iphi * coef_l3_m1 * sin_theta * (5.0*cos_theta2 - 1.0);
		sh_values[11]= -sh_values[10].conj();
		//	m=2, -2
		sh_values[12]= exp_2iphi * coef_l3_m2 * sin_theta2 * cos_theta;
		sh_values[13]=  sh_values[12].conj();
		//	m=3, -3
		sh_values[14]= exp_3iphi * coef_l3_m3 * sin_theta3;
		sh_values[15]= -sh_values[14].conj();
		if(rank==3)	return;

		//	j=4
//		complex exp_4iphi = exp(complex(0.0, 4.0*phi));
		complex exp_4iphi = exp_euler(4.0*phi);
		double cos_theta4 = cos_theta2*cos_theta2;
		double sin_theta4 = sin_theta2*sin_theta2;
		//	m=0
		sh_values[16]= complex(coef_l4_m0 * (35.0*cos_theta4 - 30.0*cos_theta2 + 3.0), 0.0);
		//	m=1
		sh_values[17]= exp_iphi * coef_l4_m1 * sin_theta * (7.0*cos_theta3 - 3.0*cos_theta);
		sh_values[18]= -sh_values[17].conj();
		//	m=2
		sh_values[19]= exp_2iphi* coef_l4_m2 * sin_theta2 * (7.0*cos_theta2 - 1.0);
		sh_values[20]=  sh_values[19].conj();
		//	m=3
		sh_values[21]= exp_3iphi* coef_l4_m3 * sin_theta3 * cos_theta;
		sh_values[22]= -sh_values[21].conj();
		//	m=4
		sh_values[23]= exp_4iphi* coef_l4_m4 * sin_theta4;
		sh_values[24]=  sh_values[23].conj();
		
		return;
	}

	inline void SphericalHarmonicsExpansion::run(const double x, const double y, const double z){
		if(rank == 0) return;
		//	j=1
		//	m=0
		sh_values[1] = coef_l1_m0*z;
		//	m=1, -1
		sh_values[2] = complex(x, y)*coef_l1_m1;
		sh_values[3] = -sh_values[2].conj();
		if(rank==1)	return;
		
		//	j=2
		//	m=0
		double xx = x*x, xy = x*y, xz = x*z;
		double yy = y*y, yz = y*z;
		double zz = z*z;
		complex x_plus_iy_2nd = complex(xx-yy, 2*xy);
		sh_values[4] = 2.0*zz-xx-yy;
		//	m=1, -1
		sh_values[5] = complex(xz, yz)*coef_l2_m1;
		sh_values[6] = -sh_values[5].conj();
		//	m=2, -2
		sh_values[7] = x_plus_iy_2nd*coef_l2_m2;
		sh_values[8] =  sh_values[7].conj();
		if(rank==2)	return;

		// j=3
		complex x_plus_iy_3rd = complex(xx*(x-3.0*y), (3.0*x-y)*yy);
		//	m=0
		sh_values[9] = z*(2.0*zz-3.0*xx-3.0*yy)*coef_l3_m0;
		//	m=1, -1
		sh_values[10]= complex(x,y)*(coef_l3_m1*(zz-xx-yy));
		sh_values[11]= -sh_values[10].conj();
		//	m=2, -2
		sh_values[12]= x_plus_iy_2nd *(coef_l3_m2*z);
		sh_values[13]=  sh_values[12].conj();
		//	m=3, -3
		sh_values[14]= x_plus_iy_3rd *coef_l3_m3;
		sh_values[15]= -sh_values[14].conj();
		if(rank==3)	return;

		//	j=4
		double rr = xx+yy+zz;
		complex x_plus_iy_4th = x_plus_iy_2nd*x_plus_iy_2nd;
		//	m=0
		sh_values[16]= (35.0*zz*zz-30.0*zz*rr+3.0*rr*rr)*coef_l4_m0;
		//	m=1
		sh_values[17]= complex(x,y)*(z*(7.0*zz-3.0*rr)*coef_l4_m1);
		sh_values[18]= -sh_values[17].conj();
		//	m=2
		sh_values[19]= x_plus_iy_2nd*((7.0*zz-rr)*coef_l4_m2);
		sh_values[20]=  sh_values[19].conj();
		//	m=3
		sh_values[21]= x_plus_iy_3rd*(z*coef_l4_m3);
		sh_values[22]= -sh_values[21].conj();
		//	m=4
		sh_values[23]= x_plus_iy_4th*coef_l4_m4;
		sh_values[24]=  sh_values[23].conj();
		
		return;
	}
};