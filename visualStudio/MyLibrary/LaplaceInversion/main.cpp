#include "LaplaceInversion.hpp"
#include <IGOR/IGORitx.hpp>

int main(){
	std::vector<double> s(11);
	std::vector<double> Fs(11);
	double b = 0.1;
	double d = 0.3;
	for(size_t i=0; i<11; ++i){
		s[i] = (b+i+1)*d;
		Fs[i] = 1.0/((s[i]+1.0)*(s[i]+1.0)+1.0);
	}
	LaplaceInversion test;
	test.set(Fs, b, d, 11);

	std::vector<double> t(128);
	std::vector<double> ft(128);
	std::vector<double> ft_real(128);
	for(size_t i=0; i<128; ++i){
		t[i] = 5.5/127*i;
		ft[i]=test.get(t[i]);
		ft_real[i] = exp(-t[i])*sin(t[i]);
	}
	IGORdata::write_itx(ft, ft_real, "example.itx", "example","example_real");
	IGORdata::write_itx(t, "tt.itx", "tt");
};