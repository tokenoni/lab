#include <AndorSIF\AndorSIF.h>
#include <IGOR\IGORitx.hpp>
void test_azuma();
void test_azuma2(const std::string, const std::string, const std::string);

int main(){
	//AndorSIF andorsif;
	test_azuma();

}
void test_azuma(){
	std::string indir = "\\\\10.232.53.161\\disk1\\thesis\\2011_Feb_Shuron_Azuma\\final_data\\THR1000_finaldata\\20100927_27kpa_3mA\\rawdata\\";
	std::string outdir= "";
	
	test_azuma2(indir, "3889.sif", outdir);	
}

void test_azuma2(const std::string indir, const std::string filename, const std::string outdir){
	AndorSIF andorsif;
//	andorsif.

}