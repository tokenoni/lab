#include <libtiff\mytiff.hpp>
int main(){
	mytiff tiffdata;
	tiffdata.read("C:\\Users\\KeisukeFujii\\Dropbox\\future_works\\interferometer_with_multipleCamera\\20140509_pointGrayCameraTest\\He_lamp.tif");
	getchar();
	return 0;
}