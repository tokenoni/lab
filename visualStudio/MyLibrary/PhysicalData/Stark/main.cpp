#include <PhysicalData\Stark\Stark_Stehle.hpp>
int main(){
	myPhysicalData::stark_Stehle stark;
	stark.read("C:\\Users\\KeisukeFujii\\Dropbox\\visual_studio\\MyLibrary\\PhysicalData\\Stark\\hydrogen_atom\\", myPhysicalData::stark_Stehle::BalmerSeries, 3);
	return 0;
}