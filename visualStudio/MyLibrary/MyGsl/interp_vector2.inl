//------------------	definitions	------------------

namespace mygsl{
//-------------------------------------------------//
//												   //
//---		1 dimensional interpolation			---//
//												   //
//-------------------------------------------------//

	inline interp_vector2::interp_vector2(void){allocated = false;}
	inline interp_vector2::interp_vector2(const size_t size1_, const size_t size2_, interp::InterpMethod method_){
		method = method_;
		allocated = false;
		Allocate(size1_, size2_);
	}

	inline interp_vector2::interp_vector2(const interp_vector2&obj){
		allocated = false;
		FreeMemory();
		if(obj.isallocated()){
			method = obj.method;
			Allocate(obj.size1(), obj.size2());
			for(size_t i=0; i<obj.size1(); ++i) data[i] = obj.data[i];
			if(obj.issetInterped())	setInterp();
		}
	}
	inline interp_vector2& interp_vector2::operator=(const interp_vector2& obj){
		FreeMemory();
		method = obj.method;
		if(obj.isallocated()){
			Allocate(obj.size1(), obj.size2());
			for(size_t i=0; i<obj.size1(); ++i) data[i] = obj.data[i];
			if(obj.issetInterped())	setInterp();
		}
		return *this;
	}
	inline interp_vector2::~interp_vector2(void){FreeMemory();	}

	inline void interp_vector2::FreeMemory(){
		if(allocated){
			for(size_t i=0; i<num1; ++i) data[i].clear();
			delete [] data;
		}
		allocated = false;
		setInterped = false;
		num1 = 0; num2=0;
	}

	inline void interp_vector2::setSize(const size_t size1_, const size_t size2_, interp::InterpMethod method_){
		if(isallocated() && ((size1() != size1_ || size2() != size2_) || method_ != method)){
			method = method_;
			FreeMemory();
			Allocate(size1_, size2_);
		}
		if(!isallocated()){
			method = method_;
			Allocate(size1_, size2_);
		}
	}
	
	inline void interp_vector2::setInterp(bool isSorted_){
/*		if(interp::method == CsplineFlatEdge)
			CsplineClampedInit(interp_obj, x_array, data, num,0.0, 0.0);
		else
		*/
		for(size_t i=0; i<size1(); ++i) data[i].setInterp(isSorted_);
		setInterped = true;
	}
	
	inline void interp_vector2::Allocate(size_t size1_, size_t size2_){
		num1 = size1_;
		num2 = size2_;
		if(!allocated){
			data = new interp [num1];
			for(size_t i=0; i<size1(); ++i) data[i].setSize(size2_, method);
		}
		allocated = true;
		setInterped = false;
	}


}