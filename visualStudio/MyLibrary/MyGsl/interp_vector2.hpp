#ifndef __MYGSL_INTERP_VECTOR2_HPP__
#define __MYGSL_INTERP_VECTOR2_HPP__

#include "interp.hpp"

namespace mygsl{
//-------------------------------------------------//
//	multi "interp" stack for an identical x        //
//	faster version than interp_vector              //
//-------------------------------------------------//
	class interp_vector2{
	public:
		//---	constructors	---
		interp_vector2(void);
		interp_vector2(const size_t size1_, const size_t size2_, interp::InterpMethod method_ = interp::Cspline);
		interp_vector2(const interp_vector2& obj);
		interp_vector2& operator = (const interp_vector2& obj);
		~interp_vector2(void);

		//--- setting function	---
//		template < class Tyv_ >
//		void set(const Tyv_& x, const Tyv_& y, InterpMethod method_ = Cspline);
		
		//---	for directly setting the contained data	---
		void setSize(const size_t size1_, const size_t size2_, interp::InterpMethod method_ = interp::Cspline);
		void setInterp(bool isSorted = false);
		void sort();
		void set(const size_t i, const size_t j, const double xj, const double yij){data[i].set(j,xj,yij);}
		void set_x(const size_t j, const double xj){ for(size_t i=0; i<size2(); ++i) data[i].set_x(j, xj);}
		void set_y(const size_t i, const size_t j, const double yij){data[i].set_y(j,yij);};
		void set_y_zero(const size_t i){data[i].set_y_zero();};
		void set_all_y_zero(){ for(size_t i=0; i<size2(); ++i) data[i].set_y_zero();}

		//---	operations	---
		interp& operator[] (size_t i){return data[i];}
		double& get_original_x(size_t j){return data[0].get_original_x(j);}
		double& get_original_y(size_t i,size_t j){return data[i].get_original_y(j);}
		std::vector<double> get_original_x(){return data[0].get_original_x();}
		std::vector<double> get_original_y(size_t i){return data[i].get_original_y();}
		void clear(){FreeMemory();};

		//---	high level functions	---
		double get(size_t i, double x){return data[i].get(x, data[0].getAcc());};							//	get value at x
		double getderiv(size_t i, double x){return data[i].getderiv(x, data[0].getAcc());};				//	get 1st derivative
		double getderiv2(size_t i, double x){return data[i].getderiv2(x, data[0].getAcc());};				//	get 2nd derivative
		double getinteg(size_t i, double a, double b){return data[i].getinteg(a,b, data[0].getAcc());};	//	get integration from [a, b]

		//---	tip functions	---
		size_t size1()const{return num1;}	//	array size of interp*
		size_t size2()const{return num2;}	//	array size of each interp
		bool isallocated()const{return allocated;}
		bool issetInterped()const{return setInterped;}
		double get_xmin()const{if(isallocated()) return data[0].get_xmin();else return NULL;}
		double get_xmax()const{if(isallocated()) return data[0].get_xmax();else return NULL;}

	//	actual data set
	protected:
		interp* data;
		bool allocated;
		bool setInterped;
		size_t num1, num2;
		void FreeMemory();
		void Allocate(size_t size1, size_t size2);
		interp::InterpMethod method;
	};
}

#include "interp_vector2.inl"
#endif