#ifndef __IGORIBW_HPP__
#define __IGORIBW_HPP__

#include "base_functions\IGORibw_base.hpp"
#include <string>

namespace IGORdata{
	template < class Ty_ >
	static bool write_ibw( const Ty_ & data, const size_t size, const std::string& file_path, const std::string& wavename, const std::string& comment = "");
	
	template < class Ty_ >
	static bool write_ibw_onedim_double( const Ty_ & data, const size_t size, const std::string& file_path, const std::string& wavename, const std::string& comment = "");

};

#include "IGORibw.inl"

#endif
