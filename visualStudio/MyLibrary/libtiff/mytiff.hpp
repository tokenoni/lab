#ifndef __MYTIFF_HPP__
#define __MYTIFF_HPP__

#include <stdlib.h>
#include <libtiff\tiffio.h>
#include <string>
#include <new>

class mytiff{
public:
	mytiff(void){ isAllocated = false;}
	~mytiff(void);
	void read(const std::string& filename_);
	unsigned short get16bitData(size_t i, size_t j) const;
	unsigned char get8bitData(size_t i, size_t j) const;
	unsigned int get32bitData(size_t i, size_t j) const;
	size_t getBps()const{ return (size_t)bps; }
	size_t size0()const{ return size_i; };
	size_t size1()const{ return size_j; };

private:
	std::string filename;
	size_t bps;
	char* data;
	bool isAllocated;
	size_t size_i, size_j;
};


#include "mytiff.inl"
#endif