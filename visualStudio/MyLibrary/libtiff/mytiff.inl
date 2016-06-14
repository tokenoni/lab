inline mytiff::~mytiff(void){
	if (isAllocated) delete[] data;
}

inline unsigned short mytiff::get16bitData(size_t i, size_t j) const{
	char tmp1 = data[2 * (i + size_i * j)];
	char tmp2 = data[2 * (i + size_i * j) + 1];

	unsigned short rslt = ~((tmp2 << 8) | (tmp1 & 0x00FF));
	return rslt;
}
inline unsigned char mytiff::get8bitData(size_t i, size_t j) const{
	return ~data[i + size_i * j];
}
inline unsigned int mytiff::get32bitData(size_t i, size_t j) const{
	return 0;
}

inline void mytiff::read(const std::string& filename_){
	filename = filename_;
	TIFF *image;
	uint16 photo, spp, fillorder;
	uint32 width;
	tsize_t stripSize;

	unsigned long imageOffset, result;
	int stripMax;
	char tempbyte;
	unsigned long bufferSize, count;


	// Open the TIFF image 
	if ((image = TIFFOpen(filename.c_str(), "r")) == NULL){
		std::string errorMessage("Could not open incoming image\n");
		throw errorMessage + "in reading" + filename + "\n";
	}
	// Check that it is of a type that we support 
	if ((TIFFGetField(image, TIFFTAG_BITSPERSAMPLE, &bps) == 0) ){
		std::string errorMessage("Either undefined or unsupported number of bits per sample\n");
		throw errorMessage + "in reading" + filename + "\n";
	}
	if ((TIFFGetField(image, TIFFTAG_SAMPLESPERPIXEL, &spp) == 0) || (spp != 1)){
		std::string errorMessage("Either undefined or unsupported number of samples per pixel\n");
		throw errorMessage + "in reading" + filename + "\n";
	}
	// Read in the possibly multiple strips 
	stripSize = TIFFStripSize(image);
	stripMax = TIFFNumberOfStrips(image);
	imageOffset = 0;
	bufferSize = TIFFNumberOfStrips(image) * stripSize * bps;
	try{
		data = new char[bufferSize];
		isAllocated = true;
	}
	catch (std::bad_alloc){
		isAllocated = false;
		std::string errorMessage("allocation failed\n");
		throw errorMessage + "in reading" + filename + "\n";
	}

	for (int stripCount = 0; stripCount < stripMax; stripCount++){
		if ((result = TIFFReadEncodedStrip(image, stripCount,
			data + imageOffset,
			stripSize)) == -1){
			std::string errorMessage("Read error on input strip number\n");
			throw errorMessage + "in reading" + filename + "\n";
		}
		imageOffset += result;
	}
	// Deal with photometric interpretations 
	if (TIFFGetField(image, TIFFTAG_PHOTOMETRIC, &photo) == 0){
		std::string errorMessage("Image has an undefined photometric interpretation\n");
		throw errorMessage + "in reading" + filename + "\n";
	}
	if (photo != PHOTOMETRIC_MINISWHITE){
		// Flip bits 
		//		printf("Fixing the photometric interpretation\n");
		for (count = 0; count < bufferSize; count++)
			data[count] = ~data[count];
	}
	// Deal with fillorder 
	if (TIFFGetField(image, TIFFTAG_FILLORDER, &fillorder) == 0){
		std::string errorMessage("Image has an undefined fillorder\n");
		throw errorMessage + "in reading" + filename + "\n";
	}
	if (fillorder != FILLORDER_MSB2LSB){
		// We need to swap bits -- ABCDEFGH becomes HGFEDCBA 
		//	printf("Fixing the fillorder\n");
		for (count = 0; count < bufferSize; count++){
			tempbyte = 0;
			if (data[count] & 128) tempbyte += 1;
			if (data[count] & 64) tempbyte += 2;
			if (data[count] & 32) tempbyte += 4;
			if (data[count] & 16) tempbyte += 8;
			if (data[count] & 8) tempbyte += 16;
			if (data[count] & 4) tempbyte += 32;
			if (data[count] & 2) tempbyte += 64;
			if (data[count] & 1) tempbyte += 128;
			data[count] = tempbyte;
		}
	}
	// Do whatever it is we do with the buffer -- we dump it in hex 
	if (TIFFGetField(image, TIFFTAG_IMAGEWIDTH, &width) == 0){
		std::string errorMessage("Image does not define its width\n");
		throw errorMessage + "in reading" + filename + "\n";
	}

	//	size
	// Find the width and height of the image 
	TIFFGetField(image, TIFFTAG_IMAGEWIDTH, &size_i);
	TIFFGetField(image, TIFFTAG_IMAGELENGTH, &size_j);
	TIFFClose(image);
};