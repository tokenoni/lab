namespace IGORdata{
	template < class Ty_ >
	static bool write_ibw_onedim_double( const Ty_ & data, const size_t size, const std::string& file_path, const std::string& wavename, const std::string& comment)
	{
		OriginalFunctions::WaveHeader2 wh;
		unsigned long now;
		CP_FILE_REF fr;
		int err;
		
		#ifdef _WIN32
			now = 0;			// It would be possible to write a Windows equivalent for the Macintosh GetDateTime function but it is not easy.
		#else
			GetDateTime(&now);
		#endif
		
		memset(&wh, 0, sizeof(wh));
		wh.type = OriginalFunctions::NT_FP64;							// Double precision floating point.
		strcpy(wh.bname, wavename.c_str());
		wh.hsA = 1.0;								// 1 per point.
		wh.hsB = 0.0e0;								// Starting from zero.
		wh.npnts = size;
		wh.modDate = now;
		
		#ifdef _WIN32
			wh.platform = 2;						// Windows
		#else
			wh.platform = 1;						// Macintosh
		#endif
		
		double *data_p = new double[size];
		// Fill in data.
		for (size_t p = 0; p < size; p++)
			data_p[p] = data[p];
			
		// Copy the first 16 bytes of data to the wData field of the WaveHeader2 structure.
		// This is necessary for the checksum done by WriteVersion2NumericWave.
		memcpy(wh.wData, data_p, 16);
		
		do {
			// Create a file.
			if (err = OriginalFunctions::CPCreateFile(file_path.c_str(), 1, 'IGR0', 'IGBW')) {
				printf("Error %d creating the file.\n", err);
				break;
			}
				
			// Open the file.
			if (err = OriginalFunctions::CPOpenFile(file_path.c_str(), 1, &fr)) {
				printf("Error %d opening the file.\n", err);
				break;
			}
		
			// Write the data.
			if (err = OriginalFunctions::WriteVersion2NumericWave(fr, &wh,(void *)data_p, comment.c_str(), strlen(comment.c_str())))
				printf("Error %d writing the file.\n", err);
		
			OriginalFunctions::CPCloseFile(fr);
		} while(0);
		
		delete [] data_p;
/*		if (err == 0)
			printf("Successfully wrote the file %s.\n", filePath);
*/		
		if(err == 0) return true;
		else return false;
		
	};
};





