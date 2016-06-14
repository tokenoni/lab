inline double AndorSIF::ExposureTime()const{  if(!Allocated) return 0.0; return 1.0*InstaImage.exposure_time;}
inline size_t AndorSIF::SequenceSize()const{  if(!Allocated) return 0  ; return sequence_size;}
inline size_t AndorSIF::VerticalSize()const{  if(!Allocated) return 0  ; return vertical_size;}
inline size_t AndorSIF::HorizontalSize()const{if(!Allocated) return 0  ; return horizontal_size;}

inline std::vector<std::vector<std::vector<double>>> AndorSIF::getAllFrames()const{
	std::vector<std::vector<std::vector<double>>> rslt(horizontal_size, std::vector<std::vector<double>>(vertical_size, std::vector<double>(sequence_size, 0.0)));
	for(size_t i=0; i<horizontal_size; ++i){
		for(size_t j=0; j<vertical_size; ++j){
			for(size_t k=0; k<sequence_size; ++k)
				rslt[i][j][k] = data[k*one_image_size + j*horizontal_size + i];
		}
	}
	return rslt;
}


inline 
double AndorSIF::get(size_t sequence_index, size_t vertical_index, size_t horizontal_index)
const{
	return data[sequence_index*one_image_size + vertical_index*horizontal_size + horizontal_index];
}

inline double AndorSIF::operator()(size_t sequence_index, size_t vertical_index, size_t horizontal_index)
const{
	return data[sequence_index*one_image_size + vertical_index*horizontal_size + horizontal_index];
}

inline 
std::vector<double > AndorSIF::get(size_t sequence_index, size_t vertical_index)
const{
	std::vector<double > data_v(horizontal_size, 0);
	size_t offset = sequence_index*one_image_size + vertical_index*horizontal_size;
	for(size_t i=0; i<horizontal_size; ++i)
		data_v[i] = data[offset + i];
	return data_v;
}
inline 
std::vector<std::vector<double >> AndorSIF::get(size_t sequence_index)
const{
	std::vector<std::vector<double >> data_m(horizontal_size, std::vector<double >(vertical_size, 0));
	size_t offset = sequence_index*one_image_size;
	for(size_t i=0; i<horizontal_size; ++i){
		for(size_t j=0; j<vertical_size; ++j){
			data_m[i][j] = data[offset + j*horizontal_size + i];
		}
	}
	return data_m;
}

inline 
std::vector<std::vector<double >> AndorSIF::getHorizontalDevelopment(size_t vertical_index)
const{
	std::vector<std::vector<double >> data_m(horizontal_size, std::vector<double >(sequence_size, 0));
	for(size_t i=0; i<sequence_size; ++i){
		size_t offset = i*one_image_size + vertical_index*horizontal_size;
		for(size_t j=0; j<horizontal_size; ++j)
			data_m[j][i] = data[offset + j];
	}
	return data_m;
}

inline 
std::vector<std::vector<double >> AndorSIF::getVerticalDevelopment(size_t horizontal_index)
const{
	std::vector<std::vector<double >> data_m(vertical_size, std::vector<double >(sequence_size, 0));
	for(size_t i=0; i<sequence_size; ++i){
		size_t offset = i*one_image_size + horizontal_index;
		for(size_t j=0; j<horizontal_size; ++j)
			data_m[j][i] = data[offset + j*horizontal_index];
	}
	return data_m;
}

inline bool AndorSIF::subBackground(const size_t bg_start, const size_t bg_num_){
	if(bg_start>SequenceSize()) return false;
	size_t bg_num=bg_num_;
	if(bg_start+bg_num > SequenceSize()) bg_num = SequenceSize() - bg_start -1;

	std::vector<std::vector<double>> bgData(horizontal_size, std::vector<double>(vertical_size, 0.0));

	//	get background by averaging
	for(size_t count=0; count<bg_num; ++count){
		size_t offset = (bg_start+count)*one_image_size;
		for(size_t i=0; i<horizontal_size; ++i){
			for(size_t j=0; j<vertical_size; ++j){
				bgData[i][j] += data[offset + j*horizontal_size + i]/(1.0*bg_num);
			}
		}
	}

	//	substract the background data
	for(size_t count=0; count<SequenceSize(); ++count){
		size_t offset = count*one_image_size;
		for(size_t i=0; i<horizontal_size; ++i){
			for(size_t j=0; j<vertical_size; ++j){
				data[offset + j*horizontal_size + i] -= bgData[i][j];
			}
		}
	}
	return true;
}


//-------------------------------------------------//
//												   //
//		functions for reading SIF file			   //
//												   //
//-------------------------------------------------//
inline 
bool AndorSIF::read(std::string filename){
	sif_name = filename;
	fin.open(filename.c_str(), std::ios::in | std::ios::binary);
	//	return false if file cannot be open.
	if(!fin) return false;

	std::string title;
	//---	if title is not appropriate	---
	SIFread::read_string(fin, title, '\n');
	if((  title.find("Andor Technology Multi-Channel File" )!=std::string::npos) && 
		( title.find("Oriel Instruments Multi-Channel File")!=std::string::npos)){
			std::cerr<<"This is not a proper SIF file: The file may be corrupt\nPress any key to continue"<<std::endl;
			return false;
	}

	//---	check version	----
	long version = SIFread::read_int(fin, ' ');

	//---	read 5 blocks	---
	for(size_t p=0; p<5; ++p){
		if(!fin.eof()){
			if(SIFread::read_int(fin, '\n')){	// (1 if true, 0 if false)
				read_instaimage();
  				read_calibimage();
  				read_image_structure();
				read_image();
			}
		}
	}
	char ch;
	fin.read(&ch, 1);
	if(fin.eof()){
		fin.close();
		return true;
	}else{
		std::cerr<<"End of file has not been reached\n";
		fin.close();
		return false;
	}
}

inline 
void AndorSIF::read_image(){
	MemoryControl();	//	<---	memory allocation	
	fin.read((char*) data, size*sizeof(float));	//	read the image and contains the data into 1 long vector
	
	//---	remember the matrix size	---
	sequence_size = Image.no_images;
	vertical_size = Image.no_subimages;
	if(vertical_size == 1)
		vertical_size= (Image.position[0].top - Image.position[0].bottom + 1)/Image.position[0].vertical_bin;
//		vertical_size  = Image.image_format.top - Image.image_format.bottom+1;
	if(sequence_size != 0)
		horizontal_size= (Image.position[0].right - Image.position[0].left + 1)/Image.position[0].horizontal_bin;
	one_image_size = vertical_size*horizontal_size;
}

inline 
void AndorSIF::MemoryControl(){
	if(Allocated && size != Image.total_length) FreeMemory();
	size = Image.total_length;
	if(!Allocated){
		data = new float [size];
		Allocated = true;
	}
}
inline 
void AndorSIF::FreeMemory(){
	if(Allocated) delete [] data;
	Allocated = false;
}

inline 
void AndorSIF::read_image_structure(){
	long version1;
	std::vector<long> version2(0);
	int j, k;

	version1				 =	    SIFread::read_int(fin, ' ');
	Image.image_format.left  = (int)SIFread::read_int(fin, ' ');
	Image.image_format.top   = (int)SIFread::read_int(fin, ' ');
	Image.image_format.right = (int)SIFread::read_int(fin, ' ');
	Image.image_format.bottom= (int)SIFread::read_int(fin, ' ');
	//	added by goto
//	channel_number[data_sequence_number] = Image.image_format.right - Image.image_format.left + 1;
	Image.no_images			 = (int)SIFread::read_int(fin, ' ');
//	time_series[data_sequence_number] = Image.no_images;
	Image.no_subimages		 = (int)SIFread::read_int(fin, ' ');
//   subimage_number[data_sequence_number] = Image.no_subimages; 
	Image.total_length		 = (int)SIFread::read_int(fin, ' ');
	Image.image_length		 = (int)SIFread::read_int(fin, '\n');
//	subimage_size[data_sequence_number] = Image.image_length;
	version2.resize(Image.no_subimages);
	Image.position.resize(Image.no_subimages);
	Image.subimage_offset.resize(Image.no_subimages);
	for(j=0;j<Image.no_subimages;j++){ //Repeat no_subimages times
		version2[j]						= (int)SIFread::read_int(fin, ' ');
		Image.position[j].left			= (int)SIFread::read_int(fin, ' ');
		Image.position[j].top			= (int)SIFread::read_int(fin, ' ');
		Image.position[j].right			= (int)SIFread::read_int(fin, ' ');
		Image.position[j].bottom		= (int)SIFread::read_int(fin, ' ');
		Image.position[j].vertical_bin	= (int)SIFread::read_int(fin, ' ');
		Image.position[j].horizontal_bin= (int)SIFread::read_int(fin, ' ');
		Image.subimage_offset[j]		= (int)SIFread::read_int(fin, '\n');
	} // End of for(j=0;j<no_subimages;j++)

	Image.time_stamps.resize(Image.no_images);
	for(k=0;k<Image.no_images;k++){
		Image.time_stamps[k]			= (int)SIFread::read_int(fin, '\n');
	}
	if ((InstaImage.version[0] == 65559) && InstaImage.head_model.find("DV420") != std::string::npos)	SIFread::read_int(fin, '\n');
}

/*****************************************************************************/
//Function: read_calibimage
//The purpose of this function is to read in the calibimage structure data
/*****************************************************************************/
inline 
void AndorSIF::read_calibimage(){
	long version;
	int len;
	version				= SIFread::read_int(fin, ' ');
	CalibImage.x_type   = SIFread::read_byte_and_skip_terminator(fin);
	CalibImage.x_unit   = SIFread::read_byte_and_skip_terminator(fin);
	CalibImage.y_type   = SIFread::read_byte_and_skip_terminator(fin);
	CalibImage.y_unit   = SIFread::read_byte_and_skip_terminator(fin);
	CalibImage.z_type   = SIFread::read_byte_and_skip_terminator(fin);
	CalibImage.z_unit   = SIFread::read_byte_and_skip_terminator(fin);
	CalibImage.x_cal[0] = SIFread::read_float(fin, ' ');
	CalibImage.x_cal[1] = SIFread::read_float(fin, ' ');
	CalibImage.x_cal[2] = SIFread::read_float(fin, ' ');
	CalibImage.x_cal[3] = SIFread::read_float(fin, '\n');
	CalibImage.y_cal[0] = SIFread::read_float(fin, ' ');
	CalibImage.y_cal[1] = SIFread::read_float(fin, ' ');
	CalibImage.y_cal[2] = SIFread::read_float(fin, ' ');
	CalibImage.y_cal[3] = SIFread::read_float(fin, '\n');
	CalibImage.z_cal[0] = SIFread::read_float(fin, ' ');
	CalibImage.z_cal[1] = SIFread::read_float(fin, ' ');
	CalibImage.z_cal[2] = SIFread::read_float(fin, ' ');
	CalibImage.z_cal[3] = SIFread::read_float(fin, '\n');

	CalibImage.rayleigh_wavelength = SIFread::read_float(fin, '\n');
	CalibImage.pixel_length = SIFread::read_float(fin, '\n');
	CalibImage.pixel_height = SIFread::read_float(fin, '\n');

	len = (int)SIFread::read_int(fin, '\n');
	CalibImage.x_text = SIFread::read_len_chars(fin, len);
	len = (int)SIFread::read_int(fin, '\n');
	CalibImage.y_text = SIFread::read_len_chars(fin, len);
	len = (int)SIFread::read_int(fin, '\n');
	CalibImage.z_text = SIFread::read_len_chars(fin, len);
//	print_calib_image(version,CalibImage.x_text,CalibImage.y_text,CalibImage.z_text);
} //End of TCalibImage
/*****************************************************************************/
//Function: read_instaimage
//The purpose of this function is to read in the instaimage structure data
/*****************************************************************************/
inline 
void AndorSIF::read_instaimage(){
	std::string instaimage_title;
//	long version[4];
	long result;
	InstaImage.version[0]		 =	(unsigned int)SIFread::read_int(fin, ' ');
	InstaImage.type  			 =	(unsigned int)SIFread::read_int(fin, ' ');
  	InstaImage.active			 =	(unsigned int)SIFread::read_int(fin, ' ');
  	InstaImage.structure_version =  (unsigned int)SIFread::read_int(fin, ' ');
	InstaImage.timedate			 =                SIFread::read_int(fin, ' ');
	InstaImage.temperature		 =				  SIFread::read_float(fin, ' ');
	InstaImage.head				 = SIFread::read_byte_and_skip_terminator(fin);
  	InstaImage.store_type		 = SIFread::read_byte_and_skip_terminator(fin);
  	InstaImage.data_type		 = SIFread::read_byte_and_skip_terminator(fin);
  	InstaImage.mode				 = SIFread::read_byte_and_skip_terminator(fin);
  	InstaImage.trigger_source 	 = SIFread::read_byte_and_skip_terminator(fin);
	InstaImage.trigger_level	 = SIFread::read_float(fin, ' ');
  	InstaImage.exposure_time	 = SIFread::read_float(fin, ' ');
  	InstaImage.delay			 = SIFread::read_float(fin, ' ');
  	InstaImage.integration_cycle_time = SIFread::read_float(fin, ' ');
	InstaImage.no_integrations   = (int)SIFread::read_int(fin, ' ');
  	InstaImage.sync				 = SIFread::read_byte_and_skip_terminator(fin);
  	InstaImage.kinetic_cycle_time= SIFread::read_float(fin, ' ');
  	InstaImage.pixel_readout_time= SIFread::read_float(fin, ' ');
	InstaImage.no_points		 = (int)SIFread::read_int(fin, ' ');
  	InstaImage.fast_track_height = (int)SIFread::read_int(fin, ' ');
	InstaImage.gain				 = (int)SIFread::read_int(fin, ' ');
	InstaImage.gate_delay		 = SIFread::read_float(fin, ' ');
  	InstaImage.gate_width		 = SIFread::read_float(fin, ' ');

  	if( (HIWORD(InstaImage.version[0]) >= 1) && (LOWORD(InstaImage.version[0]) >=6) )
		InstaImage.GateStep = SIFread::read_float(fin, ' ');
	
	InstaImage.track_height		 = (int)SIFread::read_int(fin, ' ');
	InstaImage.series_length	 = (int)SIFread::read_int(fin, ' ');
	InstaImage.read_pattern		 = SIFread::read_byte_and_skip_terminator(fin);
	InstaImage.shutter_delay	 = SIFread::read_byte_and_skip_terminator(fin);

	if( (HIWORD(InstaImage.version[0]) >= 1) && (LOWORD(InstaImage.version[0]) >=7) ) {
	  InstaImage.st_centre_row	 = (int)SIFread::read_int(fin, ' ');
	  InstaImage.mt_offset		 = (int)SIFread::read_int(fin, ' ');
	  InstaImage.operation_mode  = (int)SIFread::read_int(fin, ' ');
	}

	if( (HIWORD(InstaImage.version[0]) >= 1) && (LOWORD(InstaImage.version[0]) >=8) ) {
	  InstaImage.FlipX			 = (int)SIFread::read_int(fin, ' ');
	  InstaImage.FlipY			 = (int)SIFread::read_int(fin, ' ');
	  InstaImage.Clock			 = (int)SIFread::read_int(fin, ' ');
	  InstaImage.AClock			 = (int)SIFread::read_int(fin, ' ');
	  InstaImage.MCP			 = (int)SIFread::read_int(fin, ' ');
	  InstaImage.Prop			 = (int)SIFread::read_int(fin, ' ');
	  InstaImage.IOC			 = (int)SIFread::read_int(fin, ' ');
	  InstaImage.Freq		 = (float)1.0*SIFread::read_int(fin, ' ');
	}

	if((HIWORD(InstaImage.version[0]) >= 1) && (LOWORD(InstaImage.version[0]) >=9)) {
	  InstaImage.VertClockAmp	 = (int)SIFread::read_int(fin, ' ');
	  InstaImage.data_v_shift_speed = SIFread::read_float(fin, ' ');
	}

	if((HIWORD(InstaImage.version[0]) >= 1) && (LOWORD(InstaImage.version[0]) >=10)) {
	  InstaImage.OutputAmp		 = (int)SIFread::read_int(fin, ' ');
	  InstaImage.PreAmpGain		 = SIFread::read_float(fin, ' ');
	}

	if((HIWORD(InstaImage.version[0]) >= 1) && (LOWORD(InstaImage.version[0]) >=11)) {
	  InstaImage.Serial			 = (int)SIFread::read_int(fin, ' ');
	}

	if((HIWORD(InstaImage.version[0]) >= 1) && (LOWORD(InstaImage.version[0]) >=13)) {
	  InstaImage.NumPulses		 = (int)SIFread::read_int(fin, ' ');
	}

	if((HIWORD(InstaImage.version[0]) >= 1) && (LOWORD(InstaImage.version[0]) >=14)) {
	  InstaImage.mFrameTransferAcqMode  = (int)SIFread::read_int(fin, ' ');
	}

	if((HIWORD(InstaImage.version[0]) >= 1) && (LOWORD(InstaImage.version[0]) >=5) ){
		SIFread::read_int(fin, '\n');//
		InstaImage.head_model = SIFread::read_string(fin, '\n');
		InstaImage.detector_format_x = (int)SIFread::read_int(fin, ' ');
		InstaImage.detector_format_z = (int)SIFread::read_int(fin, ' ');
	}
	else if( (HIWORD(InstaImage.version[0]) >= 1) && (LOWORD(InstaImage.version[0]) >=3) ){
		unsigned long head_model = (unsigned long)SIFread::read_int(fin, ' ');
		InstaImage.head_model	 = SIFread::add_to_string(std::string(""), head_model);
		InstaImage.detector_format_x = (int)SIFread::read_int(fin, ' ');
		InstaImage.detector_format_z = (int)SIFread::read_int(fin, ' ');
	}
	else {
		InstaImage.head_model = "Unknown";
		InstaImage.detector_format_x = 1024u;
		InstaImage.detector_format_z = 256u;
	}

	SIFread::read_int(fin, '\n');
	InstaImage.filename = SIFread::read_string(fin, '\n');

	//Start of TUserText
	InstaImage.version[1]= SIFread::read_int(fin, ' ');
	result				 = SIFread::read_int(fin, '\n');
//	InstaImage.user_text= SIFread::read_len_chars(fin, (int)result);
	InstaImage.user_text= SIFread::read_string(fin, '\n');
	//End of TUserText

	//Start of TShutter
	if( (HIWORD(InstaImage.version[0]) >= 1) && (LOWORD(InstaImage.version[0]) >=4) ) {
		InstaImage.version[2]			 = SIFread::read_int(fin, ' ');
		InstaImage.shutter.type			 = SIFread::read_byte_and_skip_terminator(fin);
		InstaImage.shutter.mode			 = SIFread::read_byte_and_skip_terminator(fin);
		InstaImage.shutter.custom_bg_mode= SIFread::read_byte_and_skip_terminator(fin);
		InstaImage.shutter.custom_mode	 = SIFread::read_byte_and_skip_terminator(fin);
		InstaImage.shutter.closing_time  = SIFread::read_float(fin, ' ');
		InstaImage.shutter.opening_time  = SIFread::read_float(fin, '\n');
	}
	// End of TShutter
	if ((InstaImage.version[0] >= 65548) && (InstaImage.version[0] <= 65557)) {
		for (size_t i = 0; i < 2; i++)  SIFread::read_string(fin, '\n');
	}

	if (InstaImage.version[0] == 65558) {
		for (size_t i = 0; i < 5; i++) SIFread::read_string(fin, '\n');
	}

	if (InstaImage.version[0] >= 65559) {
		for (size_t i = 0; i < 9; i++)	SIFread::read_string(fin, '\n');
	}
  
	//Start of TShamrockSave
/*	if( (HIWORD(version[0])>1) || ((HIWORD(version[0])==1) && (LOWORD(version[0]) >=12)) ){
		version[3]								 = SIFread::read_int(fin, ' ');
		InstaImage.shamrock_save.IsActive		 = SIFread::read_int(fin, ' ');
		InstaImage.shamrock_save.WavePresent	 = SIFread::read_int(fin, ' ');
		InstaImage.shamrock_save.Wave			 = SIFread::read_float(fin, ' ');
		InstaImage.shamrock_save.GratPresent	 = SIFread::read_int(fin, ' ');
		InstaImage.shamrock_save.GratIndex		 = SIFread::read_int(fin, ' ');
		InstaImage.shamrock_save.GratLines		 = SIFread::read_float(fin, ' ');
		InstaImage.shamrock_save.GratBlaze		 = SIFread::read_string(fin, '\n');
		InstaImage.shamrock_save.SlitPresent	 = SIFread::read_int(fin, ' ');
		InstaImage.shamrock_save.SlitWidth		 = SIFread::read_float(fin, ' ');
		InstaImage.shamrock_save.FlipperPresent  = SIFread::read_int(fin, ' ');
		InstaImage.shamrock_save.FlipperPort	 = SIFread::read_int(fin, ' ');
		InstaImage.shamrock_save.FilterPresent	 = SIFread::read_int(fin, ' ');
		InstaImage.shamrock_save.FilterIndex	 = SIFread::read_int(fin, ' ');
		len										 = SIFread::read_int(fin, ' ');
		InstaImage.shamrock_save.FilterString	 = SIFread::read_len_chars(fin, len);
		InstaImage.shamrock_save.AccessoryPresent= SIFread::read_int(fin, ' ');
		InstaImage.shamrock_save.Port1State		 = SIFread::read_int(fin, ' ');
		InstaImage.shamrock_save.Port2State		 = SIFread::read_int(fin, ' ');
		InstaImage.shamrock_save.Port3State		 = SIFread::read_int(fin, ' ');
		InstaImage.shamrock_save.Port4State		 = SIFread::read_int(fin, ' ');
		InstaImage.shamrock_save.OutputSlitPresent = SIFread::read_int(fin, ' ');
		InstaImage.shamrock_save.OutputSlitWidth   = SIFread::read_float(fin, ' ');
	}
	//End of TShamrockSave

//	print_instaimage(version,InstaImage.user_text.text);
*/
};