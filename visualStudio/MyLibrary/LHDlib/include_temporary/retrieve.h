#ifndef _RETRIEVE_H_
#define _RETRIEVE_H_

#if defined(WIN32) || defined(WIN64)
 #ifndef _GNUC__ 
  #pragma once

// ���Ŏw�肳�ꂽ��`�̑O�ɑΏۃv���b�g�t�H�[�����w�肵�Ȃ���΂Ȃ�Ȃ��ꍇ�A�ȉ��̒�`��ύX���Ă��������B
// �قȂ�v���b�g�t�H�[���ɑΉ�����l�Ɋւ���ŐV���ɂ��ẮAMSDN ���Q�Ƃ��Ă��������B
  #ifndef WINVER				// Windows XP �ȍ~�̃o�[�W�����ɌŗL�̋@�\�̎g�p�������܂��B
   #define WINVER 0x0501		// ����� Windows �̑��̃o�[�W���������ɓK�؂Ȓl�ɕύX���Ă��������B
  #endif

  #ifndef _WIN32_WINNT		// Windows XP �ȍ~�̃o�[�W�����ɌŗL�̋@�\�̎g�p�������܂��B
   #define _WIN32_WINNT 0x0501	// ����� Windows �̑��̃o�[�W���������ɓK�؂Ȓl�ɕύX���Ă��������B
  #endif						

  #ifndef _WIN32_WINDOWS		// Windows 98 �ȍ~�̃o�[�W�����ɌŗL�̋@�\�̎g�p�������܂��B
   #define _WIN32_WINDOWS 0x0410 // ����� Windows Me �܂��͂���ȍ~�̃o�[�W���������ɓK�؂Ȓl�ɕύX���Ă��������B
  #endif

  #ifndef _WIN32_IE			// IE 6.0 �ȍ~�̃o�[�W�����ɌŗL�̋@�\�̎g�p�������܂��B
   #define _WIN32_IE 0x0600	// ����� IE. �̑��̃o�[�W���������ɓK�؂Ȓl�ɕύX���Ă��������B
  #endif

  #define WIN32_LEAN_AND_MEAN
 #endif
 #include <windows.h>
 #ifdef _MSC_VER
  #ifdef RETRIEVE_EXPORTS
   #define RETRIEVE_API __declspec(dllexport)
  #else
   #define RETRIEVE_API __declspec(dllimport)
  #endif
 #else
  #define RETRIEVE_API 
 #endif

 #define forretrieveopenwait_( x1,x2,x3,x4,x5 )			__stdcall FORRETRIEVEOPENWAIT( x1,x2,x3,x4,x5 )
 #define forretrieveopendirectwait_( x1,x2,x3,x4,x5,x6,x7 )	__stdcall FORRETRIEVEOPENDIRECTWAIT( x1,x2,x3,x4,x5,x6,x7 )
 #define forretrieveopen_( x1,x2,x3,x4 )			__stdcall FORRETRIEVEOPEN( x1,x2,x3,x4 )
 #define forretrieveopendirect_( x1,x2,x3,x4,x5,x6 )	__stdcall FORRETRIEVEOPENDIRECT( x1,x2,x3,x4,x5,x6 )
 #define forretrieveshotinfo_( x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12 )	\
													__stdcall FORRETRIEVESHOTINFO( x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12 )
 #define forretrievechinfo_( x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12 )	\
													__stdcall FORRETRIEVECHINFO( x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12 )
 #define forretrievechdata_( x1,x2,x3,x4,x5 )		__stdcall FORRETRIEVECHDATA( x1,x2,x3,x4,x5 )
 #define forretrievechparam_( x1,x2,x3,x4,x5,x6 )		__stdcall FORRETRIEVECHPARAMS( x1,x2,x3,x4,x5,x6 )
 #define forretrieveframeinfo_( x1,x2,x3,x4,x5,x6,x7 )	\
													__stdcall FORRETRIEVEFRAMEINFO( x1,x2,x3,x4,x5,x6,x7 )
 #define forretrieveframedata_( x1,x2,x3,x4,x5,x6 )	__stdcall FORRETRIEVEFRAMEDATA( x1,x2,x3,x4,x5,x6 )
 #define forretrieveclose_( x1 )					__stdcall FORRERIEVECLOSE( x1 )
 #define forretrieveversion_( )						__stdcall FORRETRIEVEVERSION( )
 #define forretrieveErrorMessage_( x1 )				__stdcall FORRETRIEVEERRORMESSAGE( x1 )

 #define forretrievechangesite_( x1 )				__stdcall FORRETRIEVECHANGESITE( x1 )
 #define forretrievegetsite_( )						__stdcall FORRETRIEVEGETSITE( )
 #define forretrievegetdtsinfo_( x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23 )	\
		__stdcall FORRETRIEVEGETDTSINFO( x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23 )
 #define forretrievegetdtsdata_( x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16 )	\
		__stdcall FORRETRIEVEGETDTSDATA( x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16 )

 #define forretrievechvolts_( x1,x2,x3,x4,x5 )		__stdcall FORRETRIEVECHVOLTS( x1,x2,x3,x4,x5 )
 #define forretrievechvoltsdbl_( x1,x2,x3,x4,x5 )	__stdcall FORRETRIEVECHVOLTSDBL( x1,x2,x3,x4,x5 )
 #define forretrieveframeinfo2_( x1,x2,x3,x4,x5,x6,x7 )	\
													__stdcall FORRETRIEVEFRAMEINFO2( x1,x2,x3,x4,x5,x6,x7 )
  typedef __int64			int64_t;
#else
 typedef void* LPVOID;  
 #define RETRIEVE_API
 #if defined(__LP64__)
  typedef long			int64_t;
 #else
  typedef long long		int64_t;
 #endif
#endif

#ifndef INDEX_ENV_UNDEFINED
 #define INDEX_ENV_UNDEFINED -101
#endif


#if defined (__linux__) || defined(__MACH__)
 #ifndef	TRUE
  #define	TRUE	true
 #endif
 #ifndef	FALSE
  #define	FALSE	false
 #endif
 #ifndef	_MAX_PATH				// 2011-06-16 s.i
  #define	_MAX_PATH	256
 #endif
 #ifndef	NAME_SEPARATOR
  #define	NAME_SEPARATOR	'/'
 #endif
#else
 #ifndef	NAME_SEPARATOR
  #define	NAME_SEPARATOR	'\\'
 #endif
#endif

#if defined(WIN64)
 typedef __int64			pvw_long;
 typedef __int64			pvw_ulong;
#elif defined(__LP64__)
 typedef long				pvw_long;
 typedef unsigned long	pvw_ulong;
#else
 typedef int				pvw_long;
 typedef unsigned int	pvw_ulong;
#endif
 typedef unsigned short  pvw_uint;
 typedef short  			pvw_int;

#ifdef __cplusplus
extern "C" {
#endif
//  ----- PV-WAVE -----
RETRIEVE_API pvw_long PvwRetrieveOpen(int argc, LPVOID **argv);
RETRIEVE_API pvw_long PvwRetrieveOpenDirect(int argc, LPVOID **argv);
RETRIEVE_API pvw_long PvwRetrieveShotInfo(int argc, LPVOID **argv);
RETRIEVE_API pvw_long PvwRetrieveChInfo(int argc, LPVOID **argv);
RETRIEVE_API pvw_long PvwRetrieveChData(int argc, LPVOID **argv);
RETRIEVE_API pvw_long PvwRetrieveChParams(int argc, LPVOID **argv);
RETRIEVE_API pvw_long PvwRetrieveFrameInfo(int argc, LPVOID **argv);
RETRIEVE_API pvw_long PvwRetrieveFrameData(int argc, LPVOID **argv);
RETRIEVE_API pvw_long PvwRetrieveClose(int argc, LPVOID **argv);
RETRIEVE_API const char* PvwRetrieveVersion(int argc, LPVOID **argv);
RETRIEVE_API const char* PvwRetrieveErrorMessage(int argc, LPVOID **argv);
RETRIEVE_API pvw_long PvwRetrieveChangeSite(int argc, LPVOID **argv);
RETRIEVE_API const char* PvwRetrieveGetSite(int argc, LPVOID **argv);

RETRIEVE_API pvw_long PvwRetrieveGetDTSdata(int argc,LPVOID **argv);
RETRIEVE_API pvw_long PvwRetrieveGetDTSdataDBL(int argc,LPVOID **argv);
RETRIEVE_API pvw_long PvwRetrieveGetDTSInfoFromParams(int argc,LPVOID **argv);	//	add 2011-03-23

RETRIEVE_API pvw_long PvwRetrieveGetDTSParameters(int argc, LPVOID **argv);	// add 2010-09-01
RETRIEVE_API pvw_long PvwRetrieveGetDTSParametersCount(int argc, LPVOID **argv);	// add 2010-09-01

RETRIEVE_API pvw_long PvwRetrieveRealShotNumber(int argc, LPVOID **argv);	// add 2009-09-14
RETRIEVE_API pvw_long PvwRetrieveChV(int argc, LPVOID **argv);				// add 2010-03-05
RETRIEVE_API pvw_long PvwRetrieveChVolts(int argc, LPVOID **argv);			// add 2010-03-05
RETRIEVE_API pvw_long PvwRetrieveChVoltsDbl(int argc, LPVOID **argv);		// add 2010-03-05

RETRIEVE_API pvw_long PvwRetrieveGetDTSdatax2(int argc,LPVOID **argv);				// add 2012-07-13
RETRIEVE_API pvw_long PvwRetrieveGetDTSinfox(int argc, LPVOID **argv);// 2012-07-13
RETRIEVE_API pvw_long PvwRetrieveGetDTSinfoFromRetrieve(int argc,LPVOID **argv);// 2012-07-13
RETRIEVE_API pvw_long PvwRetrieveGetDTSinfoFromParameters(int argc,LPVOID **argv);// 2012-07-13
RETRIEVE_API pvw_long PvwIndexNoClosing(int argc,LPVOID **argv);// 2012-07-13

#if defined(WIN64) || defined(__LP64__)
RETRIEVE_API pvw_long PvwLong64From32(int argc,LPVOID **argv);
void pvw_long64From32( pvw_long* dst,int* src,unsigned int len);
#endif
RETRIEVE_API pvw_long PvwRGBfromYUV422(int argc,LPVOID **argv);
RETRIEVE_API pvw_long PvwIntgerFromByte(int argc, LPVOID **argv);
RETRIEVE_API pvw_long PvwChannelDecode(int argc, LPVOID **argv);

//  ----- IDL -----
RETRIEVE_API int IdlRetrieveOpen(int argc, LPVOID argv[]);
RETRIEVE_API int IdlRetrieveOpenDirect(int argc, LPVOID argv[]);
RETRIEVE_API int IdlRetrieveShotInfo(int argc, LPVOID argv[]);
RETRIEVE_API int IdlRetrieveChInfo(int argc, LPVOID argv[]);
RETRIEVE_API int IdlRetrieveChData(int argc, LPVOID argv[]);
RETRIEVE_API int IdlRetrieveChParam(int argc, LPVOID argv[]);
RETRIEVE_API int IdlRetrieveFrameInfo(int argc, LPVOID argv[]);
RETRIEVE_API int IdlRetrieveFrameData(int argc, LPVOID argv[]);
RETRIEVE_API int IdlRetrieveClose(int argc, LPVOID argv[]);
RETRIEVE_API const char* IdlRetrieveVersion(int argc, LPVOID argv[]);
RETRIEVE_API const char* IdlRetrieveErrorMessage(int argc, LPVOID argv[]);
RETRIEVE_API int IdlRetrieveChangeSite(int argc, LPVOID argv[]);
RETRIEVE_API const char* IdlRetrieveGetSite(int argc, LPVOID argv[]);

RETRIEVE_API int IdlRetrieveGetDTSdata(int argc, LPVOID argv[]);
RETRIEVE_API int IdlRetrieveGetDTSdataDBL(int argc, LPVOID argv[]);
RETRIEVE_API int IdlRetrieveGetDTSInfoFromParams(int argc, LPVOID argv[]);	//	add 2011-05-12

RETRIEVE_API int IdlRetrieveGetDTSParameters(int argc, LPVOID **argv);	// add 2010-10-01
RETRIEVE_API int IdlRetrieveGetDTSParametersCount(int argc, LPVOID **argv);	// add 2010-10-01

RETRIEVE_API int IdlRetrieveRealShotNumber( int argc, LPVOID argv[] );	// add 2009-09-14
RETRIEVE_API int IdlRetrieveChV(int argc, LPVOID argv[]);				// add 2010-03-05
RETRIEVE_API int IdlRetrieveChVolts(int argc, LPVOID argv[]);			// add 2010-03-05
RETRIEVE_API int IdlRetrieveChVoltsDbl(int argc, LPVOID argv[]);		// add 2010-03-05
RETRIEVE_API int IdlIntgerFromByte(int argc, LPVOID argv[]);
RETRIEVE_API int IdlChannelDecode(int argc, LPVOID argv[]);


RETRIEVE_API int IdlRetrieveGetDTSinfox(int argc,LPVOID argv[]);				// 2012-07-13
RETRIEVE_API int IdlRetrieveGetDTSinfoFromRetrieve(int argc,LPVOID argv[]);		// 2012-07-13
RETRIEVE_API int IdlRetrieveGetDTSinfoFromParameters(int argc,LPVOID argv[]);	// 2012-07-13
RETRIEVE_API int IdlRetrieveGetDTSdatax2(int argc,LPVOID argv[]);		// 2012-07-13
RETRIEVE_API int IdlIndexNoClosing(int argc,LPVOID argv[]);		// 2012-07-13
RETRIEVE_API int IdlRetrieveChDTSparameters(int argc, LPVOID argv[]);

//  ----- Fortrun -----
RETRIEVE_API int forretrieveopenwait_( char *diag_name ,char* index_server,
							   unsigned int shot_number, unsigned short sub_shot ,int timeout_sec );
RETRIEVE_API int forretrieveopendirectwait_( char *diag_name , char* server ,char* path,
							   unsigned int shot_number, unsigned short sub_shot ,int port ,int timeout_sec);
RETRIEVE_API int forretrieveopen_( char *diag_name ,char* index_server,
							   unsigned int shot_number, unsigned short sub_shot );
RETRIEVE_API int forretrieveopendirect_( char *diag_name , char* server ,char* path,
							   unsigned int shot_number, unsigned short sub_shot ,int port);
RETRIEVE_API int forretrieveshotinfo_(int rtr_d ,unsigned int* n_channel ,
							   short* year ,short* month ,short* day ,short* hour ,short* min ,short* sec ,
							   char *management ,char* comment ,int comment_size ,char *server ); 
RETRIEVE_API int forretrievechinfo_( int rtr_d ,int ch_no ,
			unsigned int* data_length ,unsigned int* comp_length ,unsigned short* param_count ,
			short* data_type ,char* image_type ,unsigned short* value_len ,int *is_frame ,
			char* management ,char* comment ,int comment_size );
RETRIEVE_API int forretrievechdata_( int rtr_d, int ch_no ,
			char* data ,unsigned int array_size ,unsigned int* data_length );
RETRIEVE_API int forretrievechparam_( int rtr_d, int ch_no , int param_no ,
			char* param_name ,char* param_value ,int* param_type );
RETRIEVE_API int forretrieveframeinfo_( int rtr_d ,int ch_no ,int frame_no ,
			unsigned int* data_length ,char* image_type ,unsigned int* frame_x ,unsigned int* frame_y ); 
RETRIEVE_API int forretrieveframedata_( int rtr_d, int ch_no ,int frame_no ,
			char* data ,unsigned int array_size ,unsigned int* data_length );
RETRIEVE_API int forretrieveclose_(int rtr_d);
RETRIEVE_API const char* forretrieveversion_( );
RETRIEVE_API const char* forretrieveerrormessage_( int error_code );
RETRIEVE_API int forretrievechangesite_( char *site_name);
RETRIEVE_API const char* forretrievegetsite_( );

RETRIEVE_API int forretrievegetdtsinfo_
					(const char *DiagName,							// �v����
					 unsigned int AShotNumber,						// �₢���킹�V���b�g�ԍ�
					 short ASubShotNumber,							// �₢���킹�T�u�V���b�g�ԍ�
					 unsigned int *ShotNumber,						// �V���b�g�ԍ�
					 unsigned short *SubShotNumber,					// �T�u�V���b�g�ԍ�
					 unsigned short *LastChannel,					// �܂܂��`���l���f�[�^��
					 short *Year,									// �N
					 short *Month,									// ��
					 short *Day,									// ��
					 short *Hour,									// ��
					 short *Minute,									// ��
					 short *Second,									// �b

					 char  *DTSsource,
					 char  *DTShostID,
					 char  *DTSmoduleID,
					 short *DTStriggerChannel,
					 short *DTSclockChannel,
					 int  *DTSuserDefine,
					 int  *DTStimeArraySize,

					 const char *IndexServer,					// �C���f�b�N�X�T�[�o
					 short *Fastch,								// �Ŏn�`���l���ԍ�
					 short *Lastch,								// �ŏI�`���l���ԍ�
					 char *clockSource);

RETRIEVE_API int forretrievegetdtsdata_(char *IndexServer,
							short *ch,
							float *timeArray,
							char *DTSsource,
							char *DTShostID,
							char *DTSmoduleID,
							short *DTStriggerChannel,
							short *DTSclockChannel,
							int *DTSuserDefine,
							int  DTStimeArraySize,
							unsigned int shot,
							unsigned short subshot,
							char *clockSource,
							short bSilent,

							float *clockCycle,
							float *triggerTiming);

RETRIEVE_API int forretrievegetdtsdatadbl_(char *IndexServer,
							short *ch,
							double *timeArray,
							char *DTSsource,
							char *DTShostID,
							char *DTSmoduleID,
							short *DTStriggerChannel,
							short *DTSclockChannel,
							int *DTSuserDefine,
							int  DTStimeArraySize,
							unsigned int shot,
							unsigned short subshot,
							char *clockSource,
							short bSilent,

							double *clockCycle,
							double *triggerTiming);

RETRIEVE_API int forretrievechvolts_( int rtr_d, int ch_no ,
			float* data ,unsigned int array_size ,unsigned int* data_length );
RETRIEVE_API int forretrievechvoltsdbl_( int rtr_d, int ch_no ,
			double* data ,unsigned int array_size ,unsigned int* data_length );
RETRIEVE_API int forretrieveframeinfo2_( int rtr_d ,int ch_no ,int frame_no ,
			unsigned int* data_length ,char* image_type ,unsigned int* frame_x ,unsigned int* frame_y ); 


/**
 * 
 * @date 2010-11-15
 * @note int retrieveChData( ...,char* data ,..) -> void* data
 * @note int retrieveFrameData( ...,char* data ,..) -> void* data
 */
//  ----- C/C++ -----
RETRIEVE_API int retrieveOpenWait( const char *diag_name ,const char* index_server ,
							   unsigned int shot_number, unsigned short sub_shot ,int timeout_sec );
RETRIEVE_API int retrieveOpenDirectWait( const char *diag_name ,const char* server ,const char* path,
							   unsigned int shot_number, unsigned short sub_shot ,int port ,int timeout_sec );
RETRIEVE_API int retrieveOpen( const char *diag_name ,const char* index_server ,
							   unsigned int shot_number, unsigned short sub_shot );
RETRIEVE_API int retrieveOpenDirect( const char *diag_name ,const char* server ,const char* path,
							   unsigned int shot_number, unsigned short sub_shot ,int port );
RETRIEVE_API int retrieveShotInfo(int rtr_d ,unsigned int* n_channel ,
							   short* year ,short* month ,short* day ,short* hour ,short* min ,short* sec ,
							   char *management ,char* comment ,int comment_size ,char *server ); 
RETRIEVE_API int retrieveChInfo( int rtr_d ,int ch_no ,
			unsigned int* data_length ,unsigned int* comp_length ,unsigned short* param_count ,
			short* data_type ,char* image_type ,unsigned short* value_len ,int* is_nframe ,
			char* management ,char* comment ,int comment_size );
RETRIEVE_API int retrieveChData( int rtr_d, int ch_no,
			void* data, unsigned int array_size, unsigned int* data_length );
RETRIEVE_API int retrieveChParams( int rtr_d, int ch_no ,
			char** param_name ,char** param_value ,int* param_type );
RETRIEVE_API int retrieveChParam( int rtr_d, int ch_no , int param_no ,
			char* param_name ,char* param_value ,int* param_type );
RETRIEVE_API int retrieveFrameInfo( int rtr_d ,int ch_no ,int frame_no ,
			unsigned int* data_length ,char* image_type ,unsigned int* frame_x ,unsigned int* frame_y ); 
RETRIEVE_API int retrieveFrameData( int rtr_d, int ch_no ,int frame_no ,
			void* data ,unsigned int array_size ,unsigned int* data_length );

RETRIEVE_API int retrieveClose(int rtr_d);
RETRIEVE_API const char* retrieveVersion( void );
RETRIEVE_API const char* retrieveErrorMessage( int error_code );
RETRIEVE_API int retrieveChangeSite( const char* site_name );
RETRIEVE_API const char* retrieveGetSite( void );
RETRIEVE_API void retrieveResetIndexServerName( void );
RETRIEVE_API int RGBfromYUV422( unsigned char* rgb, unsigned char* yuv, int yuv_blen,int model);

RETRIEVE_API int retrieveRealShotNumber( int rtr_d, unsigned int* shot_number, unsigned short* sub_shot); // add 2009-09-14
RETRIEVE_API int retrieveChV( int rtr_d, int ch_no ,
			void* data ,unsigned int array_size ,unsigned int* data_length ,int array_type);	// add 2010-03-05
RETRIEVE_API int retrieveChVolts( int rtr_d, int ch_no ,
			float* data ,unsigned int array_size ,unsigned int* data_length );					// add 2010-03-05
RETRIEVE_API int retrieveChVoltsDbl( int rtr_d, int ch_no ,
			double* data ,unsigned int array_size ,unsigned int* data_length );					// add 2010-03-05
RETRIEVE_API const char* retrieveZversion( void );												// add 2010-04-02

RETRIEVE_API int retrieveChangeSite2( const char* site_name, char* env_string, char* env_value);	// add _ to name 2009/01/29
RETRIEVE_API const char* retrieveGetSite( void );
RETRIEVE_API void retrieveGetSite2( char* site_name );

RETRIEVE_API int retrieveFrameData2( int rtr_d, int ch_no ,int frame_no ,
			void* data ,unsigned int array_size ,unsigned int* data_length );

RETRIEVE_API int retrieveChInfo2( int rtr_d, const char* signal_name,
			unsigned int* data_length, unsigned int* comp_length, unsigned short* param_count,
			short* data_type, char* image_type, unsigned short* value_len, int* is_nframe,
			char* management, char* comment, int comment_size, int* ch_no );

RETRIEVE_API int ChannelDecode( const char* str, short* start_ch, short* end_ch);
RETRIEVE_API void retrieve_SetTestPort( int port);

RETRIEVE_API bool retrieve_check_string(const char* str);


RETRIEVE_API int IndexNoClosing( int sw );


RETRIEVE_API int retrieveGetDTSinfox
					(const char *DiagName,						// [in] �v����
					 unsigned int AShotNumber,					// [in] �₢���킹�V���b�g�ԍ�
					 short ASubShotNumber,						// [in] �₢���킹�T�u�V���b�g�ԍ�
					 short Fastch,								// [in] �Ŏn�`���l���ԍ�
					 short Lastch,								// [in] �ŏI�`���l���ԍ�
					 const char *IndexServer,					// [in] �C���f�b�N�X�T�[�o

					 char*	DTSsource,
					 char*	DTShostID,
					 char*	DTSmoduleID,
					 short*	DTStriggerChannel,
					 short*	DTSclockChannel,
					 int*	DTSuserDefine,
					 int*	DTStimeArraySize,
					 char*	ClockSource,
					 int* 	PreSamplings,
					 int* 	SamplingInterval,
					 int*	ClkDTSid,
					 char*	ClkDTSmodule
					);

// 2012-07-10 s.imazu
RETRIEVE_API int retrieveGetDTSdatax2(char* IndexServer,

						char* DTSsource,
						char* DTShostID,
						char* DTSmoduleID,
						char* DTStriggerChannel,

						char* clkSource,
						char* clkHostID,
						char* clkModuleID,
						char* DTSclockChannel,

						int* DTSuserDefine,
						int  DTStimeArraySize,
						unsigned int shot,
						unsigned short subshot,
						char* clockSource,
						char* InternalClock,

						char* strSamplingInterval,
						char* strPreSampling,
						int64_t	lastCount,
						short bSilent,
						int is_double,
						void* timeArray,
						void* clockCycle,
						void* triggerTiming
						);
RETRIEVE_API int retrieveGetDTSinfoFromRetrieve(
			int rtr_d, int type, const void* channel, int needDtslink, int* no_of_channel, 
			unsigned int* ShotNumber, unsigned short* SubShotNumber,
			char* DTSsource, char* DTShostID, char* DTSmoduleID, char* DTStriggerChannel,
			char* CLKsource, char* CLKhostID, char* CLKmoduleID, char* DTSclockChannel,
			int* DTSuserDefine, int* DTStimeArraySize, char* ExtOrInt,
			char* InternalClock, char* SamplingInterval, char* PreSampling,	int64_t* LastCount);

RETRIEVE_API int GetDTSinfoFromChParameters(
				unsigned short p_count,	char* pkey_tbl[], char* pval_tbl[],int needDtslink,
				unsigned int *ShotNumber, unsigned short *SubShotNumber,
				char* DTSsource, char* DTShostID, char* DTSmoduleID, char* DTStriggerChannel,
				char* CLKsource, char* CLKhostID, char* CLKmoduleID, char* DTSclockChannel,
				int* DTSuserDefine, int* DTStimeArraySize, char* ExtOrInt,
				char*InternalClock, char* SamplingInterval, char*PreSampling, int64_t* lastCount);

RETRIEVE_API int getParamValueFromList(
						const char *paramList,
						const char *pName,
						char *pValue);


RETRIEVE_API int retrieveGetDTSParameters(char *IndexServer,
							char *hostname,
							char *modname,
							unsigned int shot,
							unsigned short subshot,
							unsigned short channel,
							unsigned short paramterCount,
							char **param_name,
							char **param_value,
							int *param_type,
#ifdef __cplusplus
 #if defined(WIN32) || defined(WIN64)
							bool IDL_flag = FALSE);
 #else
							bool IDL_flag = false);
 #endif
#else
							bool IDL_flag);
#endif
// add 2013-05-27
RETRIEVE_API int LvGetDTSParameters(char* IndexServer,
							char* hostname,
							char* modname,
							unsigned int shot,
							unsigned short subshot,
							unsigned short channel,
							char* parameters_string,
							int parameters_size,
							int* parameters_count);

RETRIEVE_API int retrieveGetDTSParametersCount(char *IndexServer,
							char *hostname,
							char *modname,
							unsigned int shot,
							unsigned short subshot,
							short* param_count,
							int* start_ch,
							int* end_ch,
#ifdef __cplusplus
							char *true_module_name = NULL);
#else
							char *true_module_name);
#endif

RETRIEVE_API int retrieve_tInfo
					(const char *diag_name,
					 const char* index_server,
					 const char* setup_server,
					 unsigned int shot_number,
					 unsigned short sub_shot, 
					 unsigned short channel,
					 
					 char  *DTSsource,
					 char  *DTShostID,
					 char  *DTSmoduleID,
					 short *DTStriggerChannel,
					 short *DTSclockChannel,
					 unsigned int  *DTSuserDefine,
					 unsigned int  *DTStimeArraySize,

					 unsigned int	*preSamplingNum,
					 double	*trigerTiming_s,
					 unsigned int	*triggerTiming_ns,
					 double	*cycleTime_s,
					 unsigned int	*cycleTime_ns,
					 unsigned int	*cycleTime_Hz
					);

RETRIEVE_API int retrieveGetDTSParameters2(char *IndexServer,
							char *hostname,
							char *modname,
							unsigned int shot,
							unsigned short subshot,
							unsigned short channel,
							unsigned short paramterCount,
							void *param_name,
							void *param_value,
							int *param_type);
RETRIEVE_API void bufferStringCopy(char* string,void* buffer,int element_size,int number);
RETRIEVE_API void set_retrieve_t_verbose( bool sw );
RETRIEVE_API int retrieveChDTSparameters( int rtr_d ,int ch_no ,
						int need_dtslink, int buff_size,
						unsigned int* shot_no, unsigned short* sub_shot,
						char* dts_source, char* dts_hostid, char*	dts_moduleid, char* dts_trig_ch,
						char* clk_source, char* clk_hostid, char*	clk_moduleid, char* dts_clk_ch,
						int* dts_user_define, int* time_array_size,
						char* ext_or_int, char* internal_clock, char* samp_interval, char* pre_samp,
						int64_t* last_count,	int* stat);		// 2012-07-31
// �ꎞ�I�Ώ��̂��߂̊֐��@���p�֎~
RETRIEVE_API int retrieveGetDTSInfo4matlab(const char *DiagName, unsigned int AShotNumber, unsigned int ASubShotNumber,
							  int ChannelNumber, const char *IndexServer,
			char  *DTSsource, char  *DTShostID, char  *DTSmoduleID, short *DTStriggerChannel,
			short *DTSclockChannel, int  *DTSuserDefine, int  *DTStimeArraySize,
			int *Fastch,	 int *Lastch, char *clockSource );
//----------------------------------------------------------------------------------------------- 
//
// ����ȍ~�̃��W���[���͋��o�[�W�����݊��̂��߂̃��W���[���ł���B
// ����A�g�p���Ȃ��悤�Ɍ������邱�ƁB
//
//----------------------------------------------------------------------------------------------- 
RETRIEVE_API int retrieveGetDTSLastChannel(const char *DiagName,	// �v����
					 unsigned int AShotNumber,						// �₢���킹�V���b�g�ԍ�
					 short ASubShotNumber,							// �₢���킹�T�u�V���b�g�ԍ�
					 const char *IndexServer,						// �C���f�b�N�X�T�[�o
					 unsigned int *LastChannel,
#ifdef __cplusplus
					 int TimeoutSec = 0);
#else
					 int TimeoutSec);
#endif
RETRIEVE_API int retrieveGetDTSInfo(const char *DiagName,			// �v����
					 unsigned int AShotNumber,						// �₢���킹�V���b�g�ԍ�
					 short ASubShotNumber,							// �₢���킹�T�u�V���b�g�ԍ�
					 unsigned int *ShotNumber,						// �V���b�g�ԍ�
					 unsigned short *SubShotNumber,					// �T�u�V���b�g�ԍ�
					 unsigned short *LastChannel,					// �܂܂��`���l���f�[�^��
					 short *Year,									// �N
					 short *Month,									// ��
					 short *Day,									// ��
					 short *Hour,									// ��
					 short *Minute,									// ��
					 short *Second,									// �b

					 char  *DTSsource,
					 char  *DTShostID,
					 char  *DTSmoduleID,
					 short *DTStriggerChannel,
					 short *DTSclockChannel,
					 int  *DTSuserDefine,
					 int  *DTStimeArraySize,

					 const char *IndexServer,					// �C���f�b�N�X�T�[�o
					 short *Fastch,								// �Ŏn�`���l���ԍ�
					 short *Lastch,								// �ŏI�`���l���ԍ�
					 char *clockSource);

RETRIEVE_API int retrieveGetDTSData(char *IndexServer,
							short *ch,
							float *timeArray,
							char *DTSsource,
							char *DTShostID,
							char *DTSmoduleID,
							short *DTStriggerChannel,
							short *DTSclockChannel,
							int *DTSuserDefine,
							int  DTStimeArraySize,
							unsigned int shot,
							unsigned short subshot,
							char *clockSource,
#ifdef __cplusplus
							short bSilent = 0,
							float *clockCycle = NULL,
							float *triggerTiming = NULL);
#else
							short bSilent,
							float *clockCycle,
							float *triggerTiming);
#endif
RETRIEVE_API int retrieveGetDTSDataDBL(char *IndexServer,
							short *ch,
							double *timeArray,
							char *DTSsource,
							char *DTShostID,
							char *DTSmoduleID,
							short *DTStriggerChannel,
							short *DTSclockChannel,
							int *DTSuserDefine,
							int  DTStimeArraySize,
							unsigned int shot,
							unsigned short subshot,
							char *clockSource,
#ifdef __cplusplus
							short bSilent = 0,
							double *clockCycle = NULL,
							double *triggerTiming = NULL);
#else
							short bSilent,
							double *clockCycle,
							double *triggerTiming);
#endif

RETRIEVE_API int retrieveGetDTSInfoFromParams(
						unsigned short paramterCount,
						const char *paramsList,

						unsigned int *ShotNumber,						// �V���b�g�ԍ�
						unsigned short *SubShotNumber,					// �T�u�V���b�g�ԍ�

						char  *DTSsource,
						char  *DTShostID,
						char  *DTSmoduleID,
						char  *DTStriggerChannel,

						char  *CLKsource,
						char  *CLKhostID,
						char  *CLKmoduleID,
						char  *DTSclockChannel,
						int   *DTSuserDefine,
						int   *DTStimeArraySize,

						const char *IndexServer,					// �C���f�b�N�X�T�[�o
						char *clockSource,
						char *InternalClock,
						
						char *SamplingInterval,
						char *PreSampling);

RETRIEVE_API int retrieveGetDTSData2(char *IndexServer,
							short *ch,
							float *timeArray,

							char *DTSsource,
							char *DTShostID,
							char *DTSmoduleID,
							char *DTStriggerChannel,

							char *clkSource,
							char *clkHostID,
							char *clkModuleID,
							char *DTSclockChannel,

							int *DTSuserDefine,
							int  DTStimeArraySize,
							unsigned int shot,
							unsigned short subshot,
							char *clockSource,
							char *InternalClock,

							char *strSamplingInterval,
							char *strPreSampling,

							short bSilent,
							float *clockCycle,
							float *triggerTiming);

RETRIEVE_API int retrieveGetDTSDataDBL2(char *IndexServer,
							short *ch,
							double *timeArray,

							char *DTSsource,
							char *DTShostID,
							char *DTSmoduleID,
							char *DTStriggerChannel,

							char *clkSource,
							char *clkHostID,
							char *clkModuleID,
							char *DTSclockChannel,

							int *DTSuserDefine,
							int  DTStimeArraySize,
							unsigned int shot,
							unsigned short subshot,
							char *clockSource,
							char *InternalClock,

							char *strSamplingInterval,
							char *strPreSampling,

							short bSilent,
							double *clockCycle,
							double *triggerTiming);

//RETRIEVE_API int retrieveGetParameterString(int rtr_d, int ch, int *nParam, char *strParam);

//  ----- PV-WAVE -----
RETRIEVE_API pvw_long PvwRetrieveGetDTSLastChannel(int argc,LPVOID **argv);
RETRIEVE_API pvw_long PvwRetrieveGetDTSInfo(int argc,LPVOID **argv);
//  ----- IDL -----
RETRIEVE_API int IdlRetrieveGetDTSLastChannel(int argc,LPVOID argv[]);
RETRIEVE_API int IdlRetrieveGetDTSInfo(int argc, LPVOID argv[]);
#ifdef __cplusplus
}
#endif

#define DSTR(s)	_x_str(s)
#define _x_str(s)	#s

// DTS Info
#define WRONG_PARAMETER			10 
#define WRONG_PARAMETER_NUMBERS 11 
#define SIZE_IS_SHORT			12
#define INDEX_SERVER_UNCONNECT	24
#define ERROR_FILE_OPEN			44
#define SQL_EXEC_ERROR			45

#define NOT_SUPPORTED_DTS_INFO	53
#define NO_DTS_LINK_DATA		54
#define NO_DIAG_NAME			55
#define NO_DTS_HOSTID			56
#define BEFORE_SHOT_56221		57		// RESERVED PV-WAVE ( donot have DTS LINK PARAMETER )
#define ARRAY_SIZE_ZERO			58		// RESERVED PV-WAVE
#define NO_DTSTBL_DATA			60

#define INDEX_ENV			"INDEXSERVERNAME"
/**
 *
 * @date 2013-09-25
 * @version	17.0.0
 * @note	update convert to volt for PXIe5185,PXIe5186v2
 *          BUG-FIX: �t���[���\���̃T�u�V���b�g�f�[�^���擾�ł��Ȃ��B
 *
 */
#define RETRIEVE_MAJOR		17
#define RETRIEVE_MINOR		0
#ifndef RETRIEVE_TINY
 #define RETRIEVE_TINY		0
#endif
#define RETRIEVE_VERSION	DSTR(RETRIEVE_MAJOR) "." DSTR(RETRIEVE_MINOR) "." DSTR(RETRIEVE_TINY)

#ifndef RETRIEVE_T_VERSION
 //#define RETRIEVE_T_VERSION	"16.0.3"
 #define RETRIEVE_T_VERSION		RETRIEVE_VERSION
#endif

#endif /*_RETRIEVE_H_*/
