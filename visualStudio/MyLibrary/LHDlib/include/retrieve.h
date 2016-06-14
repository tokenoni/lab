#ifndef _RETRIEVE_H_
#define _RETRIEVE_H_

#if defined(WIN32) || defined(WIN64)
 #ifndef _GNUC__ 
  #pragma once

// 下で指定された定義の前に対象プラットフォームを指定しなければならない場合、以下の定義を変更してください。
// 異なるプラットフォームに対応する値に関する最新情報については、MSDN を参照してください。
  #ifndef WINVER				// Windows XP 以降のバージョンに固有の機能の使用を許可します。
   #define WINVER 0x0501		// これを Windows の他のバージョン向けに適切な値に変更してください。
  #endif

  #ifndef _WIN32_WINNT		// Windows XP 以降のバージョンに固有の機能の使用を許可します。                   
   #define _WIN32_WINNT 0x0501	// これを Windows の他のバージョン向けに適切な値に変更してください。
  #endif						

  #ifndef _WIN32_WINDOWS		// Windows 98 以降のバージョンに固有の機能の使用を許可します。
   #define _WIN32_WINDOWS 0x0410 // これを Windows Me またはそれ以降のバージョン向けに適切な値に変更してください。
  #endif

  #ifndef _WIN32_IE			// IE 6.0 以降のバージョンに固有の機能の使用を許可します。
   #define _WIN32_IE 0x0600	// これを IE. の他のバージョン向けに適切な値に変更してください。
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
#else
 typedef void* LPVOID;  
 #define RETRIEVE_API

#endif

#ifndef INDEX_ENV_UNDEFINED
 #define INDEX_ENV_UNDEFINED -101
#endif

#ifndef RETRIEVE_T_VERSION
 #define RETRIEVE_T_VERSION	"15.1.0"
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
 typedef __int64		pvw_long;
 typedef __int64		pvw_ulong;
#elif defined(__LP64__)
 typedef long			pvw_long;
 typedef unsigned long	pvw_ulong;
#else
 typedef int			pvw_long;
 typedef unsigned int	pvw_ulong;
#endif
 typedef unsigned short  pvw_uint;
 typedef short  pvw_int;

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
RETRIEVE_API pvw_long PvwRetrieveGetDTSInfo(int argc,LPVOID **argv);
RETRIEVE_API pvw_long PvwRetrieveGetDTSInfoFromParams(int argc,LPVOID **argv);	//	add 2011-03-23
RETRIEVE_API pvw_long PvwRetrieveGetDTSLastChannel(int argc,LPVOID **argv);

RETRIEVE_API pvw_long PvwRetrieveGetDTSParameters(int argc, LPVOID **argv);	// add 2010-09-01
RETRIEVE_API pvw_long PvwRetrieveGetDTSParametersCount(int argc, LPVOID **argv);	// add 2010-09-01

RETRIEVE_API pvw_long PvwRetrieveRealShotNumber(int argc, LPVOID **argv);	// add 2009-09-14
RETRIEVE_API pvw_long PvwRetrieveChV(int argc, LPVOID **argv);				// add 2010-03-05
RETRIEVE_API pvw_long PvwRetrieveChVolts(int argc, LPVOID **argv);			// add 2010-03-05
RETRIEVE_API pvw_long PvwRetrieveChVoltsDbl(int argc, LPVOID **argv);		// add 2010-03-05

#if defined(WIN64) || defined(__LP64__)
RETRIEVE_API pvw_long PvwLong64From32(int argc,LPVOID **argv);
void pvw_long64From32( pvw_long* dst,int* src,unsigned int len);
#endif
RETRIEVE_API pvw_long PvwRGBfromYUV422(int argc,LPVOID **argv);
RETRIEVE_API pvw_long PvwIntgerFromByte(int argc, LPVOID **argv);

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
RETRIEVE_API int IdlRetrieveGetDTSInfo(int argc, LPVOID argv[]);
RETRIEVE_API int IdlRetrieveGetDTSInfoFromParams(int argc, LPVOID argv[]);	//	add 2011-05-12
RETRIEVE_API int IdlRetrieveGetDTSLastChannel(int argc,LPVOID argv[]);

RETRIEVE_API int IdlRetrieveGetDTSParameters(int argc, LPVOID **argv);	// add 2010-10-01
RETRIEVE_API int IdlRetrieveGetDTSParametersCount(int argc, LPVOID **argv);	// add 2010-10-01

RETRIEVE_API int IdlRetrieveRealShotNumber( int argc, LPVOID argv[] );	// add 2009-09-14
RETRIEVE_API int IdlRetrieveChV(int argc, LPVOID argv[]);				// add 2010-03-05
RETRIEVE_API int IdlRetrieveChVolts(int argc, LPVOID argv[]);			// add 2010-03-05
RETRIEVE_API int IdlRetrieveChVoltsDbl(int argc, LPVOID argv[]);		// add 2010-03-05
RETRIEVE_API int IdlIntgerFromByte(int argc, LPVOID argv[]);
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
					(const char *DiagName,							// 計測名
					 unsigned int AShotNumber,						// 問い合わせショット番号
					 short ASubShotNumber,							// 問い合わせサブショット番号
					 unsigned int *ShotNumber,						// ショット番号
					 unsigned short *SubShotNumber,					// サブショット番号
					 unsigned short *LastChannel,					// 含まれるチャネルデータ数
					 short *Year,									// 年
					 short *Month,									// 月
					 short *Day,									// 日
					 short *Hour,									// 時
					 short *Minute,									// 分
					 short *Second,									// 秒

					 char  *DTSsource,
					 char  *DTShostID,
					 char  *DTSmoduleID,
					 short *DTStriggerChannel,
					 short *DTSclockChannel,
					 int  *DTSuserDefine,
					 int  *DTStimeArraySize,

					 const char *IndexServer,					// インデックスサーバ
					 short *Fastch,								// 最始チャネル番号
					 short *Lastch,								// 最終チャネル番号
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


RETRIEVE_API bool retrieve_check_string(const char* str);

RETRIEVE_API int retrieveGetParameterString(int rtr_d, int ch, int *nParam, char *strParam);

RETRIEVE_API int retrieveGetDTSLastChannel(const char *DiagName,	// 計測名
					 unsigned int AShotNumber,						// 問い合わせショット番号
					 short ASubShotNumber,							// 問い合わせサブショット番号
					 const char *IndexServer,						// インデックスサーバ
					 unsigned int *LastChannel,
					 int TimeoutSec = 0);

RETRIEVE_API int retrieveGetDTSInfo(const char *DiagName,			// 計測名
					 unsigned int AShotNumber,						// 問い合わせショット番号
					 short ASubShotNumber,							// 問い合わせサブショット番号
					 unsigned int *ShotNumber,						// ショット番号
					 unsigned short *SubShotNumber,					// サブショット番号
					 unsigned short *LastChannel,					// 含まれるチャネルデータ数
					 short *Year,									// 年
					 short *Month,									// 月
					 short *Day,									// 日
					 short *Hour,									// 時
					 short *Minute,									// 分
					 short *Second,									// 秒

					 char  *DTSsource,
					 char  *DTShostID,
					 char  *DTSmoduleID,
					 short *DTStriggerChannel,
					 short *DTSclockChannel,
					 int  *DTSuserDefine,
					 int  *DTStimeArraySize,

					 const char *IndexServer,					// インデックスサーバ
					 short *Fastch,								// 最始チャネル番号
					 short *Lastch,								// 最終チャネル番号
					 char *clockSource);

RETRIEVE_API int retrieveGetDTSInfoFromParams(
						unsigned short paramterCount,
						const char *paramsList,

						unsigned int *ShotNumber,						// ショット番号
						unsigned short *SubShotNumber,					// サブショット番号

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

						const char *IndexServer,					// インデックスサーバ
						char *clockSource,
						char *InternalClock,
						
						char *SamplingInterval,
						char *PreSampling);

RETRIEVE_API int getParamValueFromList(
						const char *paramList,
						const char *pName,
						char *pValue);

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
RETRIEVE_API int retrieve_t_EXE(int argc, char **argv);
RETRIEVE_API void bufferStringCopy(char* string,void* buffer,int element_size,int number);
RETRIEVE_API void set_retrieve_t_verbose( bool sw );
// 一時的対処のための関数　流用禁止
RETRIEVE_API int retrieveGetDTSInfo4matlab(const char *DiagName, unsigned int AShotNumber, unsigned int ASubShotNumber,
							  int ChannelNumber, const char *IndexServer,
			char  *DTSsource, char  *DTShostID, char  *DTSmoduleID, short *DTStriggerChannel,
			short *DTSclockChannel, int  *DTSuserDefine, int  *DTStimeArraySize,
			int *Fastch,	 int *Lastch, char *clockSource );

#ifdef __cplusplus
}
#endif

// DTS Info
#define WRONG_PARAMETER			10 
#define WRONG_PARAMETER_NUMBERS 11 
#define INDEX_SERVER_UNCONNECT	24
#define SQL_EXEC_ERROR			45

#define NOT_SUPPORTED_DTS_INFO	53
#define NO_DTS_LINK_DATA		54
#define NO_DIAG_NAME			55
#define NO_DTS_HOSTID			56

#define INDEX_ENV			"INDEXSERVERNAME"
#define RETRIEVE_VERSION	"15.1.0"

#endif /*_RETRIEVE_H_*/
