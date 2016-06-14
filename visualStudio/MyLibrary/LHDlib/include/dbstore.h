#ifndef DBSTORE_H_
#define DBSTORE_H_

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
#ifdef DBSTORE_EXPORTS
#define DBSTORE_API __declspec(dllexport)
#else
#define DBSTORE_API __declspec(dllimport)
#endif
#else
#define DBSTORE_API 
#endif

#define fordbsopen_( x1,x2,x3,x4,x5 )			__stdcall FORDBSOPEN( x1,x2,x3,x4,x5 )
#define fordbswrite_( x1,x2,x3,x4,x5,x6,x7,x8 )	__stdcall FORDBSWRITE( x1,x2,x3,x4,x5,x6,x7,x8 )
#define fordbsopenframe_( x1,x2,x3,x4,x5 )		__stdcall FORDBSOPENFRAME( x1,x2,x3,x4,x5 )
#define fordbswriteframe_( x1,x2,x3,x4,x5 )		__stdcall FORDBSWRITEFRAME( x1,x2,x3,x4,x5 )
#define fordbswriteframes_( x1,x2,x3,x4,x5,x6 )	__stdcall FORDBSWRITEFRAMES( x1,x2,x3,x4,x5,x6 )
#define fordbscloseframe_( x1,x2,x3,x4,x5 )		__stdcall FORDBSCLOSEFRAME( x1,x2,x3,x4,x5 )
#define fordbsclose_( x1 )						__stdcall FORDBSCLOSE( x1 )
#define fordbsabort_( x1 )						__stdcall FORDBSABORT( x1 )
#define fordbsversion_( )						__stdcall FORDBSVERSION( )
#define fordbserrormessage_( x1 )				__stdcall FORDBSERRORMESSAGE( x1 )

#else
typedef void* LPVOID;  
#define DBSTORE_API

#endif

#if defined(WIN64)
 typedef __int64		pvw_long;
 typedef __int64		pvw_ulong;
 #define NOT_SUPPORTED_PV_INT	"INT32"	
#elif defined(__LP64__)
 typedef long			pvw_long;
 typedef unsigned long	pvw_ulong;
 #define NOT_SUPPORTED_PV_INT	"INT32"	
#else
 typedef int			pvw_long;
 typedef unsigned int	pvw_ulong;
 #define NOT_SUPPORTED_PV_INT	"INT64"	
#endif
 typedef unsigned short  pvw_uint;
 typedef short  pvw_int;

#ifdef __cplusplus
 #include <vector>
#endif

 #ifdef __cplusplus
 extern "C" {
#endif
//  ----- PV-WAVE -----
DBSTORE_API pvw_long PvwdbsOpen(int argc, LPVOID **argv);
DBSTORE_API pvw_long PvwdbsWrite(int argc, LPVOID **argv);
DBSTORE_API pvw_long PvwdbsClose(int argc, LPVOID **argv);
DBSTORE_API pvw_long PvwdbsAbort(int argc, LPVOID **argv);
DBSTORE_API pvw_long PvwdbsOpenFrame(int argc, LPVOID **argv);
DBSTORE_API pvw_long PvwdbsWriteFrame(int argc, LPVOID **argv);
DBSTORE_API pvw_long PvwdbsWriteFrames(int argc, LPVOID **argv);
DBSTORE_API pvw_long PvwdbsCloseFrame(int argc, LPVOID **argv);
DBSTORE_API const char* PvwdbsVersion(int argc, LPVOID **argv);
DBSTORE_API const char* PvwdbsErrorMessage(int argc, LPVOID **argv);

//  ----- IDL -----
DBSTORE_API long IdldbsOpen(int argc, LPVOID argv[]);
DBSTORE_API long IdldbsWrite(int argc, LPVOID argv[]);
DBSTORE_API long IdldbsClose(int argc, LPVOID argv[]);
DBSTORE_API long IdldbsAbort(int argc, LPVOID argv[]);
DBSTORE_API long IdldbsOpenFrame(int argc, LPVOID argv[]);
DBSTORE_API long IdldbsWriteFrame(int argc, LPVOID argv[]);
DBSTORE_API long IdldbsWriteFrames(int argc, LPVOID argv[]);
DBSTORE_API long IdldbsCloseFrame(int argc, LPVOID argv[]);
DBSTORE_API const char* IdldbsVersion(int argc, LPVOID argv[]);
DBSTORE_API const char* IdldbsErrorMessage(int argc, LPVOID argv[]);


//  ----- Fortrun -----
DBSTORE_API long fordbsopen_(char* MailAddress,char *DiagName,
			unsigned long  AShotNumber, unsigned short ASubShotNumber, short DataType);
			
DBSTORE_API long fordbswrite_(long dbs, long ch_no,
			unsigned long  ParameterC,	 unsigned long  ParaLength,unsigned char *ParameterList,
			unsigned long DataSize, unsigned char *data, char *data_type );

DBSTORE_API long fordbsopenframe_(long dbs, long ch_no,
			unsigned long  x_size,	 unsigned long  y_size, char *data_type );
			
DBSTORE_API long fordbswriteframe_(long dbs, long ch_no, long frame_no,
			unsigned long DataSize, unsigned char *data);

DBSTORE_API long fordbswriteframes_(long dbs, long ch_no, long start_frame_no, long num_frame,
			unsigned long DataSize, unsigned char *data);

DBSTORE_API long fordbscloseframe_(long dbs, long ch_no,
			unsigned long ParameterC, unsigned long  ParaLength, unsigned char *ParameterList);

DBSTORE_API long fordbsclose_(long dbs);
DBSTORE_API long fordbsabort_(long dbs);
DBSTORE_API const char* fordbsversion_( );
DBSTORE_API const char* fordbserrormessage_(long error_code);


//  ----- C/C++ -----
DBSTORE_API int dbsOpen(const char *mail_address, const char *diag_name, 
			 unsigned int shot_number, unsigned short sub_shot, short data_type );

DBSTORE_API int dbsClose(int dbs_d);
DBSTORE_API int dbsAbort(int dbs_d);

DBSTORE_API int dbsWrite(int dbs_d, int ch_no,
			unsigned int param_cnt,unsigned int param_size, unsigned char *param_list,
			unsigned int data_length, unsigned char* data, char* data_type );

DBSTORE_API int dbsOpenFrame(int dbs_d, int ch_no,
			unsigned int x_size,unsigned int y_size, char* data_type );
			
DBSTORE_API int dbsWriteFrame(int dbs_d, int ch_no, int frame_no,
			unsigned int data_length, void* data);
			
DBSTORE_API int dbsWriteFrames(int dbs_d, int ch_no, int start_frame_no,int num_frame,
			unsigned int data_length, void* data);

DBSTORE_API int dbsCloseFrame(int dbs_d, int ch_no,
			unsigned int param_cnt,unsigned int param_size, unsigned char *param_list);

DBSTORE_API int dbsCloseFrame2(int dbs_d, int ch_no, 
		unsigned int plist_serial_bytes, void *plist_serial);	// add 2010-11-19 for MATLAB

DBSTORE_API unsigned int dbsSetParam2Buffer( void* buffer, unsigned int offset,
		const char* pname, const char* pvalue, const char* ptype); // add 2010-11-19 for MATLAB

DBSTORE_API const char* dbsVersion( void );
DBSTORE_API const char* dbsZversion( void );	// add 13.2.1 2010-04-02

DBSTORE_API const char* dbsErrorMessage( int error_code );

DBSTORE_API int dbsWrite2(int dbs_d, int ch_no,
			unsigned int param_cnt, char* param_item[], char* param_val[], int param_type[],
			unsigned int data_length, unsigned char* data, char* data_type );
DBSTORE_API int dbsWrite3(int dbs_d, int ch_no,
			unsigned int param_cnt, char* param_item[], char* param_val[], char* param_type[],
			unsigned int data_length, unsigned char* data, char* data_type );

DBSTORE_API int dbsWrite5(int dbs_d, int ch_no,
			unsigned int plist_serial_bytes, void *plist_serial,
			unsigned int data_length, void* data, char* data_type );

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
#if !defined(_USE_CH_PARAM_STRUCTURE)
typedef struct	channel_param_structure
{
	char* paramName;	// NULL terminate
	char* paramValue;	// NULL terminate
	int paramType;
	int	nameBuffSize;	// NULL including
	int	valueBuffSize;	// NULL including
	
}	CH_PARAM_STRUCTURE;
#define _USE_CH_PARAM_STRUCTURE
#endif

DBSTORE_API int dbsWrite4(int dbs_d, int ch_no,
			std::vector<CH_PARAM_STRUCTURE*>& params,
			unsigned int data_length, unsigned char* data, char* data_type );
#endif

#define INDEX_ENV		"INDEXSERVERNAME"
/**
 *
 * @version	14.0.0
 * @note	FTP connect timeout to sever (10sec) FTA_CONNECT_TIMEOUT
 * @version	14.1.0
 * @note	Supported MATLAB
 */
#define DBSTORE_VERSION	"15.1.0"
/***************************************************
 parameter type
****************************************************/
#define TYPE_UNDEF	0
#define TYPE_STRING	1
#define TYPE_BYTE	2
#define TYPE_SHORT	3
#define TYPE_INT	4
#define TYPE_FLOAT	5
#define TYPE_DOUBLE	6
#define TYPE_INT64	7			// add 13.2.1		2010-04-01

/***************************************************
 Error codes:
 The following error codes are peculiar to this DLL.
****************************************************/
#define ILLEGAL_DATA_TYPE			-1
#define ILLEGAL_MAIL_ADDRESS		-2
#define ILLEGAL_DIAG_NAME			-3
#define ILLEGAL_DBS_DESCRIPTOR		-4
#define ILLEGAL_PARAM_COUNT			-5
#define OVER_DBS_DESCRIPTOR			-6
#define ERROR_LOCAL_FILE_OPEN		-7
#define ALREADY_WRITE_CHANNEL		-8
#define ALREADY_WRITE_FRAME			-9	
#define APPLICATION_BUG				-10	
#define WRONG_PARAMETER_NUMBERS	-11 
#define ERROR_ZLIB_OCCURED			-12		 
#define ERROR_ADD_SYS_PARAM			-13		 
#define ERROR_ADD_USER_PARAM		-14	 
#define DUPLICATE_PARAMETER_NAME	-15	
#define ERROR_ZIP_CH_WRITE			-16	
#define ERROR_ZIP_FRAME_WRITE		-17	
#define ERROR_ZIP_SHOT_WRITE		-18	
#define ERROR_YET_OPEN_FRAME		-19	
#define ALREADY_CLOSED_FRAME		-20	
#define ERROR_ZLIB_OR_JPEG_LS		-21	 
#define ILLEGAL_FRAME_DATA_SIZE	-23		 

#define INDEX_ENV_UNDEFINED			-22
#define INDEX_SERVER_UNCONNECT		-24
#define ERROR_NOT_FRAME_TYPE		-25	
#define WRONG_PARAMETER_ITEM		-26	
#define ILLEGAL_IMAGE_TYPE			-27	
#define SQL_EXEC_ERROR				-45
#define MEDIA_ID_NOT_FOUND			-50	
#define ROOT_PATH_ID_NOT_FOUND		-51	
#define REGISTER_ROOT_NOT_FOUND	-52		
#define REGISTER_HOST_NOT_FOUND	-53
#define CYCLE_NUMBER_NOT_FOUND		-54
#define DIAG_NAME_NOT_FOUND			-55
#define SHOT_DATA_FOUND				-56
#define WRONG_PARAMETER_TYPE		-57
#define SITE_NAME_NOT_FOUND			-58
#define ERROR_NOT_EXIST_FRAME		-78	
#define ERROR_NOT_EXIST_CHANNEL	-79		
#define ERROR_FTP_000				-80	
#define ERROR_FTP_020				-100

#define INDEX_RECORD_ADD_FAILURE_0	-101
#define INDEX_RECORD_ADD_FAILURE_1	-102
#define INDEX_RECORD_ADD_FAILURE_2	-103
#define INDEX_RECORD_ADD_FAILURE_3	-104
#define INDEX_RECORD_ADD_FAILURE_4	-105
#define INDEX_RECORD_ADD_FAILURE_5	-106
#define FILE_NOT_FOUND				-67
#define DATA_FILE_OPEN_ERROR		-68
#define DATA_FILE_READ_ERROR		-69
#define PARAM_FILE_OPEN_ERROR		-70
#define PARAM_FILE_READ_ERROR		-71

#define ILLEGAL_PARAMLIST_FORMAT	-201
#define OVER_PARAMLIST_SIZE			-202
#define ILLEGAL_VALUE_TYPE			-203
#define NOT_DATA_AND_PARAM			-204
#define NOT_FRAMESIZE_PARAMETER		-205
#define ILLEGAL_ARGUMENT			-999
#define ILLEGAL_SHARED_LIBRARY		-9998
#define ERROR_CODE_END_NUMBER		-9999

#endif /*DBSTORE_H_*/

