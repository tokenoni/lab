// RTTransClientDLL.h : DLL アプリケーション用のエントリ ポイントを定義します。
//
/**
 *
 * @date	2012-08-03
 * @author	Imazu.s
 * @note	add PvwLRTGetChParam(),IdlLRTGetChParam()
 */
#ifndef RTTRANSCLIENTDLL_H
#define RTTRANSCLIENTDLL_H

#ifndef	DllClass
	#ifdef BUILD_DLL
		#define DllClass __declspec( dllexport )
	#else
		#if ( defined(WIN32)||defined(WIN64) )&&!defined(_USRSLL)
			#define DllClass __declspec( dllimport )
		#else
			#define DllClass
		#endif
	#endif // BUILD_DLL
#endif	//	DllClass

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

#if defined(WIN32) || defined(WIN64)
BOOL APIENTRY DllMain( HANDLE hModule, 
                       DWORD  ul_reason_for_call, 
                       LPVOID lpReserved);
 typedef 	__int64 int64_t;  
#else
 typedef	int 	SOCKET;
 typedef	int		BOOL;
 typedef 	void* 	LPVOID;  
#if defined(__LP64__)
 typedef	long	int64_t;
#else
 typedef	long long	int64_t;
#endif
 #define	FALSE	0
 #define	TRUE	1
 #define	INVALID_SOCKET	-1
 #define	SOCKET_ERROR	-1
 
#endif

typedef struct	lrtinf	{
	int		thinning;	//	間引き率　1/N (1～)
} LRTINF;

typedef struct	lrtparam	{
	int		samplingInterval;	//　サンプリング間隔（ナノ秒）
	int		valueLen;				//　データビット長
	int		recLen;				//　データ取得時ブロック長（Byte）
	int		packetSize;			//  パケットサイズ(Byte)
} LRTPARAM;

typedef struct	lrtcameraparam	{
	int		pixelDipth;			//　データビット長
	int		frameRate;			//　フレームレート（Frames/sec.）
	int		frameWidth;			//	フレーム幅
	int		frameHeight;			//	フレーム高さ
	int		recLen;				//　データ取得時ブロック長（Byte）
	int		packetSize;			//  パケットサイズ(Byte)
} LRTCAMERAPARAM;

struct	RT_TransClient_HANDLE {
	char	*DiagName;
	char	*ServerName;
	SOCKET	hsocket;
	SOCKET	dsocket;

	int		shotNo;
	int		subshotNo;
	int		channelNo;
	int		blockNo;

	int		thinningRate;
	int		samplingInterval;
	int		valueLen;
	int		recLen;
	int		packetSize;

	int		frameRate;
	int		frameWidth;
	int		frameHeight;

	int		tcpNo;
	int		udpNo;
	BOOL	inDataTransfer;
#ifdef _cplusplus
 #if defined(WIN32) || defined(WIN64)
	RT_TransClient_HANDLE() {DiagName=ServerName=NULL;hsocket=dsocket=NULL;inDataTransfer=FALSE;};
 #else
	RT_TransClient_HANDLE() {DiagName=ServerName=NULL;hsocket=dsocket=INVALID_SOCKET;inDataTransfer=FALSE;};
 #endif
#endif
} ;


/***************************************************
	ローカル関数群
	関数名の先頭にlc(local)がついたものは
	DLL内でのみ使用されるローカル関数
****************************************************/
#ifdef __cplusplus
#if defined(BUILD_DLL) || !( defined(WIN32) || defined(WIN64) )
u_long	lcLookupAddress(const char *servername);
SOCKET	lcEstablishConnection(u_long nRemoteAddr, u_short nPort, bool inTCP=TRUE);
int		lcGetRTTransServerName(char *diagname, char *servername);

int		lcOpenRTTransd(char *servername, int chno, int tcpno, int udpno, int *packetsize, SOCKET *sd);
int		lcCloseRTTransd(SOCKET *sd, int channel);
int		lcCloseRTTransd(RT_TransClient_HANDLE *hnd);

int		lcSetInfRTTransd(RT_TransClient_HANDLE *hnd, LRTINF inf);
int		lcGetChannelParamRTTransd(RT_TransClient_HANDLE *hnd, int *status, LRTPARAM *prm);
int		lcCameraGetChannelParamRTTransd(RT_TransClient_HANDLE *hnd, int *status, LRTCAMERAPARAM *prm);

int		lcStartTrans(RT_TransClient_HANDLE *hnd);
int		lcEndTrans(RT_TransClient_HANDLE *hnd);

int		lcGetRTTransdServicePort(void);
int		lcGetRTTransdDataPort(void);
int		lcGetRTTransdDataPortAndSocket(int *udpport, SOCKET *udpsock);
int		lcSendCommand(SOCKET sd, const char *str);
int		lcReceiveReply(SOCKET sd, char *str);

int		lcGetChannelData(RT_TransClient_HANDLE *hnd, int *status,
							int *shotno, int *subshot,int *blkno, int64_t *head, int64_t *tail,
							char *buf, int *len);
int		lcGetChParamsNum(RT_TransClient_HANDLE *hnd, int *paramsNum);
int		lcGetChAllParams(RT_TransClient_HANDLE *hnd, char **Key, char **Value);
int		lcGetChParam(RT_TransClient_HANDLE *hnd, char *Key, char *Value);
int		lcConvertYUVtoRGB(char *YUVdata, char *RGBdata, short int size_x, short int size_y);

int		lcGetParamString(char *src, char *dst);

#endif
#endif
/***************************************************
	DLL呼び出し関数群
****************************************************/
#ifdef __cplusplus
extern "C"
{
#endif
DllClass	pvw_long	PvwLRTStart(int argc, LPVOID **argv);
DllClass	pvw_long	PvwLRTEnd(int argc, LPVOID **argv);
DllClass	pvw_long	PvwLRTOpen(int argc, LPVOID **argv);
DllClass	pvw_long	PvwLRTClose(int argc, LPVOID **argv);
DllClass	pvw_long	PvwLRTSetInf(int argc, LPVOID **argv);
DllClass	pvw_long	PvwLRTGetChannelParam(int argc, LPVOID **argv);
DllClass	pvw_long	PvwLRTGetChannelData(int argc, LPVOID **argv);

DllClass	pvw_long	PvwLRTCameraStart(int argc, LPVOID **argv);
DllClass	pvw_long	PvwLRTCameraEnd(int argc, LPVOID **argv);
DllClass	pvw_long	PvwLRTCameraOpen(int argc, LPVOID **argv);
DllClass	pvw_long	PvwLRTCameraClose(int argc, LPVOID **argv);
DllClass	pvw_long	PvwLRTCameraSetInf(int argc, LPVOID **argv);
DllClass	pvw_long	PvwLRTCameraGetChannelParam(int argc, LPVOID **argv);
DllClass	pvw_long	PvwLRTCameraGetChannelData(int argc, LPVOID **argv);
DllClass	pvw_long	PvwLRTYUVCameraGetChannelData(int argc, LPVOID **argv);

DllClass	pvw_long	PvwLRTGetChParam(int argc, LPVOID **argv);

DllClass	int		IdlLRTStart(int argc, LPVOID argv[]);
DllClass	int		IdlLRTEnd(int argc, LPVOID argv[]);
DllClass	int		IdlLRTOpen(int argc, LPVOID argv[]);
DllClass	int		IdlLRTClose(int argc, LPVOID argv[]);
DllClass	int		IdlLRTSetInf(int argc, LPVOID argv[]);
DllClass	int		IdlLRTGetChannelParam(int argc, LPVOID argv[]);
DllClass	int		IdlLRTGetChannelData(int argc, LPVOID argv[]);

DllClass	int		IdlLRTCameraStart(int argc, LPVOID argv[]);
DllClass	int		IdlLRTCameraEnd(int argc, LPVOID argv[]);
DllClass	int		IdlLRTCameraOpen(int argc, LPVOID argv[]);
DllClass	int		IdlLRTCameraClose(int argc, LPVOID argv[]);
DllClass	int		IdlLRTCameraSetInf(int argc, LPVOID argv[]);
DllClass	int		IdlLRTCameraGetChannelParam(int argc, LPVOID argv[]);
DllClass	int		IdlLRTCameraGetChannelData(int argc, LPVOID argv[]);
DllClass	int		IdlLRTYUVCameraGetChannelData(int argc, LPVOID argv[]);

DllClass	int		IdlLRTGetChParam(int argc, LPVOID argv[]);

DllClass	int		LRTStart();
DllClass	int		LRTEnd();

#ifdef __cplusplus
 DllClass	int		LRTOpen(char *diagname, int chno, LRTINF inf, void *hnd, int udpp = 0, int packetsize=0);
#else
 DllClass	int		LRTOpen(char *diagname, int chno, LRTINF inf, void *hnd, int udpp, int packetsize);
#endif
DllClass	int		LRTClose(void *hnd);

DllClass	int		LRTSetInf(void *hnd,  LRTINF inf);
DllClass	int		LRTGetChannelParam(void *hnd, int *status, LRTPARAM *prm);
DllClass	int		LRTGetChannelParam2(void *hnd, int *status,
						int* samp_interval, int* value_len, int* rec_len,int* packet_size);

DllClass	int		LRTGetChannelData(void *hnd, int *status,
									  int *shotno, int *subshot,int *blkno, int64_t *head, int64_t *tail,
									  void *buf, int *len);
DllClass	int		LRTCameraStart();
DllClass	int		LRTCameraEnd();

#ifdef __cplusplus
 DllClass	int		LRTCameraOpen(char *diagname, int chno, LRTINF inf, void *hnd, int udpp = 0, int packetsize=0);
#else
 DllClass	int		LRTCameraOpen(char *diagname, int chno, LRTINF inf, void *hnd, int udpp, int packetsize);
#endif
DllClass	int		LRTCameraClose(void *hnd);

DllClass	int		LRTCameraSetInf(void *hnd,  LRTINF inf);
DllClass	int		LRTCameraGetChannelParam(void *hnd, int *status, LRTCAMERAPARAM *prm);
DllClass	int		LRTCameraGetChannelParam2(void *hnd, int *status,
		int* pixelDipth, int* frameRate, int* frameWidth, int* frameHeight, int* recLen,int* packetSize);

DllClass	int		LRTCameraGetChannelData(void *hnd, int *status,
									  int *shotno, int *subshot,int *blkno, int64_t *head, int64_t *tail,
									  void *buf, int *len);
DllClass	int		LRTYUVCameraGetChannelData(void *hnd, int *status,
									  int *shotno, int *subshot,int *blkno, int64_t *head, int64_t *tail,
									  void *buf, int *len);

DllClass	int		LRTGetChParamsNum(void *hnd, int *paramsNum);
DllClass	int		LRTGetChAllParams(void *hnd, char **Key, char **Value);
DllClass	int		LRTGetChParam(void *hnd, char *Key, char *Value);
DllClass	void	LRTset_verbose( bool sw );
DllClass	int		LRTget_last_error( void );
DllClass	void 		LRTbgrx2rgb4matlab(void *dst, void* src,int width, int height);

DllClass	int		LRTOpen2(char *diagname, int chno, LRTINF inf, void *hnd, int udpp, int packetsize, const char* host);

#ifdef __cplusplus
}
#endif

/***************************************************
 Error codes:
 The following error codes are peculiar to this DLL.  
****************************************************/
#define RT_NO_ERROR							0
#define RT_NO_DIAG_NAME						1
#define RT_CHANNEL_OVER_RANGE				2
#define RT_FAILED_CONNECTION_INDEXSERVER	3
#define RT_NOT_DEFINED_INDEXSERVERNAME	4
#define RT_ERROR_IN_SQL_EXEC				5
#define RT_CANT_FIND_SERVER_ADDRESS		6
#define RT_CANT_CONNECT_TO_SERVER			7
#define RT_ERROR_IN_SENDING_COMMAND		8
#define RT_ERROR_IN_RECEIVING_REPLY		9
#define RT_ALREADY_CLOSED_HANDLE			10
#define RT_ERROR_IN_RECEIVING_DATA			11
#define RT_INVALID_SOCKET					12
#define RT_ERROR_IN_SET_TIMEOUT			13
#define RT_TIMEOUT							14
#define RT_WRONG_PARAMETER_NUMBERS			15
#define RT_NOT_WAVE_DATA					16
#define RT_NOT_IMAGE_DATA					17
#define RT_ILLEAGAL_HANDLE					18

#define RT_REPLY_OK							100
#define RT_REPLY_GETCH_PARAM_OK			101
#define RT_REPLY_END_TRANS_OK				102
#define RT_REPLY_GETCH_PARAMSNUM_OK		103
#define RT_REPLY_GETCH_ALLPARAMS_OK		104
#define RT_REPLY_GETCH_PARAMBYKEY_OK		105

#define INDEX_SERVER_UNCONNECT				24
#define NO_DIAG_NAME							55

#endif	//	RTTRANSCLIENTDLL_H
